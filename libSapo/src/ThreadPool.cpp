#include "ThreadPool.h"

/**
 * @brief Get information about a batch
 *
 * @param[in] batch_id is the searched batch id
 * @return a reference to the batch information
 */
std::shared_ptr<ThreadPool::BatchInfo>
ThreadPool::get_batch_info(const ThreadPool::BatchId batch_id)
{
  // search for the batch id in the map
  auto found = _batches.find(batch_id);

  // if not found throw an exception
  if (found == std::end(_batches)) {
    throw std::runtime_error("Unknown batch");
  }

  // otherwise return a reference to the batch info
  return found->second;
}

/**
 * @brief Extract next task from the queue
 *
 * @return The next task to be executed
 */
bool ThreadPool::extract_next_task(ThreadPool::Task &next,
                                   std::unique_lock<std::mutex> &lock)
{
  bool owns_lock = lock.owns_lock();

  if (!owns_lock) {
    lock.lock();
  }

  // if the queue is empty
  while (_queue.empty()) {
    // otherwise wait for new task in the queue
    _waiting_task.wait(lock);

    // if the pool is about to be destroyed
    if (_terminating) {

      if (!owns_lock) {
        lock.unlock();
      }

      // return a fake task
      return false;
    }
  }

  std::swap(next, _queue.front());
  _queue.pop();

  // update the info about the number of tasks
  // of the batch in the queue
  --(get_batch_info(next.second)->in_queue);

  if (!owns_lock) {
    lock.unlock();
  }
  return true;
}

/**
 * @brief The main thread loop
 *
 * @param[in] thread_id is the thread id in the pool
 */
void ThreadPool::consumer_loop(unsigned int thread_id)
{
  (void)thread_id;
  Task task;

  std::unique_lock<std::mutex> lock(_mutex);

  // forever
  while (true) {
    // extract the next task from the queue
    extract_next_task(task, lock);

    // if the pool is about to be destroyed
    // exit from the loop
    if (_terminating) {
      lock.unlock();

      _waiting_task.notify_all();
      return;
    }

    auto info = get_batch_info(task.second);

    // update the info about the number of
    // running tasks in the batch
    ++(info->running);

    // prepare for the task run and
    // exit from the critical region
    lock.unlock();

    // Run the task
    task.first();

    // re-enter the critical region
    lock.lock();

    // update the running tasks in the batch and
    // check for task end
    if (--(info->running) == 0 && info->in_queue == 0) {
      lock.unlock();

      // notify the waiting threads that
      // this batch has been concluded
      info->waiting_end.notify_all();

      lock.lock();
    }
  }
}

/**
 * @brief Create a new Thread Pool object
 *
 * @param[in] num_of_threads is the number of thread in the pool
 */
ThreadPool::ThreadPool(const unsigned num_of_threads):
    _threads(), _queue(), _terminating(false)
{
  reinit(num_of_threads);
}

/**
 * @brief Create a pull of threads according to the hardware
 */
ThreadPool::ThreadPool(): ThreadPool(std::thread::hardware_concurrency()) {}

/**
 * @brief Initialize a new task batch
 *
 * @return The batch id of the new batch
 */
ThreadPool::BatchId ThreadPool::create_batch()
{
  std::unique_lock<std::mutex> lock(_mutex);

  BatchId batch_id;

  // search for an unused batch id
  auto first_unused = std::begin(_unused_batch_ids);

  // if any of them is available
  if (first_unused != std::end(_unused_batch_ids)) {

    // acquire it
    batch_id = *first_unused;

    // remove it from the set of the unused id
    _unused_batch_ids.erase(first_unused);
  } else { // if no unused id is available
    // get the last used id
    auto last_batch = _batches.rbegin();

    // and increase it by one
    batch_id = last_batch->first + 1;
  }

  // build an empty info for the new batch id
  _batches.emplace(batch_id, std::make_shared<BatchInfo>());

  return batch_id;
}

/**
 * @brief Close the batch and remove it from the pool
 *
 * @param[in] batch_id is the id of the batch to be closed
 */
void ThreadPool::close_batch(const unsigned int batch_id)
{
  std::unique_lock<std::mutex> lock(_mutex);

  // batch id 0 is reserved and cannot be closed
  if (batch_id == 0) {
    return;
  }

  // get batch info
  auto info = get_batch_info(batch_id);

  // if the batch has not finished yet
  while (info->unfinished_tasks() > 0) {
    // wait for its end
    info->waiting_end.wait(lock);
  }

  _batches.erase(batch_id);

  // if the pool only contains info about
  // batch 0
  if (_batches.size() == 1) {
    // clear the set of unused id
    _unused_batch_ids.clear();
  } else {
    // add the batch id to the set
    // of the unused id
    _unused_batch_ids.insert(batch_id);
  }
}

/**
 * @brief Join the thread pool and complete the batch
 *
 * @param[in] batch_id is the id of the batch that must be completed
 */
void ThreadPool::join_threads(const unsigned int batch_id)
{
  Task task;
  std::unique_lock<std::mutex> lock(_mutex);

  // repeat forever
  while (true) {
    // get batch info
    auto info = get_batch_info(batch_id);

    // if finished
    if (info->unfinished_tasks() == 0) {

      // end loop and exit
      return;
    }

    // if the queue is not empty
    if (_queue.empty()) {
      // some tasks in the batch are still running
      // wait for them

      while (info->unfinished_tasks() > 0) {
        info->waiting_end.wait(lock);
      }

      // exit
      return;
    }

    // extract a task from the queue
    std::swap(task, _queue.front());
    _queue.pop();

    // get info about the extracted task batch
    auto ex_info = get_batch_info(task.second);

    // update it
    --(ex_info->in_queue);
    ++(ex_info->running);

    // prepare for the task run and
    // exit from the critical region
    lock.unlock();

    task.first();

    // re-enter the critical region
    lock.lock();

    // decrease the number of running
    // tasks in the batch
    --(ex_info->running);

    if (ex_info->unfinished_tasks() == 0) {
      lock.unlock();

      ex_info->waiting_end.notify_all();

      lock.lock();
    }
  }
}

/**
 * @brief Terminate the pool
 */
void ThreadPool::terminate()
{
  std::unique_lock<std::mutex> lock(_mutex);

  // set the termination flag
  _terminating = true;
  while (!_queue.empty()) {
    _queue.pop();
  }

  lock.unlock();

  // notify the end of all the active batches
  // to all the main threads waiting for it
  for (auto &info: _batches) {
    info.second->running = 0;
    info.second->in_queue = 0;
    info.second->waiting_end.notify_all();
  }

  for (std::thread &t: _threads) {
    // notify all the threads waiting for
    // some new tasks in the queue and
    // wait for their termination
    _waiting_task.notify_all();
    if (t.joinable()) {
      t.join();
    }
  }
}

/**
 * @brief Re-initialize the pool
 *
 * @param[in] num_of_threads is the new number of threads
 */
void ThreadPool::reinit(const unsigned int &num_of_threads)
{
  std::unique_lock<std::mutex> lock(_mutex);

  _terminating = false;
  _batches.clear();
  _batches.emplace((BatchId)0, std::make_shared<BatchInfo>());

  _threads = std::vector<std::thread>();
  for (unsigned int i = 0; i < num_of_threads; ++i) {
    _threads.emplace_back(&ThreadPool::consumer_loop, this, _threads.size());
  }
}

/**
 * @brief Add some threads to the pool
 *
 * @param[in] num_of_new_threads is the number of new threads
 */
void ThreadPool::add_new_threads(const unsigned int num_of_new_threads)
{
  std::unique_lock<std::mutex> lock(_mutex);

  for (unsigned int i = 0; i < num_of_new_threads; ++i) {
    _threads.emplace_back(&ThreadPool::consumer_loop, this, _threads.size());
  }
}

/**
 * @brief Reset the pool
 *
 * @param[in] num_of_threads is the new number of threads
 */
void ThreadPool::reset(const unsigned int num_of_threads)
{
  terminate();
  reinit(num_of_threads);
}

/**
 * @brief Destroy the thread pool
 */
ThreadPool::~ThreadPool()
{
  terminate();
}
