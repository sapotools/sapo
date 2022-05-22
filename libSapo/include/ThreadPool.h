#ifndef THREADPOOL_H
#define THREADPOOL_H

#include <functional>
#include <map>
#include <mutex>
#include <queue>
#include <set>
#include <thread>
#include <vector>
#include <condition_variable>

/**
 * @brief A thread pool to which task can be submitted
 */
class ThreadPool
{
public:
  /**
   * @brief The type of batch identifiers
   *
   */
  typedef unsigned int BatchId;

protected:
  /**
   * @brief Information about batches
   *
   * This class stores the number of submitted and tasks of the batch.
   */
  class BatchInfo
  {
  public:
    unsigned int running;  //!< The number of running tasks in the batch
    unsigned int in_queue; //!< The number of batch tasks in the queue
    std::condition_variable waiting_end; //!< Testifying the batch conclusion

    /**
     * @brief Constructor
     */
    BatchInfo(): running(0), in_queue(0) {}

    /**
     * @brief Copy constructor
     *
     * @param[in] orig is the model for the new object
     */
    BatchInfo(const BatchInfo &orig):
        running(orig.running), in_queue(orig.in_queue)
    {
    }

    /**
     * @brief Return the number of unfinished tasks
     *
     * This method returns the number of unfinished tasks, i.e., the enqueue
     * tasks that have not been finished yet.
     *
     * @return the number of unfinished tasks
     */
    inline unsigned unfinished_tasks() const
    {
      return running + in_queue;
    }
  };

  /**
   * @brief The task type
   *
   * Every task is a function associated to a batch
   * identifier.
   */
  typedef std::pair<std::function<void()>, BatchId> Task;

  std::vector<std::thread> _threads;     //!< The thread list
  std::queue<Task> _queue;               //!< The task queue
  std::map<BatchId, BatchInfo> _batches; //!< Batch info
  std::set<BatchId> _unused_batch_ids;   //!< Not used batch ids
  std::mutex _mutex; //!< A mutex for mutual exclusive operations
  std::condition_variable _waiting_task; //!< Testify that a task is waiting
  bool _terminating;                     //!< The pool is to be destroyed

  /**
   * @brief Get information about a batch
   *
   * @param[in] batch_id is the searched batch id
   * @return a reference to the batch information
   */
  BatchInfo &get_batch_info(const BatchId batch_id);

  /**
   * @brief Extract next task from the queue
   *
   * @return The next task to be executed
   */

  /**
   * @brief Extract next task from the queue
   *
   * @param[out] next is the parameter to be assigned with the next task
   * @param[in] lock is the the lock used to extract the next task
   * @return `true` if and only if `next` has been reassigned
   */
  bool extract_next_task(ThreadPool::Task &next,
                         std::unique_lock<std::mutex> &lock);

  /**
   * @brief The main thread loop
   *
   * @param[in] thread_id is the thread id in the pool
   */
  void consumer_loop(unsigned int thread_id);

  /**
   * @brief Re-initialize the pool
   *
   * @param[in] num_of_threads is the new number of threads
   */
  void reinit(const unsigned int &num_of_thread);

public:
  /**
   * @brief Create a new Thread Pool object
   *
   * @param[in] num_of_threads is the number of thread in the pool
   */
  ThreadPool(const unsigned num_of_threads);

  /**
   * @brief Create a pull of threads according to the hardware
   */
  ThreadPool();

  /**
   * @brief Create a new task batch
   *
   * @return The batch id of the new batch
   */
  BatchId create_batch();

  /**
   * @brief Close the batch and remove it from the pool
   *
   * @param[in] batch_id is the id of the batch to be closed
   */
  void close_batch(const ThreadPool::BatchId batch_id);

  /**
   * @brief Submit a task to the pool
   *
   * @tparam T is the type of the task routine
   * @tparam Ts are the type of the task parameters
   * @param[in] routine is the task routine
   * @param[in] params is the task parameters
   */
  template<typename T, typename... Ts>
  void submit(T &&routine, Ts &&...params)
  {
    submit_to_batch(0, routine, params...);
  }

  /**
   * @brief Submit a task to a batch in the pool
   *
   * @tparam T is the type of the task routine
   * @tparam Ts are the type of the task parameters
   * @param[in] batch_id is the id of the considered batch
   * @param[in] routine is the task routine
   * @param[in] params is the task parameters
   */
  template<typename T, typename... Ts>
  void submit_to_batch(const ThreadPool::BatchId batch_id, T &&routine,
                       Ts &&...params)
  {
    {
      std::unique_lock<std::mutex> lock(_mutex);

      // increase the number of tasks in the batch
      get_batch_info(batch_id).in_queue++;

      // enqueue the task together with its batch id
      _queue.emplace(
          std::bind(std::forward<T>(routine), std::forward<Ts>(params)...),
          batch_id);
    }

    // notify that some task are present in the queue
    _waiting_task.notify_one();
  }

  /**
   * @brief Join the thread pool and complete the batch
   *
   * @param[in] batch_id is the id of the batch that must be completed
   */
  void join_threads(const ThreadPool::BatchId batch_id = 0);

  /**
   * @brief Terminate the pool
   */
  void terminate();

  /**
   * @brief Reset the pool
   *
   * @param[in] num_of_threads is the new number of threads
   */
  void reset(const unsigned int num_of_threads);

  /**
   * @brief Add some threads to the pool
   *
   * @param[in] num_of_new_threads is the number of new threads
   */
  void add_new_threads(const unsigned int num_of_new_threads);

  /**
   * @brief Destroy the thread pool
   */
  ~ThreadPool();
};

#endif // THREADPOOL_H
