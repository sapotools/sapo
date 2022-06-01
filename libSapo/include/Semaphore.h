/**
 * @file Semaphore.h
 * @author Alberto Casagrande <acasagrande@units.it>
 * @brief `Semaphore` definition
 * @version 0.1
 * @date 2021-12-10
 *
 * @copyright Copyright (c) 2021-2022
 */

#ifndef SEMAPHORE_H_
#define SEMAPHORE_H_

#include <mutex>
#include <condition_variable>
#include <chrono>

/**
 * @brief A C++11 semaphore implementation.
 */
class Semaphore
{
private:
  unsigned int _tickets; //!< The semaphore counter
  std::mutex _mutex;     //!< A mutex for mutual exclusive operations
  std::condition_variable _condition; //!< The condition variable

public:
  /**
   * @brief Build a new Semaphore object
   */
  Semaphore(): _tickets(0) {}

  /**
   * @brief Build a new Semaphore object
   *
   * @param tickets is the number of available semaphore tickers.
   */
  Semaphore(const unsigned int tickets): _tickets(tickets) {}

  /**
   * @brief Reset number of tickets
   *
   * @param tickets the new number of tickets
   */
  void reset(const unsigned int tickets)
  {
    std::unique_lock<std::mutex> lock(_mutex);

    int delta = ((int)tickets) - _tickets;

    _tickets = tickets;

    while (delta > 0) {
      _condition.notify_one();

      delta--;
    }
  }

  /**
   * @brief Get the number of available tickets.
   *
   * @return The number of available tickets, i.e.,
   *         the counter.
   */
  inline unsigned int available_tickets()
  {
    std::unique_lock<std::mutex> lock(_mutex);

    return _tickets;
  }

  /**
   * @brief Wait for a given number of tickets.
   *
   * This method does not decrement the number of available tickets:
   * it just waits until the semaphore owns at least `tickets`
   * tickets.
   *
   * @param tickets is the number of aimed available tickets.
   */
  void wait_for(const unsigned int tickets)
  {
    std::unique_lock<std::mutex> lock(_mutex);
    _condition.wait(lock, [&]() -> bool { return _tickets >= tickets; });
  }

  /**
   * @brief Wait for a ticket and reserve it.
   *
   * This method implements the P/wait function.
   */
  void reserve()
  {
    std::unique_lock<std::mutex> lock(_mutex);
    _condition.wait(lock, [&]() -> bool { return _tickets > 0; });
    --_tickets;
  }

  /**
   * @brief Release a ticket.
   *
   * This method implements the V/signal function.
   */
  void release()
  {
    std::unique_lock<std::mutex> lock(_mutex);
    ++_tickets;
    _condition.notify_one();
  }
};

#endif // SEMAPHORE_H