/**
 * @file ProgressAccounter.h
 * @author Alberto Casagrande <acasagrande@units.it>
 * @brief Computation progresses accounting
 * @version 0.1
 * @date 2022-02-17
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef PROGRESSACCOUNTER_H_
#define PROGRESSACCOUNTER_H_

#include <iostream>
#include <functional>

#include <mutex>
#include <shared_mutex>

/**
 * @brief A class to account computation progress
 *
 * The objects of this class take care of accounting
 * the progress of a Sapo computation.
 */
class ProgressAccounter
{
protected:
  unsigned int _performed; //!< number of already performed steps
  unsigned int _expected;  //!< total number of expected steps
  std::mutex _mutex;       //!< A mutex for mutual exclusive operations

public:
  /**
   * @brief A base constructor
   */
  ProgressAccounter(): _performed(0), _expected(0) {}

  /**
   * @brief A ProgressAccounter constructor
   *
   * @param expected_steps is the number of expected steps
   */
  ProgressAccounter(const unsigned int expected_steps):
      _performed(0), _expected(expected_steps)
  {
  }

  /**
   * @brief A ProgressAccounter constructor
   *
   * @param expected_steps is the number of expected steps
   * @param performed_steps is the number of already performed steps
   */
  ProgressAccounter(const unsigned int expected_steps,
                    const unsigned int performed_steps):
      _performed(performed_steps),
      _expected(expected_steps)
  {
  }

  /**
   * @brief Reset the number expected steps
   *
   * @param expected_steps is the new number of expected steps
   * @return the number of expected steps
   */
  const unsigned int &set_expected(const unsigned int expected_steps)
  {
    std::unique_lock<std::mutex> lock(_mutex);

    _expected = expected_steps;

    return _expected;
  }

  /**
   * @brief Set the number of already performed steps
   *
   * @param performed_steps is the new number of already performed steps
   * @return the number of already performed steps
   */
  const unsigned int &set_performed(const unsigned int performed_steps)
  {
    std::unique_lock<std::mutex> lock(_mutex);

    _performed = performed_steps;

    return _performed;
  }

  /**
   * @brief Get the number expected steps
   *
   * @return the number expected steps
   */
  const unsigned int &get_expected()
  {
    std::unique_lock<std::mutex> lock(_mutex);

    return _expected;
  }

  /**
   * @brief Get the number of already performed steps
   *
   * @return the number of already performed steps
   */
  const unsigned int &get_performed()
  {
    std::unique_lock<std::mutex> lock(_mutex);

    return _performed;
  }

  /**
   * @brief Increase the number of performed steps
   *
   * @param delta_steps is the number of newly performed steps
   * @return the number of already performed steps
   */
  virtual const unsigned int &increase_performed(const unsigned int delta_steps
                                                 = 1)
  {
    std::unique_lock<std::mutex> lock(_mutex);

    _performed += delta_steps;

    return _performed;
  }

  /**
   * @brief Increase the number of performed steps
   *
   * @param performed is the number of performed steps
   * @return the number of already performed steps
   */
  virtual const unsigned int &
  increase_performed_to(const unsigned int performed)
  {
    return set_performed(performed);
  }

  virtual ~ProgressAccounter() {}
};

/**
 * @brief A progress bar accounter
 */
class ProgressBar : private ProgressAccounter
{
  const unsigned int _num_of_bar_dots; //!< number of dots in the progress bar
  unsigned int _represented_steps;     //!< number of already represented_steps
  std::ostream &_bstream;              //!< output stream for the progress bar
  std::string _preamble;               //!< bar output preamble

  /**
   * @brief Increase the number of performed steps
   *
   * This method is not thread-safe
   *
   * @param performed is the number of performed steps
   * @return the number of already performed steps
   */
  const unsigned int &
  unsafe_increase_performed_to(const unsigned int performed);

public:
  /**
   * @brief A base constructor
   *
   * @param num_of_bar_dots is the total number of dots in the progress bar
   * @param bar_stream is a reference to the output stream for the progress bar
   */
  ProgressBar(const unsigned int num_of_bar_dots,
              std::ostream &bar_stream = std::ref(std::cerr)):
      ProgressAccounter(),
      _num_of_bar_dots(num_of_bar_dots), _represented_steps(0),
      _bstream(bar_stream), _preamble("")
  {
  }

  /**
   * @brief A ProgressBar constructor
   *
   * @param expected_steps is the number of expected steps
   * @param num_of_bar_dots is the total number of dots in the progress bar
   * @param bar_stream is a reference to the output stream for the progress bar
   */
  ProgressBar(const unsigned int expected_steps,
              const unsigned int num_of_bar_dots,
              std::ostream &bar_stream = std::ref(std::cerr)):
      ProgressAccounter(expected_steps),
      _num_of_bar_dots(num_of_bar_dots), _represented_steps(0),
      _bstream(bar_stream), _preamble("")
  {
  }

  /**
   * @brief Increase the number of performed steps
   *
   * @param delta_steps is the number of newly performed steps
   * @return the number of already performed steps
   */
  const unsigned int &increase_performed(const unsigned int delta_steps = 1);

  /**
   * @brief Increase the number of performed steps
   *
   * @param performed is the number of performed steps
   * @return the number of already performed steps
   */
  const unsigned int &increase_performed_to(const unsigned int performed);
};

#endif // PROGRESSACCOUNTER_H_
