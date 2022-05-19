/**
 * @file TimeInterval.h
 * @author Alberto Casagrande <acasagrande@units.it>
 * @brief Temporal interval definition
 * @version 0.1
 * @date 2022-05-04
 * 
 * @copyright Copyright (c) 2021-2022
 * 
 */

#ifndef TIMEINTERVAL_H_
#define TIMEINTERVAL_H_

#include <ostream>

/**
 * @brief Time intervals
 * 
 * This class represents time intervals
 */
class TimeInterval
{
  /**
   * @brief Print a time interval in a stream
   * 
   * @param[in] os is the output stream
   * @param[in] t_itvl is the time interval to be printed
   * @return a reference to the output stream
   */
  friend inline std::ostream &operator<<(std::ostream &os, 
                                         const TimeInterval &t_itvl)
  {
    os << "[" << t_itvl._begin << ", " << t_itvl._end << "]";

    return os;
  }

  int _begin;
  int _end;

public:
  /**
   * Constructor for the time inteval [0,0]
   */
  TimeInterval();
  /**
   * Copy constructor for time intervals
   *
   * @param[in] orig is the time interval model
   */
  TimeInterval(const TimeInterval &orig);

  /**
   * Constructor for singleton time intevals.
   *
   * This constructor builds an object that represents the time interval [time,
   * time], where time is the method parameter.
   *
   * @param[in] time is the time included in the new time interval
   */
  TimeInterval(const int time);

  /**
   * Constructor for time intevals.
   *
   * This constructor builds an object that represents the time interval
   * [begin, end], where begin and end are the method parameters.
   *
   * @param[in] begin is the begin of the new time interval
   * @param[in] end is the end of the new time interval
   */
  TimeInterval(const int begin, const int end);

  /**
   * Sets a new begin for the time interval.
   *
   * @param[in] begin is the new begin of this time interval
   * @return the time interval itself
   */
  const TimeInterval &set_begin(const int begin);

  /**
   * Sets a new end for the time interval.
   *
   * @param[in] end is the new end of this time interval
   * @return the time interval itself
   */
  const TimeInterval &set_end(const int end);

  /**
   * Returns the time interval begin
   *
   * @return the time interval begin
   */
  const int &begin() const
  {
    return this->_begin;
  }

  /**
   * Returns the time interval end
   *
   * @return the time interval end
   */
  const int &end() const
  {
    return this->_end;
  }

  /**
   * Checks whether a time is strictly included in the time interval
   *
   * @param[in] time is the time point to be compared with this interval
   * @return `true` if and only if the provided parameter is strictly
   *         included in the time interval.
   */
  bool strictly_contains(const int time) const
  {
    return (this->begin() < time) && (time < this->end());
  }

  /**
   * Checks whether a time is included in the time interval
   *
   * @param[in] time is the time point to be compared with this interval
   * @return `true` if and only if the provided parameter is included
   *         in the time interval.
   */
  bool contains(const int time) const
  {
    return (this->begin() <= time) && (time <= this->end());
  }

  /**
   * Checks whether the time interval is empty
   *
   * @return `true` if and only if the time interval is empty, i.e.,
   *         the begin of this time interval comes before the end of it
   */
  bool is_empty() const
  {
    return this->begin() > this->end();
  }
};

/**
 * Checks whether a time interval comes before a time point.
 *
 * This method verifies whether all the points in a time interval
 * come before a specified time.
 *
 * @param[in] itvl is the interval to be compared
 * @param[in] time is the time point to be compared
 * @return `true` if and only if the time point comes after
 *         all the times in the interval.
 */
inline bool operator<(const TimeInterval &itvl, const int time)
{
  return itvl.end() < time;
}

/**
 * Checks whether a time interval comes after a time point.
 *
 * This method verifies whether all the points in a time interval
 * come after a specified time.
 *
 * @param[in] itvl is the interval to be compared
 * @param[in] time is the time point to be compared
 * @return `true` if and only if the time point comes before
 *          all the times in the interval.
 */
inline bool operator>(const TimeInterval &itvl, const int time)
{
  return itvl.begin() > time;
}

/**
 * Checks whether a time point comes begin a time interval.
 *
 * This method verifies whether the specified time comes before each of
 * the points in a time interval.
 *
 * @param[in] time is the time point to be compared
 * @param[in] itvl is the interval to be compared
 * @return `true` if and only if the time point comes before
 *          all the times in the interval.
 */
inline bool operator<(const int time, const TimeInterval &itvl)
{
  return itvl > time;
}

/**
 * Checks whether a time point comes after a time interval.
 *
 * This method verifies whether the specified time comes after each of
 * the points in a time interval.
 *
 * @param[in] time is the time point to be compared
 * @param[in] itvl is the interval to be compared
 * @return `true` if and only if the time point comes after
 *          all the times in the interval.
 */
inline bool operator>(const int time, const TimeInterval &itvl)
{
  return itvl < time;
}

#endif // TIMEINTERVAL_H_
