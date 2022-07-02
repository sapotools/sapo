/**
 * @file TimeInterval.cpp
 * @author Alberto Casagrande <acasagrande@units.it>
 * @brief Temporal interval definition
 * @version 0.1
 * @date 2022-05-04
 *
 * @copyright Copyright (c) 2021-2022
 *
 */

#include "STL/TimeInterval.h"

/**
 * Constructor for the time interval [0,0]
 */
TimeInterval::TimeInterval(): _begin(0), _end(0) {}

/**
 * Copy constructor for time intervals
 *
 * @param[in] orig is the time interval model
 */
TimeInterval::TimeInterval(const TimeInterval &orig):
    _begin(orig.begin()), _end(orig.end())
{
}

/**
 * Constructor for singleton time intervals.
 *
 * This constructor builds an object that represents the time interval [time,
 * time], where time is the method parameter.
 *
 * @param[in] time is the time included in the new time interval
 */
TimeInterval::TimeInterval(const int time): _begin(time), _end(time) {}

/**
 * Constructor for time intervals.
 *
 * This constructor builds an object that represents the time interval [begin,
 * end], where begin and end are the method parameters.
 *
 * @param[in] begin is the begin of the new time interval
 * @param[in] end is the end of the new time interval
 */
TimeInterval::TimeInterval(const int begin, const int end):
    _begin(begin), _end(end)
{
}

/**
 * Sets a new begin for the time interval.
 *
 * @param[in] begin is the new begin of this time interval
 * @return the time interval itself
 */
const TimeInterval &TimeInterval::set_begin(const int begin)
{
  this->_begin = begin;

  return *this;
}

/**
 * Sets a new end for the time interval.
 *
 * @param[in] end is the new end of this time interval
 * @return the time interval itself
 */
const TimeInterval &TimeInterval::set_end(const int end)
{
  this->_end = end;

  return *this;
}
