/**
 * @file Approximation.h
 * @author Alberto Casagrande <acasagrande@units.it>
 * @brief Interval approximations of numeric types
 * @version 0.1
 * @date 2023-02-17
 *
 * @copyright Copyright (c) 2023
 */

#ifndef _APPROXIMATION_H_
#define _APPROXIMATION_H_

#include <limits>
#include <cmath>

#include "FloatingPoints.h"
#include "ErrorHandling.h"

/**
 * @brief A structure to distinguish punctual types from interval types
 *
 * This structure is only meant to distinguish punctual arithmetic types,
 * like `float`, `double`, or `int`, from interval-based types, such as
 * `Approximation`. This is required to disable `operator==` and
 * `operator!=` between interval-based vector because of their
 * puzzling semantics.
 * In order to use the `LinearAlgebra` implementation of `operator==`
 * and `operator!=` to compare vectors having a non-standard scalar
 * type, please, specialize `is_punctual`.
 *
 * @tparam T is the type to be tested
 */
template<typename T>
struct is_punctual : public std::is_arithmetic<T> {
};

/**
 * @brief Interval-based approximation for numeric types
 *
 * @tparam T is a punctual numeric type
 */
template<typename T,
         typename = typename std::enable_if<is_punctual<T>::value>::type>
class Approximation;

/**
 * @brief The over-approximated absolute value of an approximation
 *
 * @tparam T is a punctual numeric type
 * @param a is an approximation
 * @return the over-approximation of the absolute value of `a`
 */
template<typename T>
Approximation<T> abs(Approximation<T> &&a);

/**
 * @brief The over-approximated square root of an approximation
 *
 * @tparam T is a punctual numeric type
 * @param a is an approximation
 * @return the over-approximation of the square root of `a`
 */
template<typename T>
Approximation<T> sqrt(Approximation<T> &&a);

/**
 * @brief The over-approximated quotient of two approximations
 *
 * @tparam T is a punctual numeric type
 * @param a is an approximation
 * @param b is an approximation
 * @return an over-approximation of \f$a/b\f$
 */
template<typename T>
Approximation<T> operator/(const Approximation<T> &a, Approximation<T> &&b);

/**
 * @brief Approximation for floating point numbers
 *
 * @tparam T is a punctual numeric type
 *
 * @todo This template exploits `std::nextafter` to blindly lower
 *       and upper bounds the result of arithmetic operations.
 *       This approach is less accurate than using floating point
 *       approximations in combination with `std::fegetround()`
 *       and `std::fesetround()`, however, it is much faster
 *       (about 10 time faster in a MacBook Pro 2020). This
 *       performance drop is due to the contest switch required
 *       by `std::fesetround()`.
 */
template<typename T, typename>
class Approximation
{
  T _lower_bound; //!< the approximation lower bound
  T _upper_bound; //!< the approximation upper bound

  /**
   * @brief Set the approximation bounds
   *
   * @param lower_bound is the new lower bound
   * @param upper_bound is the new upper bound
   * @return a reference to the updated approximation
   */
  Approximation<T> &set_bounds(const T &lower_bound, const T &upper_bound);

  /**
   * @brief Set the approximation bounds and approximate them
   *
   * @param lower_bound is the new lower bound (to be under-approximated)
   * @param upper_bound is the new upper bound (to be over-approximated)
   * @return a reference to the updated approximation
   */
  Approximation<T> &set_bounds_and_approximate(const T &lower_bound,
                                               const T &upper_bound);

public:
  /**
   * @brief Construct an approximation for zero
   */
  Approximation();

  /**
   * @brief Construct a new approximation
   *
   * @param lower_bound is the approximation lower bound
   * @param upper_bound is the approximation upper bound
   */
  Approximation(const T &lower_bound, const T &upper_bound);

  /**
   * @brief Construct a new approximation
   *
   * @param value is the value of both upper and lower bound
   */
  Approximation(const T &value);

  /**
   * @brief Get the approximation lower bound
   *
   * @return the approximation lower bound
   */
  const T &lower_bound() const;

  /**
   * @brief Get the approximation upper bound
   *
   * @return the approximation upper bound
   */
  const T &upper_bound() const;

  /**
   * @brief Test whether an approximation contains another approximation
   *
   * @param approximation is an approximation
   * @return `true` if and only if this approximation completely contains
   *       `approximation`
   */
  bool contains(const Approximation<T> &approximation) const;

  /**
   * @brief Test whether an approximation strictly contains another
   * approximation
   *
   * @param approximation is an approximation
   * @return `true` if and only if this approximation strictly contains
   *       `approximation`
   */
  bool strictly_contains(const Approximation<T> &approximation) const;

  /**
   * @brief Test whether an approximation contains a value
   *
   * @param value is an approximation
   * @return `true` if and only if this approximation contains `value`
   */
  bool contains(const T &value) const;

  /**
   * @brief Test whether an approximation strictly contains a value
   *
   * @param value is an approximation
   * @return `true` if and only if this approximation strictly contains
   *       `value`
   */
  bool strictly_contains(const T &value) const;

  /**
   * @brief Get the approximation error
   *
   * @return the approximation error
   */
  const T error() const;

  /**
   * @brief Negate the current approximation
   *
   * This method computes the negative approximation of the current
   * object and sets the current object to it.
   *
   * @return a reference to current object whose value is set to `|*this|`
   */
  Approximation<T> &negate();

  /**
   * @brief The over-approximated in-place sum of two approximations
   *
   * This method computes the over-approximated sum of the current object and
   * another approximation and sets the current object to it.
   *
   * @param a is an approximation
   * @return a reference to the updated object whose value is set to `*this+a`
   */
  Approximation<T> &operator+=(const Approximation<T> &a);

  /**
   * @brief The over-approximated in-place difference of two approximations
   *
   * This method computes the over-approximated difference of the current
   * object and another approximation and sets the current object to it.
   *
   * @param a is an approximation
   * @return a reference to the updated object whose value is set to `*this-a`
   */
  Approximation<T> &operator-=(const Approximation<T> &a);

  /**
   * @brief The over-approximated in-place product of two approximations
   *
   * This method computes the over-approximated product of the current object
   * and another approximation and sets the current object to it.
   *
   * @param a is an approximation
   * @return a reference to the updated object whose value is set to
   * `(*this)*a`
   * @todo tight approximation
   */
  Approximation<T> &operator*=(const Approximation<T> &a);

  /**
   * @brief The over-approximated in-place quotient of two approximations
   *
   * This method computes the over-approximated quotient of the current object
   * and another approximation and sets the current object to it.
   *
   * @param a is an approximation
   * @return a reference to the updated object whose value is set to
   * `(*this)/a`
   * @todo tight approximation
   */
  Approximation<T> &operator/=(const Approximation<T> &a);

  friend Approximation<T> abs<T>(Approximation<T> &&a);
  friend Approximation<T> sqrt<T>(Approximation<T> &&a);

  friend Approximation<T> operator/<T>(const Approximation<T> &a,
                                       Approximation<T> &&b);
};

template<typename T, typename B>
inline Approximation<T> &Approximation<T, B>::set_bounds(const T &lower_bound,
                                                         const T &upper_bound)
{
  _lower_bound = lower_bound;
  _upper_bound = upper_bound;

  return *this;
}

template<typename T, typename B>
inline Approximation<T> &
Approximation<T, B>::set_bounds_and_approximate(const T &lower_bound,
                                                const T &upper_bound)
{
  if constexpr (std::is_floating_point_v<T>) {
    this->set_bounds(
        nextafter(lower_bound, -std::numeric_limits<T>::infinity()),
        nextafter(upper_bound, std::numeric_limits<T>::infinity()));
  } else {
    this->set_bounds(lower_bound, upper_bound);
  }

  return *this;
}

template<typename T, typename B>
Approximation<T, B>::Approximation(): Approximation(0, 0)
{
}

template<typename T, typename B>
Approximation<T, B>::Approximation(const T &lower_bound, const T &upper_bound):
    _lower_bound(lower_bound), _upper_bound(upper_bound)
{
  if (_upper_bound < _lower_bound) {
    SAPO_ERROR("Lower bound must be lesser than or "
               "equal to upper bound",
               std::domain_error);
  }
}

template<typename T, typename B>
Approximation<T, B>::Approximation(const T &value): Approximation(value, value)
{
}

template<typename T, typename B>
inline const T &Approximation<T, B>::lower_bound() const
{
  return _lower_bound;
}

template<typename T, typename B>
inline const T &Approximation<T, B>::upper_bound() const
{
  return _upper_bound;
}

template<typename T, typename B>
inline bool
Approximation<T, B>::contains(const Approximation<T> &approximation) const
{
  return _lower_bound <= approximation.lower_bound()
         && _upper_bound >= approximation.upper_bound();
}

template<typename T, typename B>
inline bool Approximation<T, B>::strictly_contains(
    const Approximation<T> &approximation) const
{
  return _lower_bound < approximation.lower_bound()
         && approximation.upper_bound() < _upper_bound;
}

template<typename T, typename B>
inline bool Approximation<T, B>::contains(const T &value) const
{
  return _lower_bound <= value && value <= _upper_bound;
}

template<typename T, typename B>
inline bool Approximation<T, B>::strictly_contains(const T &value) const
{
  return _lower_bound < value && value < _upper_bound;
}

template<typename T, typename B>
inline const T Approximation<T, B>::error() const
{
  return subtract(_upper_bound, _lower_bound, FE_UPWARD);
}

template<typename T, typename B>
inline Approximation<T> &Approximation<T, B>::negate()
{
  this->set_bounds(-_upper_bound, -_lower_bound);

  return *this;
}

template<typename T, typename B>
inline Approximation<T> &
Approximation<T, B>::operator+=(const Approximation<T> &a)
{
  if constexpr (!std::is_floating_point_v<T>) {
    this->_lower_bound += a._lower_bound;
    this->_upper_bound += a._upper_bound;
  } else {
    this->_lower_bound = add(_lower_bound, a._lower_bound, FE_DOWNWARD);
    this->_upper_bound = add(_upper_bound, a._upper_bound, FE_UPWARD);
  }

  return *this;
}

template<typename T, typename B>
inline Approximation<T> &
Approximation<T, B>::operator-=(const Approximation<T> &a)
{
  if (a == 0) {
    return *this;
  }

  if (this == &a) {
    this->set_bounds(0, 0);

    return *this;
  }

  if constexpr (!std::is_floating_point_v<T>) {
    this->_lower_bound -= a._upper_bound;
    this->_upper_bound += a._lower_bound;
  } else {
    this->_lower_bound = subtract(_lower_bound, a._upper_bound, FE_DOWNWARD);
    this->_upper_bound = subtract(_upper_bound, a._lower_bound, FE_UPWARD);
  }

  return *this;
}

template<typename T, typename B>
Approximation<T> &Approximation<T, B>::operator*=(const Approximation<T> &a)
{
  if (a == 0) {
    this->set_bounds(0, 0);
    return *this;
  }
  if (*this == 0) {
    return *this;
  }

  if (*this > 0) {
    if (a > 0) {
      // [*this_l, *this_u] >= 0 && [a_l , a_u] >= 0:
      // min is this_l*a_l and max is this_u*a_u
      if constexpr (!std::is_floating_point_v<T>) {
        return set_bounds(this->lower_bound() * a.lower_bound(),
                          this->upper_bound() * a.upper_bound());
      } else {
        return set_bounds(
            multiply(this->lower_bound(), a.lower_bound(), FE_DOWNWARD),
            multiply(this->upper_bound(), a.upper_bound(), FE_UPWARD));
      }
    }

    if (a.upper_bound() >= 0) {
      // [*this_l, *this_u] >= 0 && a_l < 0 &&  a_u >= 0:
      // min is this_u*a_l and max is this_u*a_u
      if constexpr (!std::is_floating_point_v<T>) {
        return set_bounds(this->upper_bound() * a.lower_bound(),
                          this->upper_bound() * a.upper_bound());
      } else {
        return set_bounds(
            multiply(this->upper_bound(), a.lower_bound(), FE_DOWNWARD),
            multiply(this->upper_bound(), a.upper_bound(), FE_UPWARD));
      }
    }

    // [*this_l, *this_u] >= 0 && [a_l , a_u] < 0:
    // min is this_u*a_l and max is this_l*a_u
    if constexpr (!std::is_floating_point_v<T>) {
      return set_bounds(this->upper_bound() * a.lower_bound(),
                        this->lower_bound() * a.upper_bound());
    } else {
      return set_bounds(
          multiply(this->upper_bound(), a.lower_bound(), FE_DOWNWARD),
          multiply(this->lower_bound(), a.upper_bound(), FE_UPWARD));
    }
  }

  if (*this < 0) {
    // compute -((-*this)*a)

    this->negate();
    *this *= a;
    return this->negate();
  }

  // this->contains(0)
  if (a > 0) {
    // *this_l < 0 && *this_u >= 0 && [a_l, a_u] >= 0:
    // min is *this_l*a_u and max is *this_u*a_u
    if constexpr (!std::is_floating_point_v<T>) {
      return set_bounds(this->lower_bound() * a.upper_bound(),
                        this->upper_bound() * a.upper_bound());
    } else {
      return set_bounds(
          multiply(this->lower_bound(), a.upper_bound(), FE_DOWNWARD),
          multiply(this->upper_bound(), a.upper_bound(), FE_UPWARD));
    }
  }

  if (a.upper_bound() >= 0) {
    // *this_l < 0 && *this_u >= 0 && a_l < 0 && a_u >= 0:
    // min is min(*this_l*a_u, *this_u*a_l) and max is max(*this_u*a_u,
    // *this_l*a_l)
    using namespace std;

    T new_lower, new_upper;
    if constexpr (!std::is_floating_point_v<T>) {
      new_lower = min(lower_bound() * a.upper_bound(),
                      upper_bound() * a.lower_bound());
      new_upper = max(lower_bound() * a.lower_bound(),
                      upper_bound() * a.upper_bound());
    } else {
      new_lower
          = min(multiply(this->lower_bound(), a.upper_bound(), FE_DOWNWARD),
                multiply(this->upper_bound(), a.lower_bound(), FE_DOWNWARD));
      new_upper
          = max(multiply(this->lower_bound(), a.lower_bound(), FE_UPWARD),
                multiply(this->upper_bound(), a.upper_bound(), FE_UPWARD));
    }
    return this->set_bounds(new_lower, new_upper);
  }

  // *this_l < 0 && *this_u >= 0 && [a_l, a_u] < 0:
  // min is *this_u*a_l and max is *this_l*a_l
  if constexpr (!std::is_floating_point_v<T>) {
    return set_bounds(this->upper_bound() * a.lower_bound(),
                      this->lower_bound() * a.lower_bound());
  } else {
    return set_bounds(
        multiply(this->upper_bound(), a.lower_bound(), FE_DOWNWARD),
        multiply(this->lower_bound(), a.lower_bound(), FE_UPWARD));
  }
}

template<typename T, typename B>
Approximation<T> &Approximation<T, B>::operator/=(const Approximation<T> &a)
{
  if (a.contains(0)) {
    // divisor contains 0
    SAPO_ERROR("division by 0", std::runtime_error);
  }

  if (*this == 0) {
    return *this;
  }

  if (this == &a) {
    return this->set_bounds(1, 1);
  }

  if (a > 0) {
    if (*this >= 0) {
      // [a_l, a_u] > 0 && [*this_l, *this_u] >= 0:
      // min is *this_l/a_u and max is *this_u/a_l
      this->set_bounds_and_approximate(lower_bound() / a.upper_bound(),
                                       upper_bound() / a.lower_bound());
    } else if (this->upper_bound() >= 0) {
      // [a_l, a_u] > 0 && *this_l < 0  && *this_u >= 0:
      // min is *this_l/a_l and max is *this_u/a_l
      this->set_bounds_and_approximate(lower_bound() / a.lower_bound(),
                                       upper_bound() / a.lower_bound());
    } else {
      // [a_l, a_u] > 0 && [*this_l, *this_u] < 0:
      // min is *this_l/a_l and max is *this_u/a_u
      this->set_bounds_and_approximate(lower_bound() / a.lower_bound(),
                                       upper_bound() / a.upper_bound());
    }
    return *this;
  }

  // !a.contains(0)&&a>0, hence, a<0

  // compute -(*this/(-a))
  *this /= (-a);

  return this->negate();
}

/**
 * @brief Test whether an approximations equals a value
 *
 * @tparam T is a punctual numeric type
 * @param a is an approximation
 * @param value is a value
 * @return `true` if and only if any possible value of `a` equals `value`
 */
template<typename T>
inline bool operator==(const Approximation<T> &a, const T &value)
{
  return a.lower_bound() == value && a.upper_bound() == value;
}

/**
 * @brief Test whether an approximations equals a value
 *
 * @tparam T is a punctual numeric type
 * @tparam K is a numeric type
 * @param a is an approximation
 * @param value is a value
 * @return `true` if and only if any possible value of `a` equals `value`
 */
template<typename T, typename K,
         typename = typename std::enable_if<
             !(std::is_same<K, Approximation<T>>::value)>::type>
inline bool operator==(const Approximation<T> &a, const K &value)
{
  return a == T(value);
}

/**
 * @brief Test whether a value equals an approximations
 *
 * @tparam T is a punctual numeric type
 * @tparam K is a numeric type
 * @param value is a value
 * @param a is an approximation
 * @return `true` if and only if `value` equals any possible value of `a`
 */
template<typename T, typename K,
         typename = typename std::enable_if<
             !(std::is_same<K, Approximation<T>>::value)>::type>
inline bool operator==(const K &value, const Approximation<T> &a)
{
  return a == T(value);
}

/**
 * @brief Test whether two approximations are the same
 *
 * @tparam T is a punctual numeric type
 * @param a is an approximation
 * @param b is an approximation
 * @return `true` if and only if `a` and `b` are the same
 */
template<typename T>
inline bool operator==(const Approximation<T> &a, const Approximation<T> &b)
{
  return a.lower_bound() == b.lower_bound()
         && a.upper_bound() == b.upper_bound();
}

/**
 * @brief Test whether an approximations differs from a value
 *
 * @tparam T is a punctual numeric type
 * @param a is an approximation
 * @param value is a value
 * @return `true` if and only if some possible values of `a` differ from
 * `value`
 */
template<typename T>
inline bool operator!=(const Approximation<T> &a, const T &value)
{
  return !(a == value);
}

/**
 * @brief Test whether an approximations differs from a value
 *
 * @tparam T is a punctual numeric type
 * @tparam K is a numeric type
 * @param a is an approximation
 * @param value is a value
 * @return `true` if and only if some possible values of `a` differ from
 * `value`
 */
template<typename T, typename K,
         typename = typename std::enable_if<
             !(std::is_same<K, Approximation<T>>::value)>::type>
inline bool operator!=(const Approximation<T> &a, const K &value)
{
  return !(a == T(value));
}

/**
 * @brief Test whether a value differs from an approximations
 *
 * @tparam T is a punctual numeric type
 * @param value is a value
 * @param a is an approximation
 * @return `true` if and only if `value` differs from some possible values of
 * `a`
 */
template<typename T>
inline bool operator!=(const T &value, const Approximation<T> &a)
{
  return !(a == value);
}

/**
 * @brief Test whether a value differs from an approximations
 *
 * @tparam T is a punctual numeric type
 * @tparam K is a numeric type
 * @param value is a value
 * @param a is an approximation
 * @return `true` if and only if `value` differs from some possible values of
 * `a`
 */
template<typename T, typename K,
         typename = typename std::enable_if<
             !(std::is_same<K, Approximation<T>>::value)>::type>
inline bool operator!=(const K &value, const Approximation<T> &a)
{
  return !(a == T(value));
}

/**
 * @brief Test whether two approximations differ
 *
 * @tparam T is a punctual numeric type
 * @param a is an approximation
 * @param b is an approximation
 * @return `true` if and only if `a` and `b` differ
 */
template<typename T>
inline bool operator!=(const Approximation<T> &a, const Approximation<T> &b)
{
  return !(a == b);
}

/**
 * @brief The "lesser than" relation between an approximation and a value
 *
 * @tparam T is a punctual numeric type
 * @param a is an approximation
 * @param value is a value
 * @return `true` if and only if any possible value of `a` is lesser than
 * `value`
 */
template<typename T>
inline bool operator<(const Approximation<T> &a, const T &value)
{
  return a.upper_bound() < value;
}

/**
 * @brief The "lesser than" relation between a value and an approximation
 *
 * @tparam T is a punctual numeric type
 * @param a is an approximation
 * @param value is a value
 * @return `true` if and only if `value` is lesser than any possible value of
 * `a`
 */
template<typename T>
inline bool operator<(const T &value, const Approximation<T> &a)
{
  return value < a.lower_bound();
}

/**
 * @brief The "lesser than or equal to" relation between an approximation and a
 * value
 *
 * @tparam T is a punctual numeric type
 * @param a is an approximation
 * @param value is a value
 * @return `true` if and only if any possible value of `a` are lesser than or
 *          equal to `value`
 */
template<typename T>
inline bool operator<=(const Approximation<T> &a, const T &value)
{
  return a.upper_bound() <= value;
}

/**
 * @brief The "lesser than or equal to" relation between a value and an
 * approximation
 *
 * @tparam T is a punctual numeric type
 * @param a is an approximation
 * @param value is a value
 * @return `true` if and only if `value` is lesser than or equal to any
 * possible value of `a`
 */
template<typename T>
inline bool operator<=(const T &value, const Approximation<T> &a)
{
  return value <= a.lower_bound();
}

/**
 * @brief The "greater than" relation between an approximation and a value
 *
 * @tparam T is a punctual numeric type
 * @param value is a value
 * @param a is an approximation
 * @return `true` if and only if and only if any possible value of `a` is
 * greater than `value`
 */
template<typename T>
inline bool operator>(const Approximation<T> &a, const T &value)
{
  return a.lower_bound() > value;
}

/**
 * @brief The "greater than" relation between a value and an approximation
 *
 * @tparam T is a punctual numeric type
 * @param value is a value
 * @param a is an approximation
 * @return `true` if and only if `value` is greater than any possible value of
 * `a`
 */
template<typename T>
inline bool operator>(const T &value, const Approximation<T> &a)
{
  return value > a.upper_bound();
}

/**
 * @brief The "greater than or equal to" relation between an approximation and
 * a value
 *
 * @tparam T is a punctual numeric type
 * @param value is a value
 * @param a is an approximation
 * @return `true` if and only if any possible value of `a` are greater than or
 *          equal to `value`
 */
template<typename T>
inline bool operator>=(const Approximation<T> &a, const T &value)
{
  return a.lower_bound() >= value;
}

/**
 * @brief The "greater than or equal to" relation between a value and an
 * approximation
 *
 * @tparam T is a punctual numeric type
 * @param value is a value
 * @param a is an approximation
 * @return `true` if and only if `value` is greater than or equal to any
 * possible value of `a`
 */
template<typename T>
inline bool operator>=(const T &value, const Approximation<T> &a)
{
  return value >= a.upper_bound();
}

/**
 * @brief The "lesser than" relation between approximations
 *
 * @tparam T is a punctual numeric type
 * @param a is an approximation
 * @param b is an approximation
 * @return `true` if and only if any possible value of `a` is lesser
 *     than any possible value of `b`
 */
template<typename T>
inline bool operator<(const Approximation<T> &a, const Approximation<T> &b)
{
  return a.upper_bound() < b.lower_bound();
}

/**
 * @brief The "lesser than or equal to" relation between approximations
 *
 * @tparam T is a punctual numeric type
 * @param a is an approximation
 * @param b is an approximation
 * @return `true` if and only if any possible value of `a` is lesser than or
 *      equal to any possible value of `b`
 */
template<typename T>
inline bool operator<=(const Approximation<T> &a, const Approximation<T> &b)
{
  return a.upper_bound() <= b.lower_bound();
}

/**
 * @brief The "greater than or equal to" relation between approximations
 *
 * @tparam T is a punctual numeric type
 * @param a is an approximation
 * @param b is an approximation
 * @return `true` if and only if any possible value of `a` is greater than or
 *     equal to any possible value of `b`
 */
template<typename T>
inline bool operator>(const Approximation<T> &a, const Approximation<T> &b)
{
  return a.lower_bound() > b.upper_bound();
}

/**
 * @brief The "greater than or equal to" relation between a value and an
 * approximation
 *
 * @tparam T is a punctual numeric type
 * @tparam K is a numeric type
 * @param value is a value
 * @param a is an approximation
 * @return `true` if and only if `value` is greater than or equal to any
 * possible value of `a`
 */
template<typename T, typename K,
         typename = typename std::enable_if_t<
             !(std::is_same<K, T>::value
               || std::is_same<K, Approximation<T>>::value)>>
inline bool operator>=(const K &value, const Approximation<T> &a)
{
  return T(value) >= a;
}

/**
 * @brief The "greater than or equal to" relation between an
 * approximation and a value
 *
 * @tparam T is a punctual numeric type
 * @tparam K is a numeric type
 * @param a is an approximation
 * @param value is a value
 * @return `true` if and only if any possible value of `a` is
 * greater than or equal to `value`
 */
template<typename T, typename K,
         typename = typename std::enable_if_t<
             !(std::is_same<K, T>::value
               || std::is_same<K, Approximation<T>>::value)>>
inline bool operator>=(const Approximation<T> &a, const K &value)
{
  return a >= T(value);
}

/**
 * @brief The "lesser than" relation between a value and an
 * approximation
 *
 * @tparam T is a punctual numeric type
 * @tparam K is a numeric type
 * @param value is a value
 * @param a is an approximation
 * @return `true` if and only if `value` is lesser than any
 * possible value of `a`
 */
template<typename T, typename K,
         typename = typename std::enable_if_t<
             !(std::is_same<K, T>::value
               || std::is_same<K, Approximation<T>>::value)>>
inline bool operator<(const K &value, const Approximation<T> &a)
{
  return T(value) < a;
}

/**
 * @brief The "lesser than" relation between an approximation
 * and a value
 *
 * @tparam T is a punctual numeric type
 * @tparam K is a numeric type
 * @param a is an approximation
 * @param value is a value
 * @return `true` if and only if any possible value of `a` is lesser than
 * `value
 */
template<typename T, typename K,
         typename = typename std::enable_if_t<
             !(std::is_same<K, T>::value
               || std::is_same<K, Approximation<T>>::value)>>
inline bool operator<(const Approximation<T> &a, const K &value)
{
  return a < T(value);
}

/**
 * @brief The "lesser than or equal to" relation between a value and an
 * approximation
 *
 * @tparam T is a punctual numeric type
 * @tparam K is a numeric type
 * @param value is a value
 * @param a is an approximation
 * @return `true` if and only if `value` is lesser than or
 *      equal to any possible value of `a`
 */
template<typename T, typename K,
         typename = typename std::enable_if_t<
             !(std::is_same<K, T>::value
               || std::is_same<K, Approximation<T>>::value)>>
inline bool operator<=(const K &value, const Approximation<T> &a)
{
  return T(value) <= a;
}

/**
 * @brief The "lesser than or equal to" relation between an
 * approximation and a value
 *
 * @tparam T is a punctual numeric type
 * @tparam K is a numeric type
 * @param a is an approximation
 * @param value is a value
 * @return `true` if and only if any possible value of `a` is lesser than or
 *      equal to `value`
 */
template<typename T, typename K,
         typename = typename std::enable_if_t<
             !(std::is_same<K, T>::value
               || std::is_same<K, Approximation<T>>::value)>>
inline bool operator<=(const Approximation<T> &a, const K &value)
{
  return a <= T(value);
}

/**
 * @brief The "greater than or equal to" relation between a value and an
 * approximation
 *
 * @tparam T is a punctual numeric type
 * @tparam K is a numeric type
 * @param value is a value
 * @param a is an approximation
 * @return `true` if and only if `value` is greater than or
 *      equal to any possible value of `a`
 */
template<typename T, typename K,
         typename = typename std::enable_if_t<
             !(std::is_same<K, T>::value
               || std::is_same<K, Approximation<T>>::value)>>
inline bool operator>(const K &value, const Approximation<T> &a)
{
  return T(value) > a;
}

/**
 * @brief The "greater than or equal to" relation between an
 * approximation and a value
 *
 * @tparam T is a punctual numeric type
 * @tparam K is a numeric type
 * @param a is an approximation
 * @param value is a value
 * @return `true` if and only if any possible value of `a` is greater than or
 *      equal to `value`
 */
template<typename T, typename K,
         typename = typename std::enable_if_t<
             !(std::is_same<K, T>::value
               || std::is_same<K, Approximation<T>>::value)>>
inline bool operator>(const Approximation<T> &a, const K &value)
{
  return a > T(value);
}

/**
 * @brief The "greater than or equal to" relation between approximations
 *
 * @tparam T is a punctual numeric type
 * @param a is an approximation
 * @param b is an approximation
 * @return `true` if and only if any possible value of `a` is greater than or
 *      equal to any possible value of `b`
 */
template<typename T>
inline bool operator>=(const Approximation<T> &a, const Approximation<T> &b)
{
  return a.lower_bound() >= b.upper_bound();
}

/**
 * @brief The over-approximated sum of two approximations
 *
 * @tparam T is a punctual numeric type
 * @param a is an approximation
 * @param b is an approximation
 * @return an over-approximation of \f$a+b\f$
 */
template<typename T>
inline Approximation<T> operator+(Approximation<T> &&a,
                                  const Approximation<T> &b)
{
  a += b;

  return std::move(a);
}

/**
 * @brief The over-approximated sum of two approximations
 *
 * @tparam T is a punctual numeric type
 * @param a is an approximation
 * @param b is an approximation
 * @return an over-approximation of \f$a+b\f$
 */
template<typename T>
inline Approximation<T> operator+(const Approximation<T> &a,
                                  Approximation<T> &&b)
{
  b += a;

  return std::move(b);
}

/**
 * @brief The over-approximated sum
 *
 * @tparam T is a punctual numeric type
 * @tparam K is a numeric type
 * @param a is an approximation
 * @param b is a value
 * @return an over-approximation of \f$a+b\f$
 */
template<typename T, typename K,
         typename = typename std::enable_if<
             !(std::is_same<K, Approximation<T>>::value)>::type>
inline Approximation<T> operator+(const Approximation<T> &a, const K &b)
{
  return a + Approximation<T>(b);
}

/**
 * @brief The over-approximated sum
 *
 * @tparam T is a punctual numeric type
 * @tparam K is a numeric type
 * @param a is a value
 * @param b is an approximation
 * @return an over-approximation of \f$a+b\f$
 */
template<typename T, typename K,
         typename = typename std::enable_if<
             !(std::is_same<K, Approximation<T>>::value)>::type>
inline Approximation<T> operator+(const K &a, const Approximation<T> &b)
{
  return Approximation<T>(a) + b;
}

/**
 * @brief The over-approximated sum of two approximations
 *
 * @tparam T is a punctual numeric type
 * @param a is an approximation
 * @param b is an approximation
 * @return an over-approximation of \f$a+b\f$
 */
template<typename T>
inline Approximation<T> operator+(const Approximation<T> &a,
                                  const Approximation<T> &b)
{
  return Approximation<T>(a) + b;
}

/**
 * @brief Negative approximation
 *
 * @tparam T is a punctual numeric type
 * @param a is an approximation
 * @return the approximation of \f$-a\f$
 */
template<typename T>
inline Approximation<T> operator-(Approximation<T> &&a)
{
  a.negate();

  return std::move(a);
}

/**
 * @brief Negative approximation
 *
 * @tparam T is a punctual numeric type
 * @param a is an approximation
 * @return the approximation of \f$-a\f$
 */
template<typename T>
inline Approximation<T> operator-(const Approximation<T> &a)
{
  return -(Approximation<T>(a));
}

/**
 * @brief The over-approximated difference of two approximations
 *
 * @tparam T is a punctual numeric type
 * @param a is an approximation
 * @param b is an approximation
 * @return an over-approximation of \f$a-b\f$
 */
template<typename T>
inline Approximation<T> operator-(const Approximation<T> &a,
                                  Approximation<T> &&b)
{
  b.negate();

  b += a;

  return std::move(b);
}

/**
 * @brief The over-approximated difference of two approximations
 *
 * @tparam T is a punctual numeric type
 * @param a is an approximation
 * @param b is an approximation
 * @return an over-approximation of \f$a-b\f$
 */
template<typename T>
inline Approximation<T> operator-(Approximation<T> &&a,
                                  const Approximation<T> &b)
{
  a -= b;

  return std::move(a);
}

/**
 * @brief The over-approximated difference
 *
 * @tparam T is a punctual numeric type
 * @tparam K is a numeric type
 * @param a is a value
 * @param b is an approximation
 * @return an over-approximation of \f$a-b\f$
 */
template<typename T, typename K,
         typename = typename std::enable_if<
             !(std::is_same<K, Approximation<T>>::value)>::type>
inline Approximation<T> operator-(const K &a, const Approximation<T> &b)
{
  return Approximation<T>(a) - b;
}

/**
 * @brief The over-approximated difference
 *
 * @tparam T is a punctual numeric type
 * @tparam K is a numeric type
 * @param a is an approximation
 * @param b is a value
 * @return an over-approximation of \f$a-b\f$
 */
template<typename T, typename K,
         typename = typename std::enable_if<
             !(std::is_same<K, Approximation<T>>::value)>::type>
inline Approximation<T> operator-(const Approximation<T> &a, const K &b)
{
  return a - Approximation<T>(b);
}

/**
 * @brief The over-approximated difference of two approximations
 *
 * @tparam T is a punctual numeric type
 * @param a is an approximation
 * @param b is an approximation
 * @return an over-approximation of \f$a-b\f$
 */
template<typename T>
Approximation<T> operator-(const Approximation<T> &a,
                           const Approximation<T> &b)
{
  if (&a == &b) {
    return Approximation<T>(0);
  }

  return (-b) + a;
}

/**
 * @brief The over-approximated product of two approximations
 *
 * @tparam T is a punctual numeric type
 * @param a is an approximation
 * @param b is an approximation
 * @return an over-approximation of \f$a*b\f$
 */
template<typename T>
inline Approximation<T> operator*(Approximation<T> &&a,
                                  const Approximation<T> &b)
{
  a *= b;

  return std::move(a);
}

/**
 * @brief The over-approximated product of two approximations
 *
 * @tparam T is a punctual numeric type
 * @param a is an approximation
 * @param b is an approximation
 * @return an over-approximation of \f$a*b\f$
 */
template<typename T>
inline Approximation<T> operator*(const Approximation<T> &a,
                                  Approximation<T> &&b)
{
  return std::move(b) * a;
}

/**
 * @brief The over-approximated product
 *
 * @tparam T is a punctual numeric type
 * @tparam K is a numeric type
 * @param a is an approximation
 * @param b is a value
 * @return an over-approximation of \f$a*b\f$
 */
template<typename T, typename K,
         typename = typename std::enable_if<
             !(std::is_same<K, Approximation<T>>::value)>::type>
inline Approximation<T> operator*(const Approximation<T> &a, const K &b)
{
  return a * Approximation<T>(b);
}

/**
 * @brief The over-approximated product
 *
 * @tparam T is a punctual numeric type
 * @tparam K is a numeric type
 * @param a is a value
 * @param b is an approximation
 * @return an over-approximation of \f$a*b\f$
 */
template<typename T, typename K,
         typename = typename std::enable_if<
             !(std::is_same<K, Approximation<T>>::value)>::type>
inline Approximation<T> operator*(const K &a, const Approximation<T> &b)
{
  return Approximation<T>(a) * b;
}

/**
 * @brief The over-approximated product of two approximations
 *
 * @tparam T is a punctual numeric type
 * @param a is an approximation
 * @param b is an approximation
 * @return an over-approximation of \f$a*b\f$
 */
template<typename T>
inline Approximation<T> operator*(const Approximation<T> &a,
                                  const Approximation<T> &b)
{
  return Approximation<T>(a) * b;
}

/**
 * @brief The over-approximated quotient of two approximations
 *
 * @tparam T is a punctual numeric type
 * @param a is an approximation
 * @param b is an approximation
 * @return an over-approximation of \f$a/b\f$
 */
template<typename T>
Approximation<T> operator/(Approximation<T> &&a, const Approximation<T> &b)
{
  a /= b;

  return std::move(a);
}

/**
 * @brief The over-approximated quotient of two approximations
 *
 * @tparam T is a punctual numeric type
 * @param a is an approximation
 * @param b is an approximation
 * @return an over-approximation of \f$a/b\f$
 */
template<typename T>
Approximation<T> operator/(const Approximation<T> &a, Approximation<T> &&b)
{
  if (b > 0) {
    if (a >= 0) {
      // [b_l, b_u] > 0 && [a_l, a_u] >= 0: min a_l/b_u and max a_u/b_l
      b.set_bounds_and_approximate(a.lower_bound() / b.upper_bound(),
                                   a.upper_bound() / b.lower_bound());
    } else if (a.upper_bound() >= 0) {
      // [b_l, b_u] > 0 && a_l < 0  && a_u >= 0:  min a_l/b_l and max a_u/b_l
      b.set_bounds_and_approximate(a.lower_bound() / b.lower_bound(),
                                   a.upper_bound() / b.lower_bound());
    } else {
      // [b_l, b_u] > 0 && [a_l, a_u] < 0: min a_l/b_l and max a_u/b_u
      b.set_bounds_and_approximate(a.lower_bound() / b.lower_bound(),
                                   a.upper_bound() / b.upper_bound());
    }
    return std::move(b);
  }

  if (b < 0) {
    return -(a / (-std::move(b)));
  }

  // divisor contains 0
  SAPO_ERROR("division by 0", std::runtime_error);
}

/**
 * @brief The over-approximated quotient
 *
 * @tparam T is a punctual numeric type
 * @tparam K is a numeric type
 * @param a is a value
 * @param b is an approximation
 * @return an over-approximation of \f$a/b\f$
 */
template<typename T, typename K,
         typename = typename std::enable_if<
             !(std::is_same<K, Approximation<T>>::value)>::type>
inline Approximation<T> operator/(const K &a, const Approximation<T> &b)
{
  return Approximation<T>(a) / b;
}

/**
 * @brief The over-approximated quotient
 *
 * @tparam T is a punctual numeric type
 * @tparam K is a numeric type
 * @param a is an approximation
 * @param b is a value
 * @return an over-approximation of \f$a/b\f$
 */
template<typename T, typename K,
         typename = typename std::enable_if<
             !(std::is_same<K, Approximation<T>>::value)>::type>
inline Approximation<T> operator/(const Approximation<T> &a, const K &b)
{
  return a / Approximation<T>(b);
}

/**
 * @brief The over-approximated quotient of two approximations
 *
 * @tparam T is a punctual numeric type
 * @param a is an approximation
 * @param b is an approximation
 * @return an over-approximation of \f$a/b\f$
 */
template<typename T>
inline Approximation<T> operator/(const Approximation<T> &a,
                                  const Approximation<T> &b)
{
  if (b.contains(0)) {
    SAPO_ERROR("division by 0", std::runtime_error);
  }

  if (&a == &b) {
    return Approximation<T>(1);
  }

  return Approximation<T>(a) / b;
}

/**
 * @brief The over-approximated absolute value of an approximation
 *
 * @tparam T is a punctual numeric type
 * @param a is an approximation
 * @return the over-approximation of the absolute value of `a`
 */
template<typename T>
Approximation<T> abs(Approximation<T> &&a)
{
  if (a > 0) {
    return std::move(a);
  }

  if (a < 0) {
    return -a;
  }

  a._upper_bound = max(a.upper_bound(), -a.lower_bound());
  a._lower_bound = 0;

  return std::move(a);
}

/**
 * @brief The over-approximated absolute value of an approximation
 *
 * @tparam T is a punctual numeric type
 * @param a is an approximation
 * @return the over-approximation of the absolute value of `a`
 */
template<typename T>
inline Approximation<T> abs(const Approximation<T> &a)
{
  return abs(Approximation<T>(a));
}

/**
 * @brief The over-approximated square root of an approximation
 *
 * @tparam T is a punctual numeric type
 * @param a is an approximation
 * @return the over-approximation of the square root of `a`
 */
template<typename T>
Approximation<T> sqrt(Approximation<T> &&a)
{
  auto lower_bound = sqrt(a.lower_bound());

  while (lower_bound * lower_bound > a.lower_bound()) {
    lower_bound = nextafter(lower_bound, -std::numeric_limits<T>::infinity());
  }

  auto upper_bound = sqrt(a.upper_bound());
  while (upper_bound * upper_bound < a.lower_bound()) {
    upper_bound = nextafter(upper_bound, std::numeric_limits<T>::infinity());
  }
  a.set_bounds(lower_bound, upper_bound);

  return std::move(a);
}

/**
 * @brief The over-approximated square root of an approximation
 *
 * @tparam T is a punctual numeric type
 * @param a is an approximation
 * @return the over-approximation of the square root of `a`
 */
template<typename T>
Approximation<T> sqrt(const Approximation<T> &a)
{
  return sqrt(Approximation<T>(a));
}

/**
 * @brief Write in a stream an approximation
 *
 * @tparam T is a punctual numeric type
 * @param os is an output stream
 * @param approximation is an approximation
 * @return a reference to the output stream
 */
template<typename T>
std::ostream &operator<<(std::ostream &os,
                         const Approximation<T> &approximation)
{
  os << "[" << approximation.lower_bound() << ","
     << approximation.upper_bound() << "]";

  return os;
}

#endif // _APPROXIMATION_H_