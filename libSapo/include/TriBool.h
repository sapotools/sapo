/**
 * @file TriBool.h
 * @author Alberto Casagrande <acasagrande@units.it>
 * @brief Tri-Boolean type
 * @version 0.1
 * @date 2023-03-04
 *
 * @copyright Copyright (c) 2023
 */

#ifndef _TRIBOOL_H_
#define _TRIBOOL_H_

#include <utility>

#include "ErrorHandling.h"

/**
 * @brief Tri-Boolean type
 *
 * This class represents ternary logic values: true, false,
 * and uncertain.
 */
struct TriBool {
  /**
   * @brief The possible values of a tri-Boolean state
   */
  enum values { FALSE, TRUE, UNCERTAIN };

private:
  values _value; //!< the tri-Boolean state

public:
  /**
   * @brief The empty constructor
   *
   * This constructor creates the tri-Boolean value `false`.
   */
  TriBool();

  /**
   * @brief A constructor
   *
   * This constructor creates a tri-Boolean value consistent
   * with the parameter.
   *
   * @param value is the value of the tri-Boolean value
   */
  TriBool(const TriBool::values value);

  /**
   * @brief A Boolean-to-tri-Boolean constructor
   *
   * This constructor creates a tri-Boolean value consistent
   * with the Boolean parameter.
   *
   * @param value is the Boolean value that the tri-Boolean
   *      object must assume
   */
  TriBool(const bool value);

  /**
   * @brief A copy constructor
   *
   * This constructor creates a tri-Boolean value consistent
   * with the parameter.
   *
   * @param orig is the tri-Boolean value that the
   *      object must assume
   */
  TriBool(const TriBool &orig);

  /**
   * @brief Assignement operator
   *
   * @param orig is the original tri-Boolean object
   * @return a reference to the updated object
   */
  TriBool &operator=(const TriBool &orig);

  /**
   * @brief Test the value of a `TriBool` object
   *
   * This method checks whether the value of `*this`
   * is exactly the parameter value.
   *
   * @param value
   * @return `true` if and only if `this->_value==value`;
   *       `false`, otherwise.
   */
  bool is(const TriBool::values &value) const;

  /**
   * @brief Check whether the value is `TriBool::TRUE`
   *
   * @return `true` if and only if the value of
   *       `this->_value` is `TriBool::TRUE`
   */
  bool is_true() const;

  /**
   * @brief Check whether the value is `TriBool::FALSE`
   *
   * @return `true` if and only if the value of
   *       `this->_value` is `TriBool::FALSE`
   */
  bool is_false() const;

  /**
   * @brief Check whether the value is `TriBool::UNCERTAIN`
   *
   * @return `true` if and only if the value of
   *       `this->_value` is `TriBool::UNCERTAIN`
   */
  bool is_uncertain() const;

  /**
   * @brief Get the value
   *
   * @return return the value of the current object
   */
  const TriBool::values &value() const;

  /**
   * @brief Equality between tri-Boolean values
   *
   * @param a is a tri-Boolean value
   * @param b is a tri-Boolean value
   * @return a `TriBool` object `res` such that
   *     `res.is_uncertain()` when either
   *     `a.is_uncertain()` or `b.is_uncertain()`.
   *     If `a._value==b._value` and `!a.is_uncertain()`,
   *     then the method returns a `TriBool` object `res`
   *     such that `res.is_true()`. When, instead,
   *     `a._value!=b._value` and `!a.is_uncertain()`
   *     and `!b.is_uncertain()`, the method
   *     returns a `TriBool` object `res` such that
   *     `res.is_false()`.
   */
  friend TriBool operator==(TriBool &&a, const TriBool &b);

  /**
   * @brief Inequality between tri-Boolean values
   *
   * @param a is a tri-Boolean value
   * @param b is a tri-Boolean value
   * @return a `TriBool` object `res` such that
   *     `res.is_uncertain()` when either
   *     `a.is_uncertain()` or `b.is_uncertain()`.
   *     If `a._value==b._value` and `!a.is_uncertain()`,
   *     then the method returns a `TriBool` object `res`
   *     such that `res.is_false()`. When, instead,
   *     `a._value!=b._value` and `!a.is_uncertain()`
   *     and `!b.is_uncertain()`, the method
   *     returns a `TriBool` object `res` such that
   *     `res.is_true()`.
   */
  friend TriBool operator!=(TriBool &&a, const TriBool &b);

  /**
   * @brief Negate a tri-Boolean value
   *
   * @param a is the tri-Boolean value to be negated
   * @return a `TriBool` object `res` such that
   *     `res.is_uncertain()` when `a.is_uncertain()`.
   *     If `a.is_true()`, this method returns
   *     a `TriBool` object `res` such that
   *     `res.is_false()`. When, instead,
   *     `a.is_false()`, this method returns a
   *     `TriBool` object `res` such that
   *     `res.is_true()`.
   */
  friend TriBool operator!(TriBool &&a);

  /**
   * @brief Tri-Boolean logic conjunction
   *
   * @param a is a tri-Boolean value
   * @param b is a tri-Boolean value
   * @return a copy of `b` when either `a.is_uncertain()`
   *     and `b.is_false()` or `a.is_true()`. In the
   *     remaining cases, this method returns `a`.
   */
  friend TriBool operator&&(TriBool &&a, const TriBool &b);

  /**
   * @brief Tri-Boolean logic exclusive disjunction
   *
   * @param a is a tri-Boolean value
   * @param b is a tri-Boolean value
   * @return `a` if `b.is_false()` or `a.is_true()`. A
   *     copy of `b` if `a.is_false()` or `b.is_true()`.
   *     In the remaining cases, a `TriBool` object `res`
   *     such that `res.is_uncertain()`.
   */
  friend TriBool operator||(TriBool &&a, const TriBool &b);
};

/**
 * @brief Equality between tri-Boolean values
 *
 * @param a is a tri-Boolean value
 * @param b is a tri-Boolean value
 * @return a `TriBool` object `res` such that
 *     `res.is_uncertain()` when either
 *     `a.is_uncertain()` or `b.is_uncertain()`.
 *     If `a._value==b._value` and `!a.is_uncertain()`,
 *     then the method returns a `TriBool` object `res`
 *     such that `res.is_true()`. When, instead,
 *     `a._value!=b._value` and `!a.is_uncertain()`
 *     and `!b.is_uncertain()`, the method
 *     returns a `TriBool` object `res` such that
 *     `res.is_false()`.
 */
TriBool operator==(TriBool &&a, const TriBool &b);

/**
 * @brief Equality between tri-Boolean values
 *
 * @param a is a tri-Boolean value
 * @param b is a tri-Boolean value
 * @return a `TriBool` object `res` such that
 *     `res.is_uncertain()` when either
 *     `a.is_uncertain()` or `b.is_uncertain()`.
 *     If `a._value==b._value` and `!a.is_uncertain()`,
 *     then the method returns a `TriBool` object `res`
 *     such that `res.is_true()`. When, instead,
 *     `a._value!=b._value` and `!a.is_uncertain()`
 *     and `!b.is_uncertain()`, the method
 *     returns a `TriBool` object `res` such that
 *     `res.is_false()`.
 */
TriBool operator==(const TriBool &a, TriBool &&b);

/**
 * @brief Equality between tri-Boolean values
 *
 * @param a is a tri-Boolean value
 * @param b is a tri-Boolean value
 * @return if `a` and `b` are reference to the same
 *      object, then a `TriBool` object `res` such that
 *     `res.is_true()`. When either `a.is_uncertain()`
 *     or `b.is_uncertain()`, a  `TriBool` object `res`
 *     such that `res.is_uncertain()`.
 *     If `a._value==b._value` and `!a.is_uncertain()`,
 *     then the method returns a `TriBool` object `res`
 *     such that `res.is_true()`. When, instead,
 *     `a._value!=b._value` and `!a.is_uncertain()`
 *     and `!b.is_uncertain()`, the method
 *     returns a `TriBool` object `res` such that
 *     `res.is_false()`.
 */
TriBool operator==(const TriBool &a, const TriBool &b);

/**
 * @brief Inequality between tri-Boolean values
 *
 * @param a is a tri-Boolean value
 * @param b is a tri-Boolean value
 * @return a `TriBool` object `res` such that
 *     `res.is_uncertain()` when either
 *     `a.is_uncertain()` or `b.is_uncertain()`.
 *     If `a._value==b._value` and `!a.is_uncertain()`,
 *     then the method returns a `TriBool` object `res`
 *     such that `res.is_false()`. When, instead,
 *     `a._value!=b._value` and `!a.is_uncertain()`
 *     and `!b.is_uncertain()`, the method
 *     returns a `TriBool` object `res` such that
 *     `res.is_true()`.
 */
TriBool operator!=(TriBool &&a, const TriBool &b);

/**
 * @brief Inequality between tri-Boolean values
 *
 * @param a is a tri-Boolean value
 * @param b is a tri-Boolean value
 * @return a `TriBool` object `res` such that
 *     `res.is_uncertain()` when either
 *     `a.is_uncertain()` or `b.is_uncertain()`.
 *     If `a._value==b._value` and `!a.is_uncertain()`,
 *     then the method returns a `TriBool` object `res`
 *     such that `res.is_false()`. When, instead,
 *     `a._value!=b._value` and `!a.is_uncertain()`
 *     and `!b.is_uncertain()`, the method
 *     returns a `TriBool` object `res` such that
 *     `res.is_true()`.
 */
TriBool operator!=(const TriBool &a, TriBool &&b);

/**
 * @brief Inequality between tri-Boolean values
 *
 * @param a is a tri-Boolean value
 * @param b is a tri-Boolean value
 * @return if `a` and `b` are reference to the same
 *      object, then a `TriBool` object `res` such that
 *     `res.is_true()`. When either `a.is_uncertain()`
 *     or `b.is_uncertain()`, a  `TriBool` object `res`
 *     such that `res.is_uncertain()`.
 *     If `a._value==b._value` and `!a.is_uncertain()`,
 *     then the method returns a `TriBool` object `res`
 *     such that `res.is_false()`. When, instead,
 *     `a._value!=b._value` and `!a.is_uncertain()`
 *     and `!b.is_uncertain()`, the method
 *     returns a `TriBool` object `res` such that
 *     `res.is_true()`.
 */
TriBool operator!=(const TriBool &a, const TriBool &b);

/**
 * @brief Negate a tri-Boolean value
 *
 * @param a is the tri-Boolean value to be negated
 * @return a `TriBool` object `res` such that
 *     `res.is_uncertain()` when `a.is_uncertain()`.
 *     If `a.is_true()`, this method returns
 *     a `TriBool` object `res` such that
 *     `res.is_false()`. When, instead,
 *     `a.is_false()`, this method returns a
 *     `TriBool` object `res` such that
 *     `res.is_true()`.
 */
TriBool operator!(TriBool &&a);

/**
 * @brief Negate a tri-Boolean value
 *
 * @param a is the tri-Boolean value to be negated
 * @return a `TriBool` object `res` such that
 *     `res.is_uncertain()` when `a.is_uncertain()`.
 *     If `a.is_true()`, this method returns
 *     a `TriBool` object `res` such that
 *     `res.is_false()`. When, instead,
 *     `a.is_false()`, this method returns a
 *     `TriBool` object `res` such that
 *     `res.is_true()`.
 */
TriBool operator!(const TriBool &a);

/**
 * @brief Tri-Boolean logic conjunction
 *
 * @param a is a tri-Boolean value
 * @param b is a tri-Boolean value
 * @return a copy of `b` when either `a.is_uncertain()`
 *     and `b.is_false()` or `a.is_true()`. In the
 *     remaining cases, this method returns `a`.
 */
TriBool operator&&(TriBool &&a, const TriBool &b);

/**
 * @brief Tri-Boolean logic conjunction
 *
 * @param a is a tri-Boolean value
 * @param b is a tri-Boolean value
 * @return a copy of `b` when either `a.is_uncertain()`
 *     and `b.is_false()` or `a.is_true()`. In the
 *     remaining cases, this method returns `a`.
 */
TriBool operator&&(const TriBool &a, TriBool &&b);

/**
 * @brief Tri-Boolean logic conjunction
 *
 * @param a is a tri-Boolean value
 * @param b is a tri-Boolean value
 * @return a copy of `b` when either `a.is_uncertain()`
 *     and `b.is_false()` or `a.is_true()`. In the
 *     remaining cases, this method returns `a`.
 */
TriBool operator&&(const TriBool &a, const TriBool &b);

/**
 * @brief Tri-Boolean logic exclusive disjunction
 *
 * @param a is a tri-Boolean value
 * @param b is a tri-Boolean value
 * @return a copy of `a` if `b.is_false()` or `a.is_true()`.
 *     A copy of `b` if `a.is_false()` or `b.is_true()`.
 *     In the remaining cases, a `TriBool` object `res` such
 *     that `res.is_uncertain()`.
 */
TriBool operator||(TriBool &&a, const TriBool &b);

/**
 * @brief Tri-Boolean logic exclusive disjunction
 *
 * @param a is a tri-Boolean value
 * @param b is a tri-Boolean value
 * @return a copy of `a` if `b.is_false()` or `a.is_true()`.
 *     A copy of `b` if `a.is_false()` or `b.is_true()`.
 *     In the remaining cases, a `TriBool` object `res` such
 *     that `res.is_uncertain()`.
 */
TriBool operator||(const TriBool &a, TriBool &&b);

/**
 * @brief Tri-Boolean logic exclusive disjunction
 *
 * @param a is a tri-Boolean value
 * @param b is a tri-Boolean value
 * @return a copy of `a` if `b.is_false()` or `a.is_true()`.
 *     A copy of `b` if `a.is_false()` or `b.is_true()`.
 *     In the remaining cases, a `TriBool` object `res` such
 *     that `res.is_uncertain()`.
 */
TriBool operator||(const TriBool &a, const TriBool &b);

/**
 * @brief A characterization for Boolean types
 *
 * @tparam T is the type to be tested
 */
template<typename T>
struct is_logic
    : public std::integral_constant<bool,
                                    std::is_same<bool, T>::value
                                        || std::is_same<TriBool, T>::value> {
};

template<class T>
inline constexpr bool is_logic_v = is_logic<T>::value;

/**
 * @brief Test whether a value is true
 *
 * @param bool_value is a Boolean value
 * @return `true` if and only if `bool_value==true`
 */
bool is_true(const bool &bool_value);

/**
 * @brief Test whether a value is false
 *
 * @param bool_value is a Boolean value
 * @return `true` if and only if `bool_value==false`
 */
bool is_false(const bool &bool_value);

/**
 * @brief Test whether a value is uncertain
 *
 * @param bool_value is a Boolean value
 * @return `false`
 */
bool is_uncertain(const bool &bool_value);

/**
 * @brief Test whether a value is true
 *
 * @param bool_value is a tri-Boolean value
 * @return `bool_value.is_true()`
 */
bool is_true(const TriBool &bool_value);

/**
 * @brief Test whether a value is false
 *
 * @param bool_value is a tri-Boolean value
 * @return `bool_value.is_false()`
 */
bool is_false(const TriBool &bool_value);

/**
 * @brief Test whether a value is uncertain
 *
 * @param bool_value is a tri-Boolean value
 * @return `bool_value.is_uncertain()`
 */
bool is_uncertain(const TriBool &bool_value);

/**
 * @brief Test the value of a `TriBool` object
 *
 * This method checks whether the value of `*this`
 * is exactly the parameter value.
 *
 * @param value
 * @return `true` if and only if `this->_value==value`;
 *       `false`, otherwise.
 */
inline bool TriBool::is(const TriBool::values &value) const
{
  return this->_value == value;
}

/**
 * @brief Check whether the value is `TriBool::TRUE`
 *
 * @return `true` if and only if the value of
 *       `this->_value` is `TriBool::TRUE`
 */
inline bool TriBool::is_true() const
{
  return is(TriBool::TRUE);
}

/**
 * @brief Check whether the value is `TriBool::FALSE`
 *
 * @return `true` if and only if the value of
 *       `this->_value` is `TriBool::FALSE`
 */
inline bool TriBool::is_false() const
{
  return is(TriBool::FALSE);
}

/**
 * @brief Check whether the value is `TriBool::UNCERTAIN`
 *
 * @return `true` if and only if the value of
 *       `this->_value` is `TriBool::UNCERTAIN`
 */
inline bool TriBool::is_uncertain() const
{
  return is(TriBool::UNCERTAIN);
}

/**
 * @brief Get the value
 *
 * @return return the value of the current object
 */
inline const TriBool::values &TriBool::value() const
{
  return this->_value;
}

/**
 * @brief Equality between tri-Boolean values
 *
 * @param a is a tri-Boolean value
 * @param b is a tri-Boolean value
 * @return a `TriBool` object `res` such that
 *     `res.is_uncertain()` when either
 *     `a.is_uncertain()` or `b.is_uncertain()`.
 *     If `a._value==b._value` and `!a.is_uncertain()`,
 *     then the method returns a `TriBool` object `res`
 *     such that `res.is_true()`. When, instead,
 *     `a._value!=b._value` and `!a.is_uncertain()`
 *     and `!b.is_uncertain()`, the method
 *     returns a `TriBool` object `res` such that
 *     `res.is_false()`.
 */
inline TriBool operator==(const TriBool &a, TriBool &&b)
{
  return std::move(b) == a;
}

/**
 * @brief Equality between tri-Boolean values
 *
 * @param a is a tri-Boolean value
 * @param b is a tri-Boolean value
 * @return a `TriBool` object `res` such that
 *     `res.is_uncertain()` when either
 *     `a.is_uncertain()` or `b.is_uncertain()`.
 *     If `a._value==b._value` and `!a.is_uncertain()`,
 *     then the method returns a `TriBool` object `res`
 *     such that `res.is_true()`. When, instead,
 *     `a._value!=b._value` and `!a.is_uncertain()`
 *     and `!b.is_uncertain()`, the method
 *     returns a `TriBool` object `res` such that
 *     `res.is_false()`.
 */
inline TriBool operator==(TriBool &&a, TriBool &&b)
{
  return std::move(a) == b;
}

/**
 * @brief Equality between tri-Boolean values
 *
 * @param a is a tri-Boolean value
 * @param b is a tri-Boolean value
 * @return if `a` and `b` are reference to the same
 *      object, then a `TriBool` object `res` such that
 *     `res.is_true()`. When either `a.is_uncertain()`
 *     or `b.is_uncertain()`, a  `TriBool` object `res`
 *     such that `res.is_uncertain()`.
 *     If `a._value==b._value` and `!a.is_uncertain()`,
 *     then the method returns a `TriBool` object `res`
 *     such that `res.is_true()`. When, instead,
 *     `a._value!=b._value` and `!a.is_uncertain()`
 *     and `!b.is_uncertain()`, the method
 *     returns a `TriBool` object `res` such that
 *     `res.is_false()`.
 */
TriBool operator==(const TriBool &a, const TriBool &b);

/**
 * @brief Inequality between tri-Boolean values
 *
 * @param a is a tri-Boolean value
 * @param b is a tri-Boolean value
 * @return a `TriBool` object `res` such that
 *     `res.is_uncertain()` when either
 *     `a.is_uncertain()` or `b.is_uncertain()`.
 *     If `a._value==b._value` and `!a.is_uncertain()`,
 *     then the method returns a `TriBool` object `res`
 *     such that `res.is_false()`. When, instead,
 *     `a._value!=b._value` and `!a.is_uncertain()`
 *     and `!b.is_uncertain()`, the method
 *     returns a `TriBool` object `res` such that
 *     `res.is_true()`.
 */
inline TriBool operator!=(const TriBool &a, TriBool &&b)
{
  return std::move(b) != a;
}

/**
 * @brief Inequality between tri-Boolean values
 *
 * @param a is a tri-Boolean value
 * @param b is a tri-Boolean value
 * @return if `a` and `b` are reference to the same
 *      object, then a `TriBool` object `res` such that
 *     `res.is_true()`. When either `a.is_uncertain()`
 *     or `b.is_uncertain()`, a  `TriBool` object `res`
 *     such that `res.is_uncertain()`.
 *     If `a._value==b._value` and `!a.is_uncertain()`,
 *     then the method returns a `TriBool` object `res`
 *     such that `res.is_false()`. When, instead,
 *     `a._value!=b._value` and `!a.is_uncertain()`
 *     and `!b.is_uncertain()`, the method
 *     returns a `TriBool` object `res` such that
 *     `res.is_true()`.
 */
inline TriBool operator!=(const TriBool &a, const TriBool &b)
{
  if (&a == &b) {
    return false;
  }

  return TriBool(a) != b;
}

/**
 * @brief Negate a tri-Boolean value
 *
 * @param a is the tri-Boolean value to be negated
 * @return a `TriBool` object `res` such that
 *     `res.is_uncertain()` when `a.is_uncertain()`.
 *     If `a.is_true()`, this method returns
 *     a `TriBool` object `res` such that
 *     `res.is_false()`. When, instead,
 *     `a.is_false()`, this method returns a
 *     `TriBool` object `res` such that
 *     `res.is_true()`.
 */
inline TriBool operator!(const TriBool &a)
{
  return !TriBool(a);
}

/**
 * @brief Tri-Boolean logic conjunction
 *
 * @param a is a tri-Boolean value
 * @param b is a tri-Boolean value
 * @return a copy of `b` when either `a.is_uncertain()`
 *     and `b.is_false()` or `a.is_true()`. In the
 *     remaining cases, this method returns `a`.
 */
inline TriBool operator&&(TriBool &&a, TriBool &&b)
{
  return std::move(a) && b;
}

/**
 * @brief Tri-Boolean logic conjunction
 *
 * @param a is a tri-Boolean value
 * @param b is a tri-Boolean value
 * @return a copy of `b` when either `a.is_uncertain()`
 *     and `b.is_false()` or `a.is_true()`. In the
 *     remaining cases, this method returns `a`.
 */
inline TriBool operator&&(const TriBool &a, TriBool &&b)
{
  return std::move(b) && a;
}

/**
 * @brief Tri-Boolean logic conjunction
 *
 * @param a is a tri-Boolean value
 * @param b is a tri-Boolean value
 * @return a copy of `b` when either `a.is_uncertain()`
 *     and `b.is_false()` or `a.is_true()`. In the
 *     remaining cases, this method returns `a`.
 */
inline TriBool operator&&(const TriBool &a, const TriBool &b)
{
  return TriBool(a) && b;
}

/**
 * @brief Tri-Boolean logic exclusive disjunction
 *
 * @param a is a tri-Boolean value
 * @param b is a tri-Boolean value
 * @return a copy of `a` if `b.is_false()` or `a.is_true()`.
 *     A copy of `b` if `a.is_false()` or `b.is_true()`.
 *     In the remaining cases, a `TriBool` object `res` such
 *     that `res.is_uncertain()`.
 */
inline TriBool operator||(const TriBool &a, TriBool &&b)
{
  return std::move(b) || a;
}

/**
 * @brief Tri-Boolean logic exclusive disjunction
 *
 * @param a is a tri-Boolean value
 * @param b is a tri-Boolean value
 * @return a copy of `a` if `b.is_false()` or `a.is_true()`.
 *     A copy of `b` if `a.is_false()` or `b.is_true()`.
 *     In the remaining cases, a `TriBool` object `res` such
 *     that `res.is_uncertain()`.
 */
inline TriBool operator||(const TriBool &a, const TriBool &b)
{
  return TriBool(a) || b;
}

/**
 * @brief Test whether a value is true
 *
 * @param bool_value is a Boolean value
 * @return `true` if and only if `bool_value==true`
 */
inline bool is_true(const bool &bool_value)
{
  return bool_value;
}

/**
 * @brief Test whether a value is false
 *
 * @param bool_value is a Boolean value
 * @return `true` if and only if `bool_value==false`
 */
inline bool is_false(const bool &bool_value)
{
  return !bool_value;
}

/**
 * @brief Test whether a value is uncertain
 *
 * @param bool_value is a Boolean value
 * @return `false`
 */
inline bool is_uncertain(const bool &bool_value)
{
  (void)bool_value; // to avoid `-Wunused-parameter`

  return false;
}

/**
 * @brief Test whether a value is true
 *
 * @param bool_value is a tri-Boolean value
 * @return `bool_value.is_true()`
 */
inline bool is_true(const TriBool &bool_value)
{
  return bool_value.is_true();
}

/**
 * @brief Test whether a value is false
 *
 * @param bool_value is a tri-Boolean value
 * @return `bool_value.is_false()`
 */
inline bool is_false(const TriBool &bool_value)
{
  return bool_value.is_false();
}

/**
 * @brief Test whether a value is uncertain
 *
 * @param bool_value is a tri-Boolean value
 * @return `bool_value.is_uncertain()`
 */
inline bool is_uncertain(const TriBool &bool_value)
{
  return bool_value.is_uncertain();
}

/**
 * @brief Print a tri-Boolean in an output stream
 *
 * @param os is the output stream
 * @param tribool is the tri-Boolean value to stream
 * @return a reference to the output stream
 */
std::ostream &operator<<(std::ostream &os, const TriBool &tribool);

#endif // _TRIBOOL_H_
