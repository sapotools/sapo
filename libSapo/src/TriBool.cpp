/**
 * @file TriBool.cpp
 * @author Alberto Casagrande <acasagrande@units.it>
 * @brief Tri-Boolean type
 * @version 0.1
 * @date 2023-03-05
 *
 * @copyright Copyright (c) 2023
 */

#include "TriBool.h"

/**
 * @brief The empty constructor
 *
 * This constructor creates the tri-Boolean value `false`.
 */
TriBool::TriBool(): _value(values::FALSE) {}

/**
 * @brief A constructor
 *
 * This constructor creates a tri-Boolean value consistent
 * with the parameter.
 *
 * @param value is the value of the tri-Boolean value
 */
TriBool::TriBool(const TriBool::values value): _value(value) {}

/**
 * @brief A Boolean-to-tri-Boolean constructor
 *
 * This constructor creates a tri-Boolean value consistent
 * with the Boolean parameter.
 *
 * @param value is the Boolean value that the tri-Boolean
 *      object must assume
 */
TriBool::TriBool(const bool value):
    _value(value ? TriBool::TRUE : TriBool::FALSE)
{
}

/**
 * @brief A copy constructor
 *
 * This constructor creates a tri-Boolean value consistent
 * with the parameter.
 *
 * @param orig is the tri-Boolean value that the
 *      object must assume
 */
TriBool::TriBool(const TriBool &orig): _value(orig._value) {}

/**
 * @brief Assignement operator
 *
 * @param orig is the original tri-Boolean object
 * @return a reference to the updated object
 */
TriBool &TriBool::operator=(const TriBool &orig)
{
  _value = orig._value;

  return *this;
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
TriBool operator==(TriBool &&a, const TriBool &b)
{
  if (a.is_uncertain()) {
    return a;
  }

  if (b.is_uncertain()) {
    a._value = TriBool::UNCERTAIN;

    return a;
  }

  a._value = (a._value == b._value ? TriBool::TRUE : TriBool::FALSE);

  return a;
}

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
TriBool operator!=(TriBool &&a, const TriBool &b)
{
  if (a.is_uncertain()) {
    return a;
  }

  if (b.is_uncertain()) {
    a._value = TriBool::UNCERTAIN;

    return a;
  }

  a._value = (a._value != b._value ? TriBool::TRUE : TriBool::FALSE);

  return a;
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
TriBool operator!(TriBool &&a)
{
  switch (a._value) {
  case TriBool::TRUE:
    a._value = TriBool::FALSE;

    return a;
  case TriBool::FALSE:
    a._value = TriBool::TRUE;

    return a;
  case TriBool::UNCERTAIN:
    return a;
  default:
    SAPO_ERROR("Unsupported TriBool value: " << a._value, std::runtime_error);
  }
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
TriBool operator&&(TriBool &&a, const TriBool &b)
{
  if (a.is_true()) {
    a._value = b._value;

    return a;
  }

  if (a.is_uncertain() && b.is_false()) {
    a._value = TriBool::FALSE;
  }

  return a;
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
TriBool operator||(TriBool &&a, const TriBool &b)
{
  switch (a._value) {
  case TriBool::TRUE:
    break;
  case TriBool::FALSE:
    a._value = b._value;
    break;
  case TriBool::UNCERTAIN:
    if (b.is_true()) {
      a._value = b._value;
    }
    break;
  default:
    SAPO_ERROR("Unsupported TriBool value: " << a._value, std::runtime_error);
  }
  return a;
}

/**
 * @brief Print a tri-Boolean in an output stream
 *
 * @param os is the output stream
 * @param tribool is the tri-Boolean value to stream
 * @return a reference to the output stream
 */
std::ostream &operator<<(std::ostream &os, const TriBool &tribool)
{
  if (tribool.is_false()) {
    return (os << "TriBool::FALSE");
  }

  if (tribool.is_true()) {
    return (os << "TriBool::TRUE");
  }

  return (os << "TriBool::UNCERTAIN");
}
