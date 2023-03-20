/**
 * @file FloatingPoints.h
 * @author Alberto Casagrande <acasagrande@units.it>
 * @brief Compute rounded floating point operations
 * @version 0.1
 * @date 2023-02-24
 *
 * @copyright Copyright (c) 2023
 */

#ifndef _FLOATING_POINTS_H_
#define _FLOATING_POINTS_H_

#include <limits>
#include <bitset>
#include <cfenv>

/**
 * @brief Add two floating point values by using the specified rounding
 *
 * This method computes the tightest floating point approximation of
 * the sum of two numbers that respects a specified rounding mode.
 *
 * @tparam T is a floating point type
 * @param a is a floating point value
 * @param b is a floating point value
 * @param rounding is one of the supported rounding mode, i.e.,
 *      `FE_UPWARD`, `FE_DOWNWARD`, `FE_NEAREST`
 * @return the tightest floating point approximation of `a+b` that respects
 *       `rounding`. When `rounding` is `FE_UPWARD`, the returned value
 *       is the least number representable in `T` that is greater that or
 *       equal to `a+b`. When `rounding` is `FE_DOWNWARD`, the returned value
 *       is the greatest number representable in `T` that is lesser that or
 *       equal to `a+b`. When `rounding` is `FE_NEAREST`, the returned value
 *       is a number, representable in `T`, whose distance from `a+b` is
 *       the least one.
 */
template<typename T>
T add(const T &a, const T &b, int rounding);

/**
 * @brief Subtract two floating point values by using the specified rounding
 *
 * This method computes the tightest floating point approximation of the
 * difference between two numbers that respects a specified rounding mode.
 *
 * @tparam T is a floating point type
 * @param a is a floating point value
 * @param b is a floating point value
 * @param rounding is one of the supported rounding mode, i.e.,
 *      `FE_UPWARD`, `FE_DOWNWARD`, `FE_NEAREST`
 * @return the tightest floating point approximation of `a-b` that respects
 *       `rounding`. When `rounding` is `FE_UPWARD`, the returned value
 *       is the least number representable in `T` that is greater that or
 *       equal to `a-b`. When `rounding` is `FE_DOWNWARD`, the returned value
 *       is the greatest number representable in `T` that is lesser that or
 *       equal to `a-b`. When `rounding` is `FE_NEAREST`, the returned value
 *       is a number, representable in `T`, whose distance from `a-b` is
 *       the least one.
 */
template<typename T>
T subtract(const T &a, const T &b, int rounding);

/**
 * @brief Multiply two floating point values by using the specified rounding
 *
 * This method computes the tightest floating point approximation of the
 * product of two numbers that respects a specified rounding mode.
 *
 * @tparam T is a floating point type
 * @param a is a floating point value
 * @param b is a floating point value
 * @param rounding is one of the supported rounding mode, i.e.,
 *      `FE_UPWARD`, `FE_DOWNWARD`, `FE_NEAREST`
 * @return the tightest floating point approximation of `a*b` that respects
 *       `rounding`. When `rounding` is `FE_UPWARD`, the returned value
 *       is the least number representable in `T` that is greater that or
 *       equal to `a*b`. When `rounding` is `FE_DOWNWARD`, the returned value
 *       is the greatest number representable in `T` that is lesser that or
 *       equal to `a*b`. When `rounding` is `FE_NEAREST`, the returned value
 *       is a number, representable in `T`, whose distance from `a*b` is
 *       the least one.
 */
template<typename T>
T multiply(const T &a, const T &b, int rounding);

/**
 * @brief Divide two floating point values by using the specified rounding
 *
 * This method computes the tightest floating point approximation of the
 * quotient of two numbers that respects a specified rounding mode.
 *
 * @tparam T is a floating point type
 * @param a is a floating point value
 * @param b is a floating point value
 * @param rounding is one of the supported rounding mode, i.e.,
 *      `FE_UPWARD`, `FE_DOWNWARD`, `FE_NEAREST`
 * @return the tightest floating point approximation of `a/b` that respects
 *       `rounding`. When `rounding` is `FE_UPWARD`, the returned value
 *       is the least number representable in `T` that is greater that or
 *       equal to `a/b`. When `rounding` is `FE_DOWNWARD`, the returned value
 *       is the greatest number representable in `T` that is lesser that or
 *       equal to `a/b`. When `rounding` is `FE_NEAREST`, the returned value
 *       is a number, representable in `T`, whose distance from `a/b` is
 *       the least one.
 */
template<typename T>
T divide(const T &a, const T &b, int rounding);

template<typename T>
struct IEEE754WorkingType {
};

template<>
struct IEEE754WorkingType<float> {
  using mantissa_type = uint32_t;
  using exponent_type = uint32_t;
  using sign_type = uint32_t;
};

template<>
struct IEEE754WorkingType<double> {
  using mantissa_type = uint64_t;
  using exponent_type = uint32_t;
  using sign_type = uint32_t;
};

// avoid optimizations because of strict aliasing
#pragma GCC push_options
#pragma GCC optimize ("O0")

/**
 * @brief A class to compute rounded floating point operations
 *
 * This class contains low-lever rounded floating point operations. They are
 * not meant to be directly usable.
 *
 * @tparam T is a floating point type
 */
template<typename T,
         typename =
             typename std::enable_if<std::is_floating_point<T>::value>::type>
class IEEE754Rounding
{
public:
  using mantissa_type = typename IEEE754WorkingType<T>::mantissa_type;
  using exponent_type = typename IEEE754WorkingType<T>::exponent_type;
  using sign_type = typename IEEE754WorkingType<T>::exponent_type;

  static constexpr size_t mantissa_size = std::numeric_limits<T>::digits - 1;
  static constexpr size_t exponent_size
      = 8 * sizeof(T) - std::numeric_limits<T>::digits;

  static constexpr exponent_type exponent_bias
      = ((exponent_type(1) << (exponent_size - 1)) - 1);

  /**
   * @brief A union to extract sign, mantissa, exponent bits from IEEE754
   * floating point values
   *
   * Based on https://stackoverflow.com/a/15685301
   */
  union fp_codec {

    fp_codec() {}
    fp_codec(const T &value): value(value) {}

    T value;
    struct __attribute__ ((packed)) {
      mantissa_type mantissa : mantissa_size;
      exponent_type exponent : exponent_size;
      sign_type negative : 1;
    } binary;
  };

private:
  // internal compile-time evaluated constants
  static constexpr mantissa_type implicit_bit
      = (mantissa_type(1) << mantissa_size);
  static constexpr mantissa_type implicit_bit_successor = (implicit_bit << 1);
  static constexpr mantissa_type implicit_mantissa_maximum
      = (implicit_bit_successor - 1);

  static constexpr size_t mantissa_type_size = (8 * sizeof(mantissa_type));

  static constexpr mantissa_type exponent_maximum_value
      = ((1 << exponent_size) - 1);

  static constexpr mantissa_type mantissa_MS_bit
      = (mantissa_type(1) << (mantissa_type_size - 1));

  static constexpr mantissa_type mantissa_mask = (implicit_bit - 1);

  static constexpr size_t half_mantissa_type_size = mantissa_type_size / 2;
  static constexpr mantissa_type low_mantissa_mask
      = (mantissa_type(1) << half_mantissa_type_size) - 1;
  static constexpr mantissa_type high_mantissa_mask
      = (low_mantissa_mask << half_mantissa_type_size);

  /**
   * @brief Return a floating point value fixing its sign
   *
   * @param a is a pointer to a floating point value
   * @param negative is a Boolean flag that specifies the return value sign
   * @return `abs(*(float *)a)` if `negative` is `false` or `-abs(*(float *)a)`
   * otherwise
   */
  static inline T get_value_with_fixed_sign(const fp_codec *a,
                                            const bool &negative)
  {
    if (negative
        == (a->binary.negative
            == 0)) { // `negative` xor the sign of `a` is positive
      return -a->value;
    }

    return a->value;
  }

  /**
   * @brief Add a subnormal value from floating point value
   *
   * @param a is an `fp_codec` pointer to a floating point number whose sign
   *       bit is irrelevant
   * @param a_negative is a Boolean flag for the actual sign of `a`
   * @param subnormal_negative is a Boolean flag for the actual sign of the
   * subnormal value
   * @param rounding is one of the supported rounding mode, i.e.,
   *      `FE_UPWARD`, `FE_DOWNWARD`, `FE_NEAREST`
   * @return the tightest floating point approximation of the sum of `a`
   *    and a subnormal value that respects a specified rounding mode
   */
  static inline T add_subnormal(const fp_codec *a, const bool a_negative,
                                const bool subnormal_negative, int rounding)
  {

    T result{get_value_with_fixed_sign(a, a_negative)};

    fp_codec *r_pointer = (fp_codec *)(&result);

    if ((!a_negative && subnormal_negative && rounding == FE_DOWNWARD)
        || (a_negative && !subnormal_negative && rounding == FE_UPWARD)) {
      r_pointer->binary.mantissa -= 1;

      if (!(r_pointer->binary.mantissa
            & implicit_bit)) {                 // if the mantissa lost
                                               // the implicit bit
        r_pointer->binary.exponent -= 1;       // decrease exponent
        if (r_pointer->binary.exponent == 0) { // if exponent underflow
          r_pointer->binary.mantissa = 0;
          if ((a_negative && subnormal_negative && rounding == FE_DOWNWARD)
              || (!a_negative && !subnormal_negative
                  && rounding == FE_UPWARD)) {
            r_pointer->binary.exponent = 1;
          }
          return result; // return subnormal
        }
        r_pointer->binary.mantissa <<= 1; // shift the mantissa
      }
      return result;
    }

    if ((a_negative && subnormal_negative && rounding == FE_DOWNWARD)
        || (!a_negative && !subnormal_negative && rounding == FE_UPWARD)) {
      if (r_pointer->binary.mantissa
          == implicit_mantissa_maximum) { // if the mantissa already contains
                                          // the largest admissible value

        r_pointer->binary.exponent += 1; // increase exponent
        if (r_pointer->binary.exponent
            == exponent_maximum_value) { // exponent overflow
          r_pointer->binary.mantissa = 0;

          return result;
        }
        r_pointer->binary.mantissa = 0;
      } else {
        ++r_pointer->binary.mantissa;
      }
    }
    return result;
  }

public:
  /**
   * @brief Subtract two floating point values by using the specified rounding
   *
   * This method computes the tightest floating point approximation of the
   * difference between two numbers that respects a specified rounding mode.
   *
   * @tparam T is a floating point type
   * @param a is an `fp_codec` pointer to a floating point number whose sign
   *       bit is irrelevant
   * @param a_negative is a Boolean flag for the actual sign of `a`
   * @param b is an `fp_codec` pointer to a floating point number whose sign
   *       bit is irrelevant
   * @param b_negative is a Boolean flag for the actual sign of `b`
   * @param rounding is one of the supported rounding mode, i.e.,
   *      `FE_UPWARD`, `FE_DOWNWARD`, `FE_NEAREST`
   * @return the tightest floating point approximation of
   *        \f$S=f(a,a_negative)-f(b,b_negative)\f$ that respects `rounding`
   *       where
   *        \f[f(x,x_\textrm{negative}) = \begin{cases}
   *        -|x\rightarrow\textrm{value}| & \textrm{if } x_\textrm{negative}\\
   *        |x| & \textrm{if } \lnot x_\textrm{negative}
   *        \end{cases}.\f]
   *       When `rounding` is `FE_UPWARD`, the returned value
   *       is the least number representable in `T` that is greater that or
   *       equal to \f$S\f$. When `rounding` is `FE_DOWNWARD`, the returned
   *       value is the greatest number representable in `T` that is lesser
   *       that or equal to \f$S\f$. When `rounding` is `FE_NEAREST`, the
   *       returned value is a number, representable in `T`, whose distance
   *       from \f$S\f$ is the least one.
   */
  static T subtract(const fp_codec *a, const bool &a_negative,
                    const fp_codec *b, const bool &b_negative, int rounding)
  {
    if (a_negative) {
      if (!b_negative) {
        return add(a, 1, b, 1, rounding);
      }
    } else {
      if (b_negative) {
        return add(a, 0, b, 0, rounding);
      }
    }

    // handle limit cases
    {
      if (a->binary.exponent == 0) {   // a is either 0 or subnormal
        if (a->binary.mantissa == 0) { // a is 0
          return -get_value_with_fixed_sign(b, b_negative);
        }

        // a is subnormal
        return add_subnormal(b, !b_negative, !a_negative, rounding);
      }

      if (b->binary.exponent == 0) {   // b is either 0 or subnormal
        if (b->binary.mantissa == 0) { // b is 0
          return get_value_with_fixed_sign(a, a_negative);
        }

        // b is subnormal
        return add_subnormal(a, a_negative, b_negative, rounding);
      }

      if (a->binary.exponent == exponent_maximum_value) {
        if (a->binary.mantissa != 0) { // a is nan
          return a->value;
        }

        // abs(a)==infinity

        if (b->binary.exponent == exponent_maximum_value) {
          if (b->binary.mantissa != 0) { // b is nan
            return b->value;
          }

          // a.negative equals b.negative, thus, return nan
          return std::numeric_limits<T>::quiet_NaN();
        }

        // abs(a)==infinity
        return get_value_with_fixed_sign(a, a_negative);
      }

      if (b->binary.exponent == exponent_maximum_value) {

        // a in R and (abs(b)==infinity or b is nan)

        return -get_value_with_fixed_sign(b, b_negative);
      }
    }

    mantissa_type r_mantissa;

    exponent_type exp_delta;
    mantissa_type lost_bits;

    T result;
    fp_codec *r_pointer = (fp_codec *)(&result);

    if (a->binary.exponent > b->binary.exponent) {
      r_mantissa = a->binary.mantissa | implicit_bit;

      exp_delta = a->binary.exponent - b->binary.exponent;
      lost_bits = b->binary.mantissa | implicit_bit;

      r_pointer->value = a->value;
      r_pointer->binary.negative = (a_negative ? 1 : 0);
    } else {
      r_mantissa = b->binary.mantissa | implicit_bit;

      exp_delta = b->binary.exponent - a->binary.exponent;

      if (exp_delta == 0) {
        return (get_value_with_fixed_sign(a, a_negative)
                - get_value_with_fixed_sign(b, b_negative));
      }
      lost_bits = a->binary.mantissa | implicit_bit;

      r_pointer->value = b->value;
      r_pointer->binary.negative = (b_negative ? 0 : 1);
    }

    if (exp_delta < mantissa_size) { // some of the bits in lost_bits
                                     // will overlap r_mantissa
      if (exp_delta == 0) {
        r_mantissa -= lost_bits;
        lost_bits = 0;
      } else {
        r_mantissa -= (lost_bits >> exp_delta);

        lost_bits <<= (mantissa_type_size - exp_delta);
        if (lost_bits) {
          lost_bits = mantissa_type(0) - lost_bits;
        }
      }
      while (!(r_mantissa & implicit_bit)) {   // mantissa lost implicit bit
        if (r_pointer->binary.exponent == 0) { // exponent underflow
          r_pointer->binary.mantissa = 0;
          if ((r_pointer->binary.negative == 0 && rounding == FE_UPWARD)
              || (r_pointer->binary.negative == 1
                  && rounding == FE_DOWNWARD)) {
            r_pointer->binary.exponent = 1;
          }
          return result;
        }
        r_pointer->binary.exponent -= 1; // decrease exponent
        r_mantissa = ((r_mantissa << 1)  // shift r_mantissa by one
                      | ((mantissa_MS_bit & lost_bits) == 0
                             ? 0
                             : 1)); // and recover one bit
                                    // from lost_bits
        lost_bits <<= 1;
      }
    }

    if (lost_bits != 0) { // some of the lost bits are not 0s and the
                          // candidate answer is larger than the absolute
                          // value of the real sum among a and b
      if ((r_pointer->binary.negative && rounding == FE_UPWARD)
          || (!r_pointer->binary.negative && rounding == FE_DOWNWARD)) {

        --r_mantissa;                            // reduce the mantissa
        if (!(r_mantissa & implicit_bit)) {      // if the mantissa lost
                                                 // the implicit bit
          r_pointer->binary.exponent -= 1;       // decrease exponent
          if (r_pointer->binary.exponent == 0) { // if exponent underflow
            r_pointer->binary.mantissa = 1;

            return result; // return subnormal
          }
          r_mantissa <<= 1; // shift the mantissa
        }
      }
    }
    r_pointer->binary.mantissa = r_mantissa;

    return result;
  }

  /**
   * @brief Add two floating point values by using the specified rounding
   *
   * This method computes the tightest floating point approximation of the
   * sum of two numbers that respects a specified rounding mode.
   *
   * @tparam T is a floating point type
   * @param a is an `fp_codec` pointer to a floating point number whose sign
   *       bit is irrelevant
   * @param a_negative is a Boolean flag for the actual sign of `a`
   * @param b is an `fp_codec` pointer to a floating point number whose sign
   *       bit is irrelevant
   * @param b_negative is a Boolean flag for the actual sign of `b`
   * @param rounding is one of the supported rounding mode, i.e.,
   *      `FE_UPWARD`, `FE_DOWNWARD`, `FE_NEAREST`
   * @return the tightest floating point approximation of
   *        \f$S=f(a,a_negative)+f(b,b_negative)\f$ that respects `rounding`
   *       where
   *        \f[f(x,x_\textrm{negative}) = \begin{cases}
   *        -|x\rightarrow\textrm{value}| & \textrm{if } x_\textrm{negative}\\
   *        |x| & \textrm{if } \lnot x_\textrm{negative}
   *        \end{cases}.\f]
   *       When `rounding` is `FE_UPWARD`, the returned value
   *       is the least number representable in `T` that is greater that or
   *       equal to \f$S\f$. When `rounding` is `FE_DOWNWARD`, the returned
   *       value is the greatest number representable in `T` that is lesser
   *       that or equal to \f$S\f$. When `rounding` is `FE_NEAREST`, the
   *       returned value is a number, representable in `T`, whose distance
   *       from \f$S\f$ is the least one.
   */
  static T add(const fp_codec *a, const bool &a_negative, const fp_codec *b,
               const bool &b_negative, int rounding)
  {
    if (a_negative) {
      if (!b_negative) {
        return subtract(b, 0, a, 0, rounding);
      }
    } else {
      if (b_negative) {
        return subtract(a, 0, b, 0, rounding);
      }
    }

    // handle limit cases
    {
      if (a->binary.exponent == 0) {   // a is either 0 or subnormal
        if (a->binary.mantissa == 0) { // a is 0
          return get_value_with_fixed_sign(b, b_negative);
        }

        // a is subnormal
        return add_subnormal(b, b_negative, a_negative, rounding);
      }

      if (b->binary.exponent == 0) {   // b is either 0 or subnormal
        if (b->binary.mantissa == 0) { // b is 0
          return get_value_with_fixed_sign(a, a_negative);
        }

        // b is subnormal
        return add_subnormal(a, a_negative, b_negative, rounding);
      }

      if (a->binary.exponent == exponent_maximum_value) {
        if (a->binary.mantissa != 0) { // a is nan
          return a->value;
        }

        // abs(a)==infinity

        if (b->binary.exponent == exponent_maximum_value) {
          if (b->binary.mantissa != 0) { // b is nan
            return b->value;
          }

          // a.negative equals b.negative, thus, return a
        }

        // abs(a)==infinity
        return get_value_with_fixed_sign(a, a_negative);
      }

      if (b->binary.exponent == exponent_maximum_value) {

        // a in R and (abs(b)==infinity or b is nan)
        return get_value_with_fixed_sign(b, b_negative);
      }
    }

    mantissa_type r_mantissa;

    exponent_type exp_delta;
    mantissa_type lost_bits;

    T result;
    fp_codec *r_pointer = (fp_codec *)(&result);

    if (a->binary.exponent > b->binary.exponent) {
      r_mantissa = a->binary.mantissa | implicit_bit;

      exp_delta = a->binary.exponent - b->binary.exponent;
      lost_bits = b->binary.mantissa | implicit_bit;

      r_pointer->value = a->value;
      r_pointer->binary.negative = (a_negative ? 1 : 0);
    } else {
      r_mantissa = b->binary.mantissa | implicit_bit;

      exp_delta = b->binary.exponent - a->binary.exponent;
      lost_bits = a->binary.mantissa | implicit_bit;

      r_pointer->value = b->value;
      r_pointer->binary.negative = (b_negative ? 1 : 0);
    }

    if (exp_delta < mantissa_size) {
      if (exp_delta == 0) {
        r_mantissa += lost_bits;
        lost_bits = 0;
      } else {
        r_mantissa += (lost_bits >> exp_delta);

        lost_bits <<= (mantissa_type_size - exp_delta);
      }
      if (r_mantissa
          & implicit_bit_successor) {    // mantissa overflow implicit bit
        r_pointer->binary.exponent += 1; // increase exponent
        if (r_pointer->binary.exponent
            == exponent_maximum_value) { // exponent overflow
          r_pointer->binary.mantissa = 0;

          return result;
        }
        lost_bits
            += r_mantissa
               | 1; // add the least significant bit of r_mantissa to lost_bit
        r_mantissa >>= 1; // shift r_mantissa by one
      }
    }

    if (lost_bits != 0) { // some of the lost bits are not 0s and the
                          // candidate answer is lesser than the absolute
                          // value of the real sum among a and b
      if ((!r_pointer->binary.negative && rounding == FE_UPWARD)
          || (r_pointer->binary.negative && rounding == FE_DOWNWARD)) {

        if (r_mantissa == implicit_mantissa_maximum) { // if the mantissa
                                                       // already contains the
                                                       // largest admissible
                                                       // value
          r_pointer->binary.exponent += 1;             // increase exponent
          if (r_pointer->binary.exponent
              == exponent_maximum_value) { // exponent overflow
            r_pointer->binary.mantissa = 0;

            return result;
          }
          r_mantissa = 0;
        } else {
          ++r_mantissa;
        }
      }
    }

    r_pointer->binary.mantissa = r_mantissa;

    return result;
  }

  static inline mantissa_type high_part(const mantissa_type &mantissa)
  {
    return mantissa >> half_mantissa_type_size;
  }

  static inline mantissa_type low_part(const mantissa_type &mantissa)
  {
    return (mantissa << half_mantissa_type_size) >> half_mantissa_type_size;
  }

  /**
   * @brief Multiply two floating point values by using the specified rounding
   *
   * This method computes the tightest floating point approximation of the
   * sum of two numbers that respects a specified rounding mode. This goal is
   * achieved by implementing the Karatsuba's algorithm.
   *
   * @tparam T is a floating point type
   * @param a is an `fp_codec` pointer to a floating point number
   * @param b is an `fp_codec` pointer to a floating point number
   * @param rounding is one of the supported rounding mode, i.e.,
   *      `FE_UPWARD`, `FE_DOWNWARD`, `FE_NEAREST`
   * @return the tightest floating point approximation of \f$a*b\f$ that
   *       respects `rounding`.
   *       When `rounding` is `FE_UPWARD`, the returned value is the least
   *       number representable in `T` that is greater that or equal to
   *       \f$a*b\f$. When `rounding` is `FE_DOWNWARD`, the returned value
   *       is the greatest number representable in `T` that is lesser that
   *       or equal to \f$a*b\f$. When `rounding` is `FE_NEAREST`, the
   *       returned value is a number, representable in `T`, whose distance
   *       from \f$a*b\f$ is the least one.
   */
  static T multiply(const fp_codec *a, const fp_codec *b, int rounding)
  {
    if ((a->value == 0) || (b->value == 0)) {
      return 0;
    }

    bool negative = a->binary.negative ^ b->binary.negative;
    if (a->binary.exponent == exponent_maximum_value) {
      // either std::abs(a) == inf or a == nan
      return (negative ? -a->value : a->value);
    }
    if (b->binary.exponent == exponent_maximum_value) {
      // either std::abs(b) == inf or b == nan
      return (negative ? -b->value : b->value);
    }

    exponent_type new_exponent = a->binary.exponent + b->binary.exponent;

    if (new_exponent <= exponent_bias) {
      if (!negative && rounding == FE_UPWARD) {
        return std::numeric_limits<T>::min();
      }
      if (negative && rounding == FE_DOWNWARD) {
        return -std::numeric_limits<T>::min();
      }

      return 0;
    }

    new_exponent -= exponent_bias;
    if (new_exponent >= exponent_maximum_value) { // overflow
      return (negative ? -std::numeric_limits<T>::infinity()
                       : std::numeric_limits<T>::infinity());
    }

    const mantissa_type A = a->binary.mantissa | implicit_bit;
    const mantissa_type B = b->binary.mantissa | implicit_bit;

    const mantissa_type a_high{high_part(A)}, a_low{low_part(A)},
        b_high{high_part(B)}, b_low{low_part(B)};

    mantissa_type H = a_high * b_high;
    mantissa_type L = a_low * b_low;
    mantissa_type M = (a_high + a_low) * (b_high + b_low) - H - L;

    H += high_part(M);
    M = low_part(M) + high_part(L);
    H += high_part(M);
    M = low_part(M);
    L = low_part(L);

    constexpr size_t H_shift = (mantissa_type_size - (mantissa_size));

    H = (H << H_shift) | (M >> (half_mantissa_type_size - H_shift));
    M = (M << (half_mantissa_type_size + H_shift)) | L << 1;

    if (H & implicit_bit_successor) {
      M = ((H | 1) << (mantissa_type_size - 1)) | (M > 1);
      H >>= 1;

      new_exponent++;

      if (new_exponent >= exponent_maximum_value) { // exponent overflow
        return (negative ? -std::numeric_limits<T>::infinity()
                         : std::numeric_limits<T>::infinity());
      }
    }

    if (M != 0) {
      if ((!negative && rounding == FE_UPWARD)
          || (negative && rounding == FE_DOWNWARD)) {
        H++;
        if (H & implicit_bit_successor) {
          H >>= 1;
          new_exponent++;

          if (new_exponent >= exponent_maximum_value) { // exponent overflow
            return (negative ? -std::numeric_limits<T>::infinity()
                             : std::numeric_limits<T>::infinity());
          }
        }
      }
    }

    T result{0};
    fp_codec *r_pointer = (fp_codec *)(&result);

    r_pointer->binary.negative = (negative ? 1 : 0);
    r_pointer->binary.exponent = new_exponent;
    r_pointer->binary.mantissa = H;

    return result;
  }
};

template<typename T>
inline T add(const T &a, const T &b, int rounding)
{
  using fp_codec = typename IEEE754Rounding<T>::fp_codec;
  const fp_codec *fp_a = reinterpret_cast<const fp_codec *>(&a);
  const fp_codec *fp_b = reinterpret_cast<const fp_codec *>(&b);

  return IEEE754Rounding<T>::add(fp_a, fp_a->binary.negative,
                                 fp_b, fp_b->binary.negative,
      rounding);
}

template<typename T>
inline T subtract(const T &a, const T &b, int rounding)
{
  using fp_codec = typename IEEE754Rounding<T>::fp_codec;

  const fp_codec *fp_a = reinterpret_cast<const fp_codec *>(&a);
  const fp_codec *fp_b = reinterpret_cast<const fp_codec *>(&b);

  return IEEE754Rounding<T>::subtract(fp_a, fp_a->binary.negative,
                                      fp_b, fp_b->binary.negative,
      rounding);
}

template<typename T>
inline T multiply(const T &a, const T &b, int rounding)
{
  using fp_codec = typename IEEE754Rounding<T>::fp_codec;

  const fp_codec *fp_a = reinterpret_cast<const fp_codec *>(&a);
  const fp_codec *fp_b = reinterpret_cast<const fp_codec *>(&b);

  return IEEE754Rounding<T>::multiply(fp_a, fp_b, rounding);
}

template<typename T>
std::ostream &operator<<(std::ostream &os,
                         const typename IEEE754Rounding<T>::fp_codec &a)
{
  os << "value: " << a.value << std::endl
     << "  negative: " << std::bitset<1>(a.binary.negative) << std::endl
     << "  exponent: "
     << std::bitset<IEEE754Rounding<T>::exponent_size>(a.binary.exponent)
     << std::endl
     << "  mantissa: "
     << std::bitset<IEEE754Rounding<T>::mantissa_size>(a.binary.mantissa)
     << std::endl;

  return os;
}

#pragma GCC pop_options

#endif // _FLOATING_POINTS_H_