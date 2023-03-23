/**
 * @file SymbolicAlgebra.h
 * @author Alberto Casagrande <acasagrande@units.it>
 * @brief A simple symbolic algebra module
 * @version 0.1
 * @date 2021-12-03
 *
 * @copyright Copyright (c) 2021-2022
 */

#ifndef SYMBOLIC_ALGEBRA_H_
#define SYMBOLIC_ALGEBRA_H_

#include <iostream>
#include <list>
#include <map>
#include <set>
#include <vector>
#include <iterator>
#include <type_traits>

#ifdef WITH_THREADS
#include <mutex>

#endif // WITH_THREADS

#include "ErrorHandling.h"

/*!
 *  \addtogroup SymbolicAlgebra
 *  @{
 */

//! Symbolic algebra namespace
namespace SymbolicAlgebra
{

/**
 * @brief Expression types
 */
typedef enum {
  CONSTANT,   // Constant
  SYMBOL,     // Symbol
  FINITE_SUM, // Finite sum
  FINITE_PROD // Finite product
} ExpressionType;

namespace low_level
{
template<typename C>
class constant_type;

template<typename C>
class symbol_type;

template<typename C>
class finite_sum_type;

template<typename C>
class finite_prod_type;

template<typename C>
class base_expression_type;

}
template<typename C>
class Expression;

/**
 * @brief A class to represent algebraic symbols
 *
 * This class represents algebraic symbols. It uniquely
 * associates any symbol name to an identifier.
 *
 * @tparam C is the type of numeric constants
 */
template<typename C = double>
class Symbol : public Expression<C>
{
public:
  typedef unsigned int SymbolIdType; //!< The type of symbol identificators

#ifdef WITH_THREADS

  static std::mutex _mutex;

#endif // WITH_THREADS

protected:
  static std::map<std::string, SymbolIdType>
      _declared_symbols; //!< Name to symbol identificator map
  static std::vector<std::string>
      _symbol_names; //!< Identificator to symbol name map

public:
  /**
   * @brief Build a new empty Symbol
   */
  Symbol();

  /**
   * @brief Build a new symbol object
   *
   * If the symbol name has never been used, then this method reserves a new id
   * for the name and updates both the map of the declared symbols, that
   * associates a symbol name to its id, and the symbol name map, that relates
   * a symbol id to its name.
   *
   * @param name is the symbol name
   */
  Symbol(const char *name);

  /**
   * @brief Build a new symbol object
   *
   * If the symbol name has never been used, then this method reserves a new id
   * for the name and updates both the map of the declared symbols, that
   * associates a symbol name to its id, and the symbol name map, that relates
   * a symbol id to its name.
   *
   * @param name is the symbol name
   */
  Symbol(const std::string &name);

  /**
   * @brief Copy constructor
   *
   * @param orig is the Symbol object to be copied
   */
  Symbol(const Symbol<C> &orig);

  /**
   * @brief Assign a symbol
   *
   * @param orig is the symbol to be assigned
   * @return A reference to the updated object
   */
  Symbol<C> &operator=(const Symbol<C> &orig);

  /**
   * @brief Assign a symbol
   *
   * @param orig is the symbol to be assigned
   * @return A reference to the updated object
   */
  Symbol<C> &operator=(Symbol<C> &&orig);

  /**
   * @brief Get the Symbol name
   *
   * @return the symbol name
   */
  const std::string &get_name() const
  {
    return Symbol<C>::get_symbol_name(get_id());
  }

  /**
   * @brief Test whether two symbols are the same
   *
   * @param a is a symbol
   * @return `true` if and only if `a` and `*this`
   *    represent the same symbol
   *
   */
  inline bool operator==(const Symbol<C> &a) const
  {
    return get_id() == a.get_id();
  }

  /**
   * @brief Get the Symbol id
   *
   * @return The id of the represented symbol
   */
  const SymbolIdType &get_id() const
  {
    return ((low_level::symbol_type<C> *)(this->_ex))->get_id();
  }

  /**
   * @brief Get the name associated to a symbol id
   *
   * @param id is the query symbol id
   * @return the name corresponding to the query symbol id
   */
  static const std::string &get_symbol_name(const SymbolIdType &id)
  {
    return _symbol_names[id];
  }
};

/**
 * @brief Test whether two symbols are the same
 *
 * @tparam C is the type of constants
 * @param a is a symbol
 * @param b is a symbol name
 * @return `true` if and only if `a` and `b`
 *    represent the same symbol
 *
 */
template<typename C>
inline bool operator==(const Symbol<C> &a, const std::string &b)
{
  return a == Symbol<C>(b);
}

/**
 * @brief Test whether two symbols are the same
 *
 * @tparam C is the type of constants
 * @param a is a symbol name
 * @param b is a symbol
 * @return `true` if and only if `a` and `b`
 *    represent the same symbol
 *
 */
template<typename C>
inline bool operator==(const std::string &a, const Symbol<C> &b)
{
  return b == Symbol<C>(a);
}
}

/**
 * @brief A specialization of the std::less structure
 *
 * This specialization of the std::less structure is meant
 * to compare symbols.
 *
 * @tparam C is the constant value type
 */
template<typename C>
struct std::less<SymbolicAlgebra::Symbol<C>> {

  /**
   * @brief Compare two symbols
   *
   * This method is used to distinguish and artificially order
   * two symbols in containers such as sets. It is not meant to
   * semantically compare two symbols.
   *
   * @param a is a symbol
   * @param b is a symbol
   * @return `true` if and only if the index of `a` is smaller, than
   *         that of `b`.
   */
  constexpr bool operator()(const SymbolicAlgebra::Symbol<C> &a,
                            const SymbolicAlgebra::Symbol<C> &b) const
  {
    return a.get_id() < b.get_id();
  }
};

namespace SymbolicAlgebra
{

#ifdef WITH_THREADS
template<typename C>
std::mutex Symbol<C>::_mutex;
#endif // WITH_THREADS

template<typename C>
std::map<std::string, typename Symbol<C>::SymbolIdType>
    Symbol<C>::_declared_symbols;
template<typename C>
std::vector<std::string> Symbol<C>::_symbol_names;

/**
 * @brief A class to represent algebraic symbolic expressions
 *
 * This class is a wrapper class: the expression representation is achieved by
 * the `base_expression_type` class and its hierarchy.
 *
 * @tparam C is the type of numeric constants
 */
template<typename C = double>
class Expression
{
protected:
  low_level::base_expression_type<C>
      *_ex; //!< A pointer to a base expression object.

  /**
   * @brief Build a new Expression object
   *
   * @param _ex is the expression representation
   */
  Expression(low_level::base_expression_type<C> *_ex);

public:
  typedef std::map<Symbol<C>, Expression<C>>
      replacement_type; //!< Replacement type
  typedef std::map<Symbol<C>, C>
      interpretation_type; //!< Symbol interpretation type

  /**
   * @brief Build an empty Expression object
   */
  Expression();

  /**
   * @brief Build a constant expression
   *
   * @param value is the value of the expression
   */
  Expression(const int value);

  /**
   * @brief Build a constant expression
   *
   * @param value is the value of the expression
   */
  template<typename T = C, std::enable_if_t<!std::is_same_v<T, double>>>
  Expression(const double value):
      _ex(new low_level::constant_type<C>(static_cast<C>(value)))
  {
  }

  /**
   * @brief Build a constant expression
   *
   * @param value is the value of the expression
   */
  Expression(const C value);

  /**
   * @brief Copy constructor
   *
   * @param orig is the model for the new object
   */
  Expression(const Expression<C> &orig);

  /**
   * @brief Replace symbol occurrences by using expressions
   *
   * This method replaces any occurrence of the symbols in the domain
   * of `replacement` with the corresponding expressions. This method
   * updates this object.
   *
   * @param replacements associates symbols to their replacements
   * @return A reference to the updated object
   */
  Expression<C> &replace(const replacement_type &replacements);

  /**
   * @brief Turn the expression into a sum of products
   *
   * This method turns this expression into an algebraic equivalent
   * sum of products.
   *
   * @return A reference to the updated object
   */
  Expression<C> &expand();

  /**
   * @brief Numerically evaluate the expression
   *
   * @return The numeric evaluation of the expression
   */
  inline C evaluate() const
  {
    const C value = _ex->evaluate();

    return value == 0 ? 0 : value;
  }

  /**
   * @brief Apply a symbolic interpretation to an expression
   *
   * This method evaluates the expression using the symbol interpretation
   * provided as parameter.
   *
   * @param interpretation is an interpretation for the symbols
   * @return the expression obtained from the current object
   *         and `interpretation`
   */
  Expression<C> apply(const interpretation_type &interpretation) const;

  /**
   * @brief Get the rational polynomial representation of an expression
   *
   * @return the rational polynomial representation of this expression
   */
  Expression<C> get_rational_form() const;

  /**
   * @brief Establish whether is a polynomial expression
   *
   * @return `true` if and only if the expression is a polynomial
   */
  bool is_a_polynomial() const;

  /**
   * @brief Establish whether is in rational polynomial form
   *
   * @return `true` if and only if the expression is either
   *        polynomial or the ratio between two polynomials
   */
  bool is_in_rational_form() const;

  /**
   * @brief Establish whether this expression is a constant
   *
   * @return `true` if and only if the expression is a constant
   */
  bool is_a_constant() const;

  /**
   * @brief Get the numerator of the expression
   *
   * @return the numerator of the expression
   */
  Expression<C> get_numerator() const;

  /**
   * @brief Get the denominator of the expression
   *
   * @return the denominator of the expression
   */
  Expression<C> get_denominator() const;

  /**
   * @brief Get the derivative of the expression
   *
   * @param symbol is the symbol with respect to the derivative should be
   * computed
   * @return the derivative of the current expression with respect to
   *        `symbol`
   */
  Expression<C> get_derivative_wrt(const Symbol<C> &symbol) const;

  /**
   * @brief Get the coefficient of a term having a given degree
   *
   * @param symbol is the symbol whose term coefficient is aimed
   * @param degree is the degree of the aimed term
   * @return The coefficient of the term \f$\texrm{symbol}^\textrm{degree}\f$
   */
  inline Expression<C> get_coeff(const Symbol<C> &symbol,
                                 const int degree) const
  {
    return Expression<C>(_ex->get_coeff(symbol.get_id(), degree));
  }

  /**
   * @brief Get the coefficients of the polynomial expression
   *
   * @param symbol is the symbol whose term coefficient is aimed
   * @return The coefficients of the polynomial expression on `symbol`
   */
  std::map<int, Expression<C>> get_coeffs(const Symbol<C> &symbol) const
  {
    std::map<int, Expression<C>> coeffs;

    auto base_coeffs = _ex->get_coeffs(symbol.get_id());
    for (auto it = std::begin(base_coeffs); it != std::end(base_coeffs);
         ++it) {
      coeffs[it->first] = Expression<C>(it->second);
    }

    return coeffs;
  }

  /**
   * @brief Get the degree of a symbol
   *
   * @param symbol is the symbol whose degree is aimed
   * @return the degree of the parameter in the expression
   */
  inline int degree(const Symbol<C> &symbol) const
  {
    return _ex->degree(
        ((const low_level::symbol_type<C> *)(symbol._ex))->get_id());
  }

  /**
   * @brief Get the set of the expression symbols
   *
   * @return the set of the expression symbols
   */
  std::set<Symbol<C>> get_symbols() const
  {
    std::set<Symbol<C>> result;

    for (typename Symbol<C>::SymbolIdType sid: _ex->get_symbol_ids()) {
      result.emplace(Symbol<C>::get_symbol_name(sid));
    }

    return result;
  }

  /**
   * @brief Get a subset of the expression symbols
   *
   * @param max_cardinality is the maximum number of symbols to be returned
   * @return A set, whose cardinality is `max_cardinality` at most, of the
   * symbols in the expression
   */
  std::set<Symbol<C>> get_symbols(const unsigned int max_cardinality) const
  {
    std::set<Symbol<C>> result;

    for (typename Symbol<C>::SymbolIdType sid:
         _ex->get_symbol_ids(max_cardinality)) {
      result.emplace(Symbol<C>::get_symbol_name(sid));
    }

    return result;
  }

  /**
   * @brief Get the constant term of the expression
   *
   * @return The constant term of the expression
   */
  Expression<C> get_constant_term() const
  {
    Expression<C>::replacement_type replacements;
    for (Symbol<C> symbol: this->get_symbols()) {
      replacements[symbol] = 0;
    }

    return Expression<C>(*this).replace(replacements);
  }

  /**
   * @brief Establish whether the expression contains symbols
   *
   * @return `true` if and only if the expression has contains
   *         one symbol at least
   */
  bool has_symbols() const
  {
    return _ex->has_symbols();
  }

  /**
   * @brief Destroy the Expression object
   */
  ~Expression();

  /**
   * @brief Assign an expression
   *
   * @param rhs is the expression to be assigned
   * @return A reference to the updated object
   */
  Expression<C> &operator=(const Expression<C> &rhs);

  /**
   * @brief Assign an expression
   *
   * @param rhs is the expression to be assigned
   * @return A reference to the updated object
   */
  Expression<C> &operator=(Expression<C> &&rhs);

  /**
   * @brief Add an expression to the current one
   *
   * @param rhs is the expression to be added
   * @return A reference to the updated object
   */
  const Expression<C> &operator+=(const Expression<C> &rhs);

  /**
   * @brief Subtract an expression to the current one
   *
   * @param rhs is the expression to be subtracted
   * @return A reference to the updated object
   */
  const Expression<C> &operator-=(const Expression<C> &rhs);

  /**
   * @brief Multiply an expression to the current one
   *
   * @param rhs is the expression to be multiplied
   * @return A reference to the updated object
   */
  const Expression<C> &operator*=(const Expression<C> &rhs);

  /**
   * @brief Divide an expression to the current one
   *
   * @param rhs is the expression to be divided
   * @return A reference to the updated object
   */
  const Expression<C> &operator/=(const Expression<C> &rhs);

  /**
   * @brief Add an expression to the current one
   *
   * @param rhs is the expression to be added
   * @return A reference to the updated object
   */
  const Expression<C> &operator+=(Expression<C> &&rhs);

  /**
   * @brief Subtract an expression to the current one
   *
   * @param rhs is the expression to be subtracted
   * @return A reference to the updated object
   */
  const Expression<C> &operator-=(Expression<C> &&rhs);

  /**
   * @brief Multiply an expression to the current one
   *
   * @param rhs is the expression to be multiplied
   * @return A reference to the updated object
   */
  const Expression<C> &operator*=(Expression<C> &&rhs);

  /**
   * @brief Divide an expression to the current one
   *
   * @param rhs is the expression to be divided
   * @return A reference to the updated object
   */
  const Expression<C> &operator/=(Expression<C> &&rhs);

  /**
   * @brief Add a numeric value to the current expression
   *
   * @tparam T is the type of the numeric value
   * @param value is the numeric value to be added
   * @return A reference to the updated object
   */
  template<typename T,
           std::enable_if_t<std::is_arithmetic_v<T>>>
  const Expression<C> &operator+=(const T value)
  {
    _ex = _ex->add(static_cast<C>(value));

    return *this;
  }

  /**
   * @brief Subtract a numeric value from the current expression
   *
   * @tparam T is the type of the numeric value
   * @param value is the numeric value to be subtracted
   * @return A reference to the updated object
   */
  template<typename T,
           std::enable_if_t<std::is_arithmetic_v<T>>>
  const Expression<C> &operator-=(const T value)
  {
    _ex = _ex->subtract(static_cast<C>(value));

    return *this;
  }

  /**
   * @brief Multiply a numeric value and the current expression
   *
   * @tparam T is the type of the numeric value
   * @param value is the numeric value to be multiplied
   * @return A reference to the updated object
   */
  template<typename T,
           std::enable_if_t<std::is_arithmetic_v<T>>>
  const Expression<C> &operator*=(const T value)
  {
    _ex = _ex->multiply(static_cast<C>(value));

    return *this;
  }

  /**
   * @brief Divide an expression to the current one
   *
   * @tparam T is the type of the numeric value
   * @param value is the numeric divisor
   * @return A reference to the updated object
   */
  template<typename T,
           std::enable_if_t<std::is_arithmetic_v<T>>>
  const Expression<C> &operator/=(const T value)
  {
    _ex = _ex->divided_by(static_cast<C>(value));

    return *this;
  }

  /**
   * @brief Complement the expression
   *
   * @return the complement expression
   */
  Expression<C> operator-() const
  {
    return Expression<C>((_ex->clone())->complement());
  }

  /**
   * @brief Check whether an expression is semantically equivalent to a
   * constant value
   *
   * @tparam C is the numeric type of constants
   * @param lhs is an expression
   * @param rhs is a constant value
   * @return `true` if and only if the expression is equivalent to the
   *         constant value
   */
  template<typename T>
  friend bool operator==(const Expression<T> &lhs, const T rhs);

  /**
   * @brief Get the sum between two expressions
   *
   * This method builds an expression that represents the sum
   * between two expressions.
   *
   * @tparam T is the type of the constant values in the expression
   * @param lhs is an expression
   * @param rhs is another expression
   * @return the expression that represents `lhs`+`rhs`
   */
  template<typename T>
  friend Expression<T> operator+(const Expression<T> &lhs,
                                 const Expression<T> &rhs);

  /**
   * @brief Get the subtraction between two expressions
   *
   * This method builds an expression that represents the subtractions
   * between two expressions.
   *
   * @tparam T is the type of the constant values in the expression
   * @param lhs is an expression
   * @param rhs is another expression
   * @return the expression that represents `lhs`-`rhs`
   */
  template<typename T>
  friend Expression<T> operator-(const Expression<T> &lhs,
                                 const Expression<T> &rhs);

  /**
   * @brief Get the division between two expressions
   *
   * This method builds an expression that represents the division between
   * two expressions.
   *
   * @tparam T is the type of the constant values in the expression
   * @param lhs is the dividend expression
   * @param rhs is the divisor expression
   * @return the expression that represents `lhs`/`rhs`
   */
  template<typename T>
  friend Expression<T> operator*(const Expression<T> &lhs,
                                 const Expression<T> &rhs);

  /**
   * @brief Get the division between two expressions
   *
   * This method builds an expression that represents the division between
   * two expressions.
   *
   * @tparam T is the type of the constant values in the expression
   * @param lhs is the dividend expression
   * @param rhs is the divisor expression
   * @return the expression that represents `lhs`/`rhs`
   */
  template<typename T>
  friend Expression<T> operator/(const Expression<T> &lhs,
                                 const Expression<T> &rhs);

  /**
   * @brief Get the sum between two expressions
   *
   * This method builds an expression that represents the sum
   * between two expressions.
   *
   * @tparam T is the type of the constant values in the expression
   * @param lhs is an expression
   * @param rhs is another expression
   * @return the expression that represents `lhs`+`rhs`
   */
  template<typename T>
  friend Expression<T> operator+(Expression<T> &&lhs,
                                 const Expression<T> &rhs);

  /**
   * @brief Get the subtraction between two expressions
   *
   * This method builds an expression that represents the subtractions
   * between two expressions.
   *
   * @tparam T is the type of the constant values in the expression
   * @param lhs is an expression
   * @param rhs is another expression
   * @return the expression that represents `lhs`-`rhs`
   */
  template<typename T>
  friend Expression<T> operator-(Expression<T> &&lhs,
                                 const Expression<T> &rhs);

  /**
   * @brief Get the multiplication between two expressions
   *
   * This method builds an expression that represents the multiplication
   * between two expressions.
   *
   * @tparam T is the type of the constant values in the expression
   * @param lhs is an expression
   * @param rhs is another expression
   * @return the expression that represents `lhs`*`rhs`
   */
  template<typename T>
  friend Expression<T> operator*(Expression<T> &&lhs,
                                 const Expression<T> &rhs);

  /**
   * @brief Get the division between two expressions
   *
   * This method builds an expression that represents the division between
   * two expressions.
   *
   * @tparam T is the type of the constant values in the expression
   * @param lhs is the dividend expression
   * @param rhs is the divisor expression
   * @return the expression that represents `lhs`/`rhs`
   */
  template<typename T>
  friend Expression<T> operator/(Expression<T> &&lhs,
                                 const Expression<T> &rhs);

  /**
   * @brief Get the sum between two expressions
   *
   * This method builds an expression that represents the sum
   * between two expressions.
   *
   * @tparam T is the type of the constant values in the expression
   * @param lhs is an expression
   * @param rhs is another expression
   * @return the expression that represents `lhs`+`rhs`
   */
  template<typename T>
  friend Expression<T> operator+(const Expression<T> &lhs,
                                 Expression<T> &&rhs);

  /**
   * @brief Get the subtraction between two expressions
   *
   * This method builds an expression that represents the subtractions
   * between two expressions.
   *
   * @tparam T is the type of the constant values in the expression
   * @param lhs is an expression
   * @param rhs is another expression
   * @return the expression that represents `lhs`-`rhs`
   */
  template<typename T>
  friend Expression<T> operator-(const Expression<T> &lhs,
                                 Expression<T> &&rhs);

  /**
   * @brief Get the multiplication between two expressions
   *
   * This method builds an expression that represents the multiplication
   * between two expressions.
   *
   * @tparam T is the type of the constant values in the expression
   * @param lhs is an expression
   * @param rhs is another expression
   * @return the expression that represents `lhs`*`rhs`
   */
  template<typename T>
  friend Expression<T> operator*(const Expression<T> &lhs,
                                 Expression<T> &&rhs);

  /**
   * @brief Get the division between two expressions
   *
   * This method builds an expression that represents the division between
   * two expressions.
   *
   * @tparam T is the type of the constant values in the expression
   * @param lhs is the dividend expression
   * @param rhs is the divisor expression
   * @return the expression that represents `lhs`/`rhs`
   */
  template<typename T>
  friend Expression<T> operator/(const Expression<T> &lhs,
                                 Expression<T> &&rhs);

  /**
   * @brief Sum an expression and a numeric value
   *
   * This method builds an expression that represents the sum of the two
   * expressions passed as parameters.
   *
   * @tparam T is the type of numeric value
   * @param lhs is an expression
   * @param value is an expression
   * @return An expression that represents the sum `lhs`+`value`
   */
  template<typename T,
           typename
           = typename std::enable_if_t<std::is_arithmetic_v<T>, T>>
  friend Expression<C> operator+(const Expression<C> &lhs, const T value)
  {
    if (value == 0) {
      return lhs;
    }

    low_level::base_expression_type<C> *lhs_ex = lhs._ex->clone();

    return Expression<C>(lhs_ex->add(static_cast<C>(value)));
  }

  /**
   * @brief Subtract an expression from a numeric value
   *
   * This method builds an expression that represents the subtraction
   * of an expression from a numeric value.
   *
   * @tparam T is the type of numeric value
   * @param value is an expression
   * @param rhs is an expression
   * @return An expression that represents the subtraction `value`-`rhs`
   */
  template<typename T,
           typename
           = typename std::enable_if_t<std::is_arithmetic_v<T>, T>>
  friend inline Expression<C> operator-(const T value,
                                        const Expression<C> &rhs)
  {
    return ((-rhs) + value);
  }

  /**
   * @brief Multiply an expression and a numeric value
   *
   * This method builds an expression that represents the product
   * of an expression and a numeric value.
   *
   * @tparam T is the type of numeric value
   * @param value is an expression
   * @param rhs is an expression
   * @return An expression that represents the subtraction `value`*`rhs`
   */
  template<typename T,
           typename
           = typename std::enable_if_t<std::is_arithmetic_v<T>, T>>
  friend Expression<C> operator*(const T value, const Expression<C> &rhs)
  {
    if (value == 0) {
      return 0;
    }

    if (value == 1) {
      return rhs;
    }

    low_level::base_expression_type<C> *rhs_ex = rhs._ex->clone();

    return Expression<C>(rhs_ex->multiply(static_cast<C>(value)));
  }

  /**
   * @brief Divide a constant value by an expression
   *
   * This method builds an expression that represents the division between
   * a constant value and an expression.
   *
   * @tparam T is the type of numeric value
   * @param value is a constant value
   * @param rhs is an expression
   * @return An expression that represents the division `lhs`/`value`
   */
  template<typename T,
           typename
           = typename std::enable_if_t<std::is_arithmetic_v<T>, T>>
  friend Expression<C> operator/(const T value, const Expression<C> &rhs)
  {
    if (value == 0) {
      return 0;
    }

    return Expression<C>(value) / rhs;
  }

  /**
   * @brief Sum an expression and a numeric value
   *
   * This method builds an expression that represents the sum of the two
   * expressions passed as parameters.
   *
   * @tparam T is the type of numeric value
   * @param value is a value
   * @param rhs is an expression
   * @return An expression that represents the sum `value`+`rhs`
   */
  template<typename T,
           typename
           = typename std::enable_if_t<std::is_arithmetic_v<T>, T>>
  friend inline Expression<C> operator+(const T value,
                                        const Expression<C> &rhs)
  {
    return rhs + value;
  }

  /**
   * @brief Subtract a numeric value from an expression
   *
   * This method builds an expression that represents the subtraction
   * of a numeric value from an expression.
   *
   * @tparam T is the type of numeric value
   * @param lhs is an expression
   * @param value is an expression
   * @return An expression that represents the subtraction `lhs`-`value`
   */
  template<typename T,
           typename
           = typename std::enable_if_t<std::is_arithmetic_v<T>, T>>
  friend inline Expression<C> operator-(const Expression<C> &lhs,
                                        const T value)
  {
    return (lhs + (-value));
  }

  /**
   * @brief Multiply an expression and a numeric value
   *
   * This method builds an expression that represents the product
   * of an expression and a numeric value.
   *
   * @tparam T is the type of numeric value
   * @param lhs is an expression
   * @param value is an expression
   * @return An expression that represents the subtraction `lhs`*`value`
   */
  template<typename T,
           typename
           = typename std::enable_if_t<std::is_arithmetic_v<T>, T>>
  friend inline Expression<C> operator*(const Expression<C> &lhs,
                                        const T value)
  {
    return value * lhs;
  }

  /**
   * @brief Divide an expression by a constant
   *
   * This method builds an expression that represents the division between an
   * expression and a constant value.
   *
   * @tparam C1 is the type of numeric constants in the expression
   * @tparam T is the type of numeric value
   * @param lhs is an expression
   * @param value is a constant value
   * @return An expression that represents the division `lhs`/`value`
   */
  template<typename T,
           typename
           = typename std::enable_if_t<std::is_arithmetic_v<T>, T>>
  friend inline Expression<C> operator/(const Expression<C> &lhs,
                                        const T value)
  {
    if (value == 1) {
      return lhs;
    }

    return lhs * static_cast<C>(1 / static_cast<C>(value));
  }

  /**
   * @brief Get the sum between two expressions
   *
   * This method builds an expression that represents the sum
   * between two expressions.
   *
   * @tparam T is the type of the constant values in the expression
   * @param lhs is an expression
   * @param rhs is another expression
   * @return the expression that represents `lhs`+`rhs`
   */
  template<typename T>
  friend Expression<T> operator+(Expression<T> &&lhs, Expression<T> &&rhs);

  /**
   * @brief Get the subtraction between two expressions
   *
   * This method builds an expression that represents the subtractions
   * between two expressions.
   *
   * @tparam T is the type of the constant values in the expression
   * @param lhs is an expression
   * @param rhs is another expression
   * @return the expression that represents `lhs`-`rhs`
   */
  template<typename T>
  friend Expression<T> operator-(Expression<T> &&lhs, Expression<T> &&rhs);

  /**
   * @brief Get the multiplication between two expressions
   *
   * This method builds an expression that represents the multiplication
   * between two expressions.
   *
   * @tparam T is the type of the constant values in the expression
   * @param lhs is an expression
   * @param rhs is another expression
   * @return the expression that represents `lhs`*`rhs`
   */
  template<typename T>
  friend Expression<T> operator*(Expression<T> &&lhs, Expression<T> &&rhs);

  /**
   * @brief Get the division between two expressions
   *
   * This method builds an expression that represents the division between
   * two expressions.
   *
   * @tparam T is the type of the constant values in the expression
   * @param lhs is the dividend expression
   * @param rhs is the divisor expression
   * @return the expression that represents `lhs`/`rhs`
   */
  template<typename T>
  friend Expression<T> operator/(Expression<T> &&lhs, Expression<T> &&rhs);

  /**
   * @brief Swap two expressions
   *
   * @tparam T is the type of the constant values in the expression
   * @param a is one of the expressions to be swapped
   * @param b is the other one expression to be swapped
   */
  template<typename T>
  friend void std::swap(Expression<T> &a, Expression<T> &b);

  /**
   * @brief Print an expression in a stream
   *
   * @tparam T is the type of the constant values in the expression
   * @param os is the output stream
   * @param ex is the expression to be printed
   * @return the output stream
   */
  template<typename T>
  friend std::ostream &std::operator<<(std::ostream &os,
                                       const Expression<T> &ex);
};

namespace low_level
{
/**
 * @brief A low-level base class to represent algebraic expression
 *
 * @tparam C is the type of numeric constants
 */
template<typename C>
class base_expression_type
{
public:
  typedef typename Symbol<C>::SymbolIdType SymbolIdType;

  /**
   * @brief Build an empty expression
   */
  base_expression_type() {}

  /**
   * @brief Get the expression type
   *
   * @return The type of the expression
   */
  virtual ExpressionType type() const = 0;

  /**
   * @brief Add an expression to the current one
   *
   * This method builds an expression that represents the sum among the
   * current object and the parameter. This is done by merging the two
   * representations and deallocating the non-necessary components. The
   * original expressions depicted by the current object and by the
   * parameter will not be available anymore after the call.
   *
   * @param op is the expression to be added
   * @return a pointer to an expression that represents the sum among
   *         the current object and the parameter
   */
  virtual base_expression_type<C> *add(base_expression_type<C> *op)
  {
    if (!(this->has_symbols()) && this->evaluate() == 0) {
      delete this;

      return op;
    }
    switch (op->type()) {
    case CONSTANT: {
      C op_value = op->evaluate();
      if (op_value == 0) {
        delete op;

        return this;
      }

      finite_sum_type<C> *sum = new finite_sum_type<C>();

      sum->_constant = op_value;

      delete op;

      return sum->add(this);
    }
    case SYMBOL:
    case FINITE_PROD: {
      finite_sum_type<C> *sum = new finite_sum_type<C>();

      sum->_sum.push_back(op);
      return sum->add(this);
    }
    case FINITE_SUM:
    default:
      return op->add(this);
    }
  }

  /**
   * @brief Subtract an expression from the current one
   *
   * This method builds an expression that represents the subtraction
   * of the parameter from the current object. This is done by merging
   * the two representations and deallocating the non-necessary components.
   * The original expressions depicted by the current object and by the
   * parameter will not be available anymore after the execution.
   *
   * @param op is the expression to be subtracted
   * @return a pointer to an expression that represents the subtraction
   *         of the parameter from the current object
   */
  virtual base_expression_type<C> *subtract(base_expression_type<C> *op)
  {
    if (!(op->has_symbols()) && op->evaluate() == 0) {
      delete op;

      return this;
    }
    return this->add(op->complement());
  }

  /**
   * @brief Multiply an expression to the current one
   *
   * This method builds an expression that represents the multiplication
   * of the parameter and the current object. This is done by merging
   * the two representations and deallocating the non-necessary components.
   * The original expressions depicted by the current object and by the
   * parameter will not be available anymore after the call.
   *
   * @param op is the expression that multiplies the current object
   * @return a pointer to an expression that represents the multiplication
   *         of the parameter and the current object
   */
  virtual base_expression_type<C> *multiply(base_expression_type<C> *op)
  {
    if (!(this->has_symbols()) && this->evaluate() == 1) {
      delete this;

      return op;
    }
    if (!(this->has_symbols()) && this->evaluate() == 0) {
      delete op;

      return this;
    }
    switch (op->type()) {
    case CONSTANT: {
      const C op_value = op->evaluate();
      if (op_value == 1) {
        delete op;

        return this;
      }

      finite_prod_type<C> *prod = new finite_prod_type<C>();

      prod->_constant = op_value;

      delete op;

      return prod->multiply(this);
    }
    case SYMBOL:
    case FINITE_SUM: {
      finite_prod_type<C> *prod = new finite_prod_type<C>();

      prod->_numerator.push_back(op);
      return prod->multiply(this);
    }
    case FINITE_PROD:
    default:
      return op->multiply(this);
    }
  }

  /**
   * @brief Divide the current expression by another one
   *
   * This method builds an expression that represents the division
   * of the current object by the parameter. This is done by merging
   * the two representations and deallocating the non-necessary components.
   * The original expressions depicted by the current object and by the
   * parameter will not be available anymore after the call.
   *
   * @param op is the expression that divides the current object
   * @return a pointer to an expression that represents the division
   *         of the current object by the parameter
   */
  virtual base_expression_type<C> *divided_by(base_expression_type<C> *op)
  {
    if (!(this->has_symbols()) && this->evaluate() == 0) {
      delete op;

      return this;
    }
    switch (op->type()) {
    case CONSTANT: {
      const C op_value = op->evaluate();
      if (op_value == 1) {
        delete op;

        return this;
      }
      finite_prod_type<C> *prod = new finite_prod_type<C>();

      prod->_constant = 1 / op_value;

      delete op;

      return prod->multiply(this);
    }
    case SYMBOL:
    case FINITE_SUM: {
      finite_prod_type<C> *prod = new finite_prod_type<C>();

      prod->_denominator.push_back(op);
      return prod->multiply(this);
    }
    case FINITE_PROD:
    default:
      finite_prod_type<C> *prod = new finite_prod_type<C>();
      finite_prod_type<C> *op_prod = (finite_prod_type<C> *)op;

      std::swap(prod->_denominator, op_prod->_numerator);
      std::swap(prod->_numerator, op_prod->_denominator);
      prod->_constant = 1 / op_prod->_constant;
      return prod->multiply(this);
    }
  }

  /**
   * @brief Add a value to the current expression
   *
   * This method builds an expression that represents the sum among the
   * current object and a constant value. The original expression
   * depicted by the current object will not be available anymore after
   * the execution.
   *
   * @param value is the constant value to be added
   * @return a pointer to an expression that represents the sum among
   *         the current object and the constant value
   */
  virtual base_expression_type<C> *add(const C &value)
  {
    finite_sum_type<C> *sum = new finite_sum_type<C>();

    sum->_constant = value;
    sum->_sum.push_back(this);

    return sum;
  }

  /**
   * @brief Subtract a value to the current expression
   *
   * This method builds an expression that represents the subtraction
   * of a constant value from the current expression. The original
   * expression depicted by the current object will not be available
   * anymore after the call.
   *
   * @param value is the constant value to be subtracted
   * @return a pointer to an expression that represents the subtraction
   *         among the current object and the constant value
   */
  virtual inline base_expression_type<C> *subtract(const C &value)
  {
    const C compl_value = -value;

    return this->add(compl_value);
  }

  /**
   * @brief Multiply a value to the current expression
   *
   * This method builds an expression that represents the product
   * between the current object and a constant value. The original
   * expression depicted by the current object will not be available
   * anymore after the call.
   *
   * @param value is the constant value to be multiplied
   * @return a pointer to an expression that represents the product
   *         between the current object and the constant value
   */
  virtual base_expression_type<C> *multiply(const C &value)
  {
    finite_prod_type<C> *prod = new finite_prod_type<C>();

    prod->_constant = value;
    prod->_numerator.push_back(this);

    return prod;
  }

  /**
   * @brief Divide the current expression by a constant value
   *
   * This method builds an expression that represents the division
   * of the current object by the parameter. The original
   * expression depicted by the current object will not be
   * available anymore after the call.
   *
   * @param value is the value that divides the current object
   * @return a pointer to an expression that represents the division
   *         of the current object by the constant value
   */
  virtual base_expression_type<C> *divided_by(const C &value)
  {
    const C rec_value = 1 / value;

    return this->multiply(rec_value);
  }

  /**
   * @brief Get the rational polynomial representation of an expression
   *
   * @return the rational polynomial representation of this expression
   */
  virtual base_expression_type<C> *get_rational_form() const = 0;

  /**
   * @brief Establish whether is a polynomial expression
   *
   * @return `true` if and only if the expression is a polynomial
   */
  virtual bool is_a_polynomial() const = 0;

  /**
   * @brief Establish whether is in rational polynomial form
   *
   * @return `true` if and only if the expression is either
   *        polynomial or the ratio between two polynomials
   */
  virtual bool is_in_rational_form() const = 0;

  /**
   * @brief Complement the expression
   *
   * @return the complement expression
   */
  virtual base_expression_type<C> *complement() = 0;

  /**
   * @brief Multiply all the products in a list by the current object
   *
   * This method multiplies all the product in a list by the current object and
   * returns the list of the products. The list of products will be not
   * available anymore after the call.
   *
   * @param prods is a list of products
   * @return the list of the products between the current object and the
   *         products in the `prods`
   */
  virtual std::list<base_expression_type<C> *>
  multiply(std::list<base_expression_type<C> *> &prods)
  {
    std::list<base_expression_type<C> *> result;

    for (auto p_it = std::begin(prods); p_it != std::end(prods); ++p_it) {
      result.push_back((*p_it)->multiply(this->clone()));
    }
    prods.clear();

    return result;
  }

  /**
   * @brief Turn the expression into a sum of products, a constant, or a
   * symbol.
   *
   * This method returns an expression that is sum of products, a constant,
   * or a symbol and that is algebraic equivalent to the current object.
   *
   * @return A sum of products, a constant, or a symbol
   */
  virtual base_expression_type<C> *expand() const
  {
    return this->clone();
  }

  /**
   * @brief Clone the current object
   *
   * @return A clone of the current object
   */
  virtual base_expression_type<C> *clone() const = 0;

  /**
   * @brief Get the ids of the symbols in the expression
   *
   * @return The set of the symbol ids in the expression
   */
  virtual std::set<SymbolIdType> get_symbol_ids() const = 0;

  /**
   * @brief Get the ids of the symbols in the expression
   *
   * @param max_cardinality is the maximum number of symbol ids to be returned
   * @return A set, whose cardinality is `max_cardinality` at most, of the
   * symbol identifiers in the expression
   */
  virtual std::set<SymbolIdType>
  get_symbol_ids(const unsigned int max_cardinality) const = 0;

  /**
   * @brief Test the presence of any symbol
   *
   * @return `true` if and only if the expression as one symbol at least
   */
  virtual bool has_symbols() const = 0;

  /**
   * @brief Replace symbol occurrences by using expressions
   *
   * This method replaces any occurrence of the symbols whose ids are in
   * the domain of `replacement` with a clone of the corresponding expressions.
   * This method updates this object.
   *
   * @param replacements associates symbol ids to their replacements
   * @return A pointer to the updated object
   */
  virtual base_expression_type<C> *replace(
      const std::map<SymbolIdType, base_expression_type<C> *> &replacements)
      = 0;

  /**
   * @brief Numerically evaluate the expression
   *
   * @return The numeric evaluation of the expression
   */
  virtual C evaluate() const = 0;

  /**
   * @brief Apply a symbolic interpretation to an expression
   *
   * This method evaluates the expression using the symbol interpretation
   * provided as parameter.
   *
   * @param interpretation is an interpretation for the symbols
   * @return the expression obtained from the current object
   *         and `interpretation`
   */
  virtual base_expression_type<C> *
  apply(const std::map<SymbolIdType, C> &interpretation) const
      = 0;

  /**
   * @brief Get the coefficient of a term having a given degree
   *
   * @param symbol_id is the symbol id whose term coefficient is aimed
   * @param degree is the degree of the aimed term
   * @return The coefficient of the term in which the symbol having id
   *         `symbol_id` has degree `degree`
   */
  virtual base_expression_type<C> *get_coeff(const SymbolIdType &symbol_id,
                                             const int &degree) const
      = 0;

  /**
   * @brief Get the coefficients of the polynomial expression
   *
   * @param symbol_id is the symbol id whose term coefficient is aimed
   * @return The coefficients of the polynomial expression on the
   *         symbol whose id is `symbol_id`
   */
  virtual std::map<int, base_expression_type<C> *>
  get_coeffs(const SymbolIdType &symbol_id) const = 0;

  /**
   * @brief Get the degree of a symbol
   *
   * @param symbol_id is the id of the symbol whose degree is aimed
   * @return the degree of the parameter in the expression
   */
  virtual int degree(const SymbolIdType &symbol_id) const = 0;

  /**
   * @brief Print the expression in an output stream
   *
   * @param os is the output stream in which the expression must be printed
   */
  virtual void print(std::ostream &os) const = 0;

  /**
   * @brief Check whether this is the constant 0
   *
   * @return `true` if and only if this is the constant 0
   */
  bool is_zero() const
  {
    return (type() == CONSTANT
            && ((constant_type<C> *)this)->get_value() == 0);
  }

  /**
   * @brief Compute the derivative of an expression
   *
   * This method computes the derivative of the current expression
   * with respect to a symbol provided as the parameter.
   *
   * @param symbol_id is the id of the symbol with respect to the
   *               derivative must be computed
   * @return the derivative of the current expression with respect
   *         to the symbol associated to `symbol_id`
   */
  virtual base_expression_type<C> *
  get_derivative_wrt(const SymbolIdType &symbol_id) const
      = 0;

  /**
   * @brief Destroy the base expression type object
   */
  virtual ~base_expression_type() {}
};

/**
 * @brief A class to represent constants
 *
 * @tparam C is the type of numeric constants
 */
template<typename C>
class constant_type : public base_expression_type<C>
{
  C _value; //!< The numerical value of the constant.
public:
  typedef typename Symbol<C>::SymbolIdType SymbolIdType;

  /**
   * @brief Build a new constant expression
   *
   * @param value is the value of the constant
   */
  constant_type(const C value): base_expression_type<C>(), _value(value) {}

  /**
   * @brief Copy constructor
   *
   * @param orig is the model for the new constant
   */
  constant_type(const constant_type<C> &orig):
      base_expression_type<C>(), _value(orig._value)
  {
  }

  /**
   * @brief Get the constant value
   *
   * @return the value of the constant
   */
  const C &get_value() const
  {
    return _value;
  }

  /**
   * @brief Get the expression type
   *
   * @return The value CONSTANT
   */
  ExpressionType type() const
  {
    return CONSTANT;
  }

  /**
   * @brief Add an expression to the current one
   *
   * This method builds an expression that represents the sum among the
   * current object and the parameter. This is done by merging the two
   * representations and deallocating the non-necessary components. The
   * original expressions depicted by the current object and by the
   * parameter will not be available anymore after the call.
   *
   * @param op is the expression to be added
   * @return a pointer to an expression that represents the sum among
   *         the current object and the parameter
   */
  base_expression_type<C> *add(base_expression_type<C> *op)
  {
    if (this->evaluate() == 0) {
      delete this;

      return op;
    }
    switch (op->type()) {
    case CONSTANT:
      _value += op->evaluate();

      delete op;

      return this;

    default:
      return op->add(this);
    }
  }

  /**
   * @brief Subtract an expression from the current one
   *
   * This method builds an expression that represents the subtraction
   * of the parameter from the current object. This is done by merging
   * the two representations and deallocating the non-necessary components.
   * The original expressions depicted by the current object and by the
   * parameter will not be available anymore after the call.
   *
   * @param op is the expression to be subtracted
   * @return a pointer to an expression that represents the subtraction
   *         of the parameter from the current object
   */
  base_expression_type<C> *subtract(base_expression_type<C> *op)
  {
    if (this->evaluate() == 0) {
      delete this;

      constant_type<C> *minus1 = new constant_type<C>(-1);

      return minus1->multiply(op);
    }
    switch (op->type()) {
    case CONSTANT:
      _value -= op->evaluate();

      delete op;

      return this;

    default:
      constant_type<C> *minus1 = new constant_type<C>(-1);
      return this->add(minus1->multiply(op));
    }
  }

  /**
   * @brief Multiply an expression to the current one
   *
   * This method builds an expression that represents the multiplication
   * of the parameter and the current object. This is done by merging
   * the two representations and deallocating the non-necessary components.
   * The original expressions depicted by the current object and by the
   * parameter will not be available anymore after the call.
   *
   * @param op is the expression that multiplies the current object
   * @return a pointer to an expression that represents the multiplication
   *         of the parameter and the current object
   */
  base_expression_type<C> *multiply(base_expression_type<C> *op)
  {
    if (this->evaluate() == 0) {
      delete op;

      return this;
    }
    if (this->evaluate() == 1) {
      delete this;

      return op;
    }
    switch (op->type()) {
    case CONSTANT:
      _value *= op->evaluate();

      delete op;

      return this;

    default:
      return op->multiply(this);
    }
  }

  /**
   * @brief Divide the current expression by another one
   *
   * This method builds an expression that represents the division
   * of the current object by the parameter. This is done by merging
   * the two representations and deallocating the non-necessary components.
   * The original expressions depicted by the current object and by the
   * parameter will not be available anymore after the call.
   *
   * @param op is the expression that divides the current object
   * @return a pointer to an expression that represents the division
   *         of the current object by the parameter
   */
  base_expression_type<C> *divided_by(base_expression_type<C> *op)
  {
    if (this->evaluate() == 0) {
      delete op;

      return this;
    }
    switch (op->type()) {
    case CONSTANT:
      _value /= op->evaluate();

      delete op;

      return this;

    default: {
      base_expression_type<C> *prod = new finite_prod_type<C>();

      prod = prod->divided_by(op);

      return prod->multiply(this);
    }
    }
  }

  /**
   * @brief Add a value to the current expression
   *
   * This method builds an expression that represents the sum among the
   * current object and a constant value. The original expression
   * depicted by the current object will not be available anymore after
   * the execution.
   *
   * @param value is the constant value to be added
   * @return a pointer to an expression that represents the sum among
   *         the current object and the constant value
   */
  base_expression_type<C> *add(const C &value)
  {
    _value += value;

    return this;
  }

  /**
   * @brief Subtract a value to the current expression
   *
   * This method builds an expression that represents the subtraction
   * of a constant value from the current expression. The original
   * expression depicted by the current object will not be available
   * anymore after the call.
   *
   * @param value is the constant value to be subtracted
   * @return a pointer to an expression that represents the subtraction
   *         of the constant value from the current object
   */
  base_expression_type<C> *subtract(const C &value)
  {
    _value -= value;

    return this;
  }

  /**
   * @brief Multiply a value to the current expression
   *
   * This method builds an expression that represents the product
   * between the current object and a constant value. The original
   * expression depicted by the current object will not be available
   * anymore after the call.
   *
   * @param value is the constant value to be multiplied
   * @return a pointer to an expression that represents the product
   *         between the current object and the constant value
   */
  base_expression_type<C> *multiply(const C &value)
  {
    _value *= value;

    return this;
  }

  /**
   * @brief Divide the current expression by a constant value
   *
   * This method builds an expression that represents the division
   * of the current object by the parameter. The original
   * expression depicted by the current object will not be
   * available anymore after the call.
   *
   * @param value is the value that divides the current object
   * @return a pointer to an expression that represents the division
   *         of the current object by the constant value
   */
  base_expression_type<C> *divided_by(const C &value)
  {
    _value /= value;

    return this;
  }

  /**
   * @brief Get the rational polynomial representation of an expression
   *
   * @return the rational polynomial representation of this expression
   */
  base_expression_type<C> *get_rational_form() const
  {
    return this->clone();
  }

  /**
   * @brief Establish whether is a polynomial expression
   *
   * @return `true` if and only if the expression is a polynomial
   */
  bool is_a_polynomial() const
  {
    return true;
  }

  /**
   * @brief Establish whether is in rational polynomial form
   *
   * @return `true` if and only if the expression is either
   *        polynomial or the ratio between two polynomials
   */
  bool is_in_rational_form() const
  {
    return true;
  }

  /**
   * @brief Complement the expression
   *
   * @return the complement expression
   */
  base_expression_type<C> *complement()
  {
    _value = -_value;

    return this;
  }

  /**
   * @brief Replace symbol occurrences by using expressions
   *
   * This method replaces any occurrence of the symbols whose ids are in
   * the domain of `replacement` with a clone of the corresponding expressions.
   * This method updates this object.
   *
   * @param replacements associates symbol ids to their replacements
   * @return A pointer to the updated object
   */
  base_expression_type<C> *replace(
      const std::map<SymbolIdType, base_expression_type<C> *> &replacements)
  {
    (void)replacements;

    return this;
  }

  /**
   * @brief Clone the current object
   *
   * @return A clone of the current object
   */
  base_expression_type<C> *clone() const
  {
    return new constant_type<C>(*this);
  }

  /**
   * @brief Get the ids of the symbols in the expression
   *
   * @return The set of the symbol ids in the expression
   */
  std::set<SymbolIdType> get_symbol_ids() const
  {
    return std::set<SymbolIdType>();
  }

  /**
   * @brief Get the ids of the symbols in the expression
   *
   * @param max_cardinality is the maximum number of symbol ids to be returned
   * @return A set, whose cardinality is `max_cardinality` at most, of the
   * symbol identifiers in the expression
   */
  std::set<SymbolIdType>
  get_symbol_ids(const unsigned int max_cardinality) const
  {
    (void)max_cardinality;

    return std::set<SymbolIdType>();
  }

  /**
   * @brief Test the presence of any symbol
   *
   * @return `true` if and only if the expression as one symbol at least
   */
  bool has_symbols() const
  {
    return false;
  }

  /**
   * @brief Numerically evaluate the expression
   *
   * @return The numeric evaluation of the expression
   */
  C evaluate() const
  {
    return _value;
  }

  /**
   * @brief Apply a symbolic interpretation to an expression
   *
   * This method evaluates the expression using the symbol interpretation
   * provided as parameter.
   *
   * @param interpretation is an interpretation for the symbols
   * @return the expression obtained from the current object
   *         and `interpretation`
   */
  base_expression_type<C> *
  apply(const std::map<SymbolIdType, C> &interpretation) const
  {
    (void)interpretation;

    return this->clone();
  }

  /**
   * @brief Get the coefficient of a term having a given degree
   *
   * @param symbol_id is the symbol id whose term coefficient is aimed
   * @param degree is the degree of the aimed term
   * @return The coefficient of the term in which the symbol having id
   *         `symbol_id` has degree `degree`
   */
  base_expression_type<C> *get_coeff(const SymbolIdType &symbol_id,
                                     const int &degree) const
  {
    (void)symbol_id;

    return new constant_type(degree == 0 ? _value : 0);
  }

  /**
   * @brief Get the coefficients of the polynomial expression
   *
   * @param symbol_id is the symbol id whose term coefficient is aimed
   * @return The coefficients of the polynomial expression on the
   *         symbol whose id is `symbol_id`
   */
  std::map<int, base_expression_type<C> *>
  get_coeffs(const SymbolIdType &symbol_id) const
  {
    (void)symbol_id;

    std::map<int, base_expression_type<C> *> res;

    res[0] = new constant_type<C>(_value);

    return res;
  }

  /**
   * @brief Get the degree of a symbol
   *
   * @param symbol_id is the id of the symbol whose degree is aimed
   * @return the degree of the parameter in the expression, i.e., 0
   */
  int degree(const SymbolIdType &symbol_id) const
  {
    (void)symbol_id;

    return 0;
  }

  /**
   * @brief Compute the derivative of an expression
   *
   * This method computes the derivative of the current expression
   * with respect to a symbol provided as the parameter.
   *
   * @param symbol_id is the id of the symbol with respect to the
   *               derivative must be computed
   * @return the derivative of the current expression with respect
   *         to the symbol associated to `symbol_id`
   */
  base_expression_type<C> *
  get_derivative_wrt(const SymbolIdType &symbol_id) const
  {
    (void)symbol_id;

    return new constant_type<C>(0);
  }

  /**
   * @brief Print the expression in an output stream
   *
   * @param os is the output stream in which the expression must be printed
   */
  void print(std::ostream &os) const
  {
    os << (_value == 0 ? 0 : _value);
  }

  ~constant_type() {}
};

/**
 * @brief A class to represent finite sums
 *
 * @tparam C is the type of numeric constants
 */
template<typename C>
class finite_sum_type : public base_expression_type<C>
{
  C _constant;
  std::list<base_expression_type<C> *>
      _sum; //!< The finite list of non-constant expressions.
public:
  typedef typename Symbol<C>::SymbolIdType SymbolIdType;

  /**
   * @brief Build an empty finite sum
   */
  finite_sum_type(): base_expression_type<C>(), _constant(0), _sum() {}

  /**
   * @brief Copy constructor
   *
   * @param orig is the model for the new finite sum
   */
  finite_sum_type(const finite_sum_type<C> &orig):
      base_expression_type<C>(), _constant(orig._constant), _sum()
  {
    for (auto it = std::begin(orig._sum); it != std::end(orig._sum); ++it) {
      _sum.push_back((*it)->clone());
    }
  }

  /**
   * @brief Build a finite sum of products
   *
   * @param orig is a list of products that must be added
   */
  finite_sum_type(std::list<base_expression_type<C> *> &model):
      base_expression_type<C>(), _constant(0), _sum()
  {
    for (auto it = std::begin(model); it != std::end(model); ++it) {
      _sum.push_back(*it);
    }

    model.clear();
  }

  /**
   * @brief Get the expression type
   *
   * @return The value FINITE_SUM
   */
  ExpressionType type() const
  {
    return FINITE_SUM;
  }

  /**
   * @brief Add an expression to the current one
   *
   * This method builds an expression that represents the sum among the
   * current object and the parameter. This is done by merging the two
   * representations and deallocating the non-necessary components. The
   * original expressions depicted by the current object and by the
   * parameter will not be available anymore after the call.
   *
   * @param op is the expression to be added
   * @return a pointer to an expression that represents the sum among
   *         the current object and the parameter
   */
  base_expression_type<C> *add(base_expression_type<C> *op)
  {
    switch (op->type()) {
    case FINITE_SUM: {
      finite_sum_type<C> *op_sum = (finite_sum_type<C> *)op;

      _constant += op_sum->_constant;
      _sum.splice(std::end(_sum), op_sum->_sum);

      delete op;

      return this;
    }
    case CONSTANT:
      _constant += op->evaluate();

      delete op;

      if (_sum.size() == 0) {
        constant_type<C> *res = new constant_type<C>(_constant);

        delete this;

        return res;
      }
      return this;

    default:
      _sum.push_back(op);
    }

    return this;
  };

  /**
   * @brief Add a value to the current expression
   *
   * This method builds an expression that represents the sum among the
   * current object and a constant value. The original expression
   * depicted by the current object will not be available anymore after
   * the execution.
   *
   * @param value is the constant value to be added
   * @return a pointer to an expression that represents the sum among
   *         the current object and the constant value
   */
  base_expression_type<C> *add(const C &value)
  {
    _constant += value;

    return this;
  }

  /**
   * @brief Subtract a value to the current expression
   *
   * This method builds an expression that represents the subtraction
   * of a constant value from the current expression. The original
   * expression depicted by the current object will not be available
   * anymore after the call.
   *
   * @param value is the constant value to be subtracted
   * @return a pointer to an expression that represents the sum among
   *         the current object and the constant value
   */
  base_expression_type<C> *subtract(const C &value)
  {
    _constant -= value;

    return this;
  }

  /**
   * @brief Multiply all the products in a list by the current object
   *
   * This method multiplies all the product in a list by the current object and
   * returns the list of the products. The list of products will be not
   * available anymore after the call.
   *
   * @param prods is a list of products
   * @return the list of the products between the current object and the
   *         products in the `prods`
   */
  std::list<base_expression_type<C> *>
  multiply(std::list<base_expression_type<C> *> &prods)
  {
    std::list<base_expression_type<C> *> result;

    for (auto s_it = std::begin(_sum); s_it != std::end(_sum); ++s_it) {
      for (auto p_it = std::begin(prods); p_it != std::end(prods); ++p_it) {
        base_expression_type<C> *p_clone = (*p_it)->clone();

        result.push_back(p_clone->multiply((*s_it)->clone()));
      }
      delete *s_it;
    }
    _sum.clear();

    if (_constant != 0) {
      for (auto p_it = std::begin(prods); p_it != std::end(prods); ++p_it) {
        result.push_back((*p_it)->multiply(_constant));
      }
    }
    prods.clear();

    return result;
  }

  /**
   * @brief Get the rational polynomial representation of an expression
   *
   * @return the rational polynomial representation of this expression
   */
  base_expression_type<C> *get_rational_form() const
  {
    std::vector<base_expression_type<C> *> numerators;
    std::vector<base_expression_type<C> *> denominators;

    numerators.reserve(_sum.size() + 1);
    denominators.reserve(_sum.size() + 1);

    numerators.push_back(new constant_type<C>(_constant));
    denominators.push_back(new constant_type<C>(1));

    unsigned int real_denoms = 0;

    for (const auto &term: _sum) {

      base_expression_type<C> *rational_term = term->get_rational_form();

      if (rational_term->type() == FINITE_PROD) {
        auto r_term_prod = static_cast<finite_prod_type<C> *>(rational_term);

        if (r_term_prod->_denominator.size() > 0) {
          auto den = new finite_prod_type<C>(1);

          std::swap(den->_numerator, r_term_prod->_denominator);

          denominators.push_back(den);
          real_denoms += 1;
        } else {
          denominators.push_back(new constant_type<C>(1));
        }
        if (r_term_prod->_numerator.size() == 1
            && r_term_prod->_constant == 1) {
          numerators.push_back(r_term_prod->_numerator.front());

          r_term_prod->_numerator.clear();

          delete rational_term;
        } else {
          numerators.push_back(r_term_prod);
        }
      } else {
        numerators.push_back(rational_term);
        denominators.push_back(new constant_type<C>(1));
      }
    }

    if (real_denoms > 0) {
      auto prod = new finite_prod_type<C>();

      for (size_t j = 0; j < denominators.size(); ++j) {
        for (size_t i = 0; i < numerators.size(); ++i) {
          if (i != j) {
            numerators[i] = numerators[i]->multiply(denominators[j]->clone());
          }
        }

        prod = static_cast<finite_prod_type<C> *>(
            prod->divided_by(denominators[j]));
      }

      base_expression_type<C> *num_sum = new finite_sum_type<C>();
      for (base_expression_type<C> *num: numerators) {
        num_sum = num_sum->add(num);
      }

      prod->_numerator.push_back(num_sum);

      return prod;
    }

    return this->clone();
  }

  /**
   * @brief Establish whether is a polynomial expression
   *
   * @return `true` if and only if the expression is a polynomial
   */
  bool is_a_polynomial() const
  {
    for (const auto term: _sum) {
      if (!term->is_a_polynomial()) {
        return false;
      }
    }

    return true;
  }

  /**
   * @brief Establish whether is in rational polynomial form
   *
   * @return `true` if and only if the expression is either
   *        polynomial or the ratio between two polynomials
   */
  bool is_in_rational_form() const
  {
    return this->is_a_polynomial();
  }

  /**
   * @brief Complement the expression
   *
   * @return the complement expression
   */
  base_expression_type<C> *complement()
  {
    _constant = -_constant;
    for (auto it = std::begin(_sum); it != std::end(_sum); ++it) {
      (*it) = (*it)->complement();
    }

    return this;
  }

  /**
   * @brief Turn the expression into a sum of products, a constant, or a
   * symbol.
   *
   * This method returns an expression that is sum of products, a constant,
   * or a symbol and that is algebraic equivalent to the current object.
   *
   * @return A sum of products, a constant, or a symbol
   */
  base_expression_type<C> *expand() const
  {
    finite_sum_type<C> *new_obj = new finite_sum_type<C>();
    std::list<base_expression_type<C> *> &new_sum = new_obj->_sum;
    C &constant = new_obj->_constant;

    constant = _constant;

    for (auto it = std::begin(_sum); it != std::end(_sum); ++it) {
      base_expression_type<C> *ex_it = (*it)->expand();

      switch (ex_it->type()) {
      case FINITE_SUM: {
        finite_sum_type<C> *ex_it_sum = (finite_sum_type<C> *)ex_it;
        for (auto e_it = std::begin(ex_it_sum->_sum);
             e_it != std::end(ex_it_sum->_sum); ++e_it) {
          if ((*e_it)->type() == CONSTANT) {
            constant += ((constant_type<C> *)(*e_it))->get_value();

            delete *e_it;
          } else {
            new_sum.push_back(*e_it);
          }
        }
        ex_it_sum->_sum.clear();

        delete ex_it_sum;
        break;
      }
      case CONSTANT: {
        constant += ((constant_type<C> *)ex_it)->get_value();

        delete ex_it;
        break;
      }
      default:
        new_sum.push_back(ex_it);
      }
    }

    return new_obj;
  }

  /**
   * @brief Clone the current object
   *
   * @return A clone of the current object
   */
  base_expression_type<C> *clone() const
  {
    return new finite_sum_type<C>(*this);
  }

  /**
   * @brief Get the ids of the symbols in the expression
   *
   * @return The set of the symbol ids in the expression
   */
  std::set<SymbolIdType> get_symbol_ids() const
  {
    std::set<SymbolIdType> ids;

    for (auto it = std::begin(_sum); it != std::end(_sum); ++it) {
      std::set<SymbolIdType> it_ids = (*it)->get_symbol_ids();

      ids.insert(std::begin(it_ids), std::end(it_ids));
    }

    return ids;
  }

  /**
   * @brief Get the ids of the symbols in the expression
   *
   * @param max_cardinality is the maximum number of symbol ids to be returned
   * @return A set, whose cardinality is `max_cardinality` at most, of the
   * symbol identifiers in the expression
   */
  std::set<SymbolIdType>
  get_symbol_ids(const unsigned int max_cardinality) const
  {
    std::set<SymbolIdType> ids;

    unsigned int missing = max_cardinality;
    for (auto it = std::begin(_sum); it != std::end(_sum); ++it) {
      std::set<SymbolIdType> it_ids = (*it)->get_symbol_ids(missing);

      ids.insert(std::begin(it_ids), std::end(it_ids));

      missing = max_cardinality - ids.size();
      if (missing == 0) {
        return ids;
      }
    }

    return ids;
  }

  /**
   * @brief Test the presence of any symbol
   *
   * @return `true` if and only if the expression as one symbol at least
   */
  bool has_symbols() const
  {
    for (auto it = std::begin(_sum); it != std::end(_sum); ++it) {
      if ((*it)->has_symbols()) {
        return true;
      }
    }

    return false;
  }

  /**
   * @brief Replace symbol occurrences by using expressions
   *
   * This method replaces any occurrence of the symbols whose ids are in
   * the domain of `replacement` with a clone of the corresponding expressions.
   * This method updates this object.
   *
   * @param replacements associates symbol ids to their replacements
   * @return A pointer to the updated object
   */
  base_expression_type<C> *replace(
      const std::map<SymbolIdType, base_expression_type<C> *> &replacements)
  {
    std::list<base_expression_type<C> *> new_sum;
    for (auto it = std::begin(_sum); it != std::end(_sum); ++it) {
      base_expression_type<C> *r_it = (*it)->replace(replacements);
      switch (r_it->type()) {
      case FINITE_SUM: {
        finite_sum_type<C> *r_it_sum = (finite_sum_type<C> *)r_it;
        _constant += r_it_sum->_constant;

        new_sum.splice(std::end(new_sum), r_it_sum->_sum);

        delete r_it;
        break;
      }
      case CONSTANT:
        _constant += ((constant_type<C> *)r_it)->get_value();

        delete r_it;
        break;
      default:
        new_sum.push_back(r_it);
      }
    }

    _sum.clear();

    if (new_sum.size() > 1) {
      swap(_sum, new_sum);

      return this;
    }

    if (_constant == 0) {
      if (new_sum.size() == 1) {
        delete this;

        return new_sum.front();
      }
    }

    if (new_sum.size() == 1) {
      if (_constant != 0) {
        swap(_sum, new_sum);

        return this;
      }

      // if the sum consists in at most one element
      // we do not need a sum
      delete this;

      return new_sum.front();
    }

    // if the sum consists in at most one element
    // we do not need a sum

    constant_type<C> *result = new constant_type<C>(_constant);
    delete this;

    return result;
  }

  /**
   * @brief Numerically evaluate the expression
   *
   * @return The numeric evaluation of the expression
   */
  C evaluate() const
  {
    C total = _constant;
    for (auto it = std::begin(_sum); it != std::end(_sum); ++it) {
      total += (*it)->evaluate();
    }

    return total;
  }

  /**
   * @brief Apply a symbolic interpretation to an expression
   *
   * This method evaluates the expression using the symbol interpretation
   * provided as parameter.
   *
   * @param interpretation is an interpretation for the symbols
   * @return the expression obtained from the current object
   *         and `interpretation`
   */
  base_expression_type<C> *
  apply(const std::map<SymbolIdType, C> &interpretation) const
  {
    base_expression_type<C> *total = new constant_type<C>(_constant);
    for (auto it = std::begin(_sum); it != std::end(_sum); ++it) {
      total = total->add((*it)->apply(interpretation));
    }

    return total;
  }

  /**
   * @brief Get the coefficient of a term having a given degree
   *
   * @param symbol_id is the symbol id whose term coefficient is aimed
   * @param degree is the degree of the aimed term
   * @return The coefficient of the term in which the symbol having id
   *         `symbol_id` has degree `degree`
   */
  base_expression_type<C> *get_coeff(const SymbolIdType &symbol_id,
                                     const int &degree) const
  {
    base_expression_type<C> *total_coeff = new finite_sum_type<C>();

    if (degree == 0) {
      ((finite_sum_type<C> *)total_coeff)->_constant = this->_constant;
    }

    for (auto it = std::begin(_sum); it != std::end(_sum); ++it) {
      auto coeff_it = (*it)->get_coeff(symbol_id, degree);

      if (coeff_it != nullptr) {
        total_coeff = total_coeff->add(coeff_it);
      }
    }

    return total_coeff;
  }

  /**
   * @brief Get the coefficients of the polynomial expression
   *
   * @param symbol_id is the symbol id whose term coefficient is aimed
   * @return The coefficients of the polynomial expression on the
   *         symbol whose id is `symbol_id`
   */
  std::map<int, base_expression_type<C> *>
  get_coeffs(const SymbolIdType &symbol_id) const
  {
    std::map<int, base_expression_type<C> *> res;

    res[0] = new constant_type<C>(this->_constant);

    for (auto it = std::begin(_sum); it != std::end(_sum); ++it) {
      auto map_it = (*it)->get_coeffs(symbol_id);
      for (auto it_coeff = std::begin(map_it); it_coeff != std::end(map_it);
           ++it_coeff) {
        if (res.find(it_coeff->first) == std::end(res)) {
          res[it_coeff->first] = it_coeff->second;
        } else {
          res[it_coeff->first] = res[it_coeff->first]->add(it_coeff->second);
        }
      }
    }

    return res;
  }

  /**
   * @brief Get the degree of a symbol
   *
   * @param symbol_id is the id of the symbol whose degree is aimed
   * @return the degree of the parameter in the expression, i.e,
   *         the maximum degree among the addends
   */
  int degree(const SymbolIdType &symbol_id) const
  {
    if (_sum.size() == 0) {
      return 0;
    }

    auto it = std::begin(_sum);
    int max_degree = (*it)->degree(symbol_id);

    for (++it; it != std::end(_sum); ++it) {
      int new_degree = (*it)->degree(symbol_id);

      if (new_degree > max_degree) {
        max_degree = new_degree;
      }
    }

    return max_degree;
  }

  /**
   * @brief Compute the derivative of an expression
   *
   * This method computes the derivative of the current expression
   * with respect to a symbol provided as the parameter.
   *
   * @param symbol_id is the id of the symbol with respect to the
   *               derivative must be computed
   * @return the derivative of the current expression with respect
   *         to the symbol associated to `symbol_id`
   */
  base_expression_type<C> *
  get_derivative_wrt(const SymbolIdType &symbol_id) const
  {
    base_expression_type<C> *res = new constant_type<C>(0);

    for (const auto &term: _sum) {
      res = res->add(term->get_derivative_wrt(symbol_id));
    }

    return res;
  }

  /**
   * @brief Print the expression in an output stream
   *
   * @param os is the output stream in which the expression must be printed
   */
  void print(std::ostream &os) const
  {
    if (_constant != 0 || _sum.size() == 0) {
      os << (_constant == 0 ? 0 : _constant);
    }

    for (auto it = std::begin(_sum); it != std::end(_sum); ++it) {
      if (it != std::begin(_sum) || _constant != 0) {
        os << " + ";
      }
      switch ((*it)->type()) {
      case FINITE_SUM:
        os << "(";
        (*it)->print(os);
        os << ")";
        break;
      default:
        (*it)->print(os);
      }
    }
  }

  ~finite_sum_type()
  {
    for (auto it = std::begin(_sum); it != std::end(_sum); ++it) {
      delete *it;
    }

    _sum.clear();
  }

  template<typename T>
  friend class base_expression_type;

  template<typename T>
  friend struct less;
};

/**
 * @brief A class to represent finite products
 *
 * @tparam C is the type of numeric constants
 */
template<typename C>
class finite_prod_type : public base_expression_type<C>
{
  C _constant; //!< The numeric constant that multiplies this product.
  std::list<base_expression_type<C> *>
      _numerator; //!< The numerator list of this product.
  std::list<base_expression_type<C> *>
      _denominator; //!< The denominator list of this product.

public:
  typedef typename Symbol<C>::SymbolIdType SymbolIdType;

private:
  static void print_list(std::ostream &os,
                         const std::list<base_expression_type<C> *> &_list,
                         const C constant = 1)
  {
    if (constant != 1 || _list.size() == 0) {
      os << (constant == 0 ? 0 : constant);
    }

    for (auto it = std::begin(_list); it != std::end(_list); ++it) {
      if (it != std::begin(_list) || constant != 1) {
        os << "*";
      }
      switch ((*it)->type()) {
      case FINITE_SUM:
        os << "(";
        (*it)->print(os);
        os << ")";
        break;
      default:
        (*it)->print(os);
      }
    }
  }

  /**
   * @brief Constructor
   *
   * @param constant is a constant to initialize the new product
   */
  finite_prod_type(const C &constant):
      base_expression_type<C>(), _constant(constant), _numerator(),
      _denominator()
  {
  }

  /**
   * @brief Compute the derivative of a partial product list
   *
   * @param list of terms whose derivative must be computed
   * @param symbol_id is the id of the symbol with respect to the
   *               derivative must be computed
   * @param pos initial position in the list to be considered
   * @return the derivative with respect to the symbol associated to
   *    `symbol_id` of the product of the expressions in the list
   */
  static base_expression_type<C> *_prod_list_derivative(
      const std::list<base_expression_type<C> *> &list,
      const SymbolIdType &symbol_id,
      typename std::list<base_expression_type<C> *>::const_iterator pos)
  {
    base_expression_type<C> *term1 = (*pos)->get_derivative_wrt(symbol_id);

    if (std::next(pos) == std::end(list)) {
      return term1;
    }

    for (auto pos2 = std::next(pos); pos2 != std::end(list); ++pos2) {
      term1 = term1->multiply((*pos2)->clone());
    }

    base_expression_type<C> *term2 = (*pos)->clone();
    term2 = term2->multiply(_prod_list_derivative(list, symbol_id, ++pos));

    return term1->add(term2);
  }

  /**
   * @brief Compute the derivative of a product list
   *
   * @param list of terms whose derivative must be computed
   * @param symbol_id is the id of the symbol with respect to the
   *               derivative must be computed
   * @param pos initial position in the list to be considered
   * @return the derivative with respect to the symbol associated to
   *    `symbol_id` of the product of the expressions in the list
   */
  static base_expression_type<C> *
  _prod_list_derivative(const std::list<base_expression_type<C> *> &list,
                        const SymbolIdType &symbol_id)
  {
    if (list.size() == 0) {
      return new constant_type<C>(0);
    }

    if (list.size() == 1) {
      return list.front()->get_derivative_wrt(symbol_id);
    }

    return _prod_list_derivative(list, symbol_id, std::begin(list));
  }

public:
  /**
   * @brief Build an empty finite product
   */
  finite_prod_type():
      base_expression_type<C>(), _constant(1), _numerator(), _denominator()
  {
  }

  /**
   * @brief Copy constructor
   *
   * @param orig is a model for the new finite product
   */
  finite_prod_type(const finite_prod_type<C> &orig):
      base_expression_type<C>(), _constant(orig._constant), _numerator(),
      _denominator()
  {
    for (auto it = std::begin(orig._numerator);
         it != std::end(orig._numerator); ++it) {
      _numerator.push_back((*it)->clone());
    }

    for (auto it = std::begin(orig._denominator);
         it != std::end(orig._denominator); ++it) {
      _denominator.push_back((*it)->clone());
    }
  }

  /**
   * @brief Get the expression type
   *
   * @return The value FINITE_PROD
   */
  ExpressionType type() const
  {
    return FINITE_PROD;
  }

  /**
   * @brief Multiply an expression to the current one
   *
   * This method builds an expression that represents the multiplication
   * of the parameter and the current object. This is done by merging
   * the two representations and deallocating the non-necessary components.
   * The original expressions depicted by the current object and by the
   * parameter will not be available anymore after the call.
   *
   * @param op is the expression that multiplies the current object
   * @return a pointer to an expression that represents the multiplication
   *         of the parameter and the current object
   */
  base_expression_type<C> *multiply(base_expression_type<C> *op)
  {
    switch (op->type()) {
    case FINITE_PROD: {
      finite_prod_type<C> *op_prod = (finite_prod_type<C> *)op;

      _constant *= op_prod->_constant;

      _numerator.splice(std::end(_numerator), op_prod->_numerator);
      _denominator.splice(std::end(_denominator), op_prod->_denominator);

      delete op_prod;

      break;
    }
    case CONSTANT: {
      constant_type<C> *op_const = (constant_type<C> *)op;

      _constant *= op_const->get_value();

      delete op_const;

      break;
    }
    default:
      _numerator.push_back(op);
    }

    if (_denominator.size() == 0) {
      if (_constant == 0 || _numerator.size() == 0) {
        constant_type<C> *result = new constant_type<C>(_constant);

        delete this;

        return result;
      }

      if (_constant == 1 && _numerator.size() == 1) {

        base_expression_type<C> *result = _numerator.front();

        _numerator.clear();

        delete this;

        return result;
      }
    }

    return this;
  }

  /**
   * @brief Divide the current expression by another one
   *
   * This method builds an expression that represents the division
   * of the current object by the parameter. This is done by merging
   * the two representations and deallocating the non-necessary components.
   * The original expressions depicted by the current object and by the
   * parameter will not be available anymore after the call.
   *
   * @param op is the expression that divides the current object
   * @return a pointer to an expression that represents the division
   *         of the current object by the parameter
   */
  base_expression_type<C> *divided_by(base_expression_type<C> *op)
  {
    switch (op->type()) {
    case FINITE_PROD: {
      finite_prod_type<C> *op_prod = (finite_prod_type<C> *)op;

      _constant /= op_prod->_constant;

      _numerator.splice(std::end(_numerator), op_prod->_denominator);
      _denominator.splice(std::end(_denominator), op_prod->_numerator);

      delete op;

      break;
    }
    case CONSTANT: {
      constant_type<C> *op_const = (constant_type<C> *)op;

      _constant /= op_const->get_value();

      delete op_const;

      break;
    }
    default:
      _denominator.push_back(op);
    }

    return this;
  }

  /**
   * @brief Multiply a value to the current expression
   *
   * This method builds an expression that represents the product
   * between the current object and a constant value. The original
   * expression depicted by the current object will not be available
   * anymore after the call.
   *
   * @param value is the constant value to be multiplied
   * @return a pointer to an expression that represents the product
   *         between the current object and the constant value
   */
  base_expression_type<C> *multiply(const C &value)
  {
    _constant *= value;

    return this;
  }

  /**
   * @brief Divide the current expression by a constant value
   *
   * This method builds an expression that represents the division
   * of the current object by the parameter. The original
   * expression depicted by the current object will not be
   * available anymore after the call.
   *
   * @param value is the value that divides the current object
   * @return a pointer to an expression that represents the division
   *         of the current object by the constant value
   */
  base_expression_type<C> *divided_by(const C &value)
  {
    _constant /= value;

    return this;
  }

  /**
   * @brief Get the rational polynomial representation of an expression
   *
   * @return the rational polynomial representation of this expression
   */
  base_expression_type<C> *get_rational_form() const
  {
    base_expression_type<C> *prod = new finite_prod_type<C>(_constant);

    for (const auto &term: _numerator) {
      prod = prod->multiply(term->get_rational_form());
    }

    for (const auto &term: _denominator) {
      prod = prod->divided_by(term->get_rational_form());
    }

    return prod;
  }

  /**
   * @brief Establish whether is a polynomial expression
   *
   * @return `true` if and only if the expression is a polynomial
   */
  bool is_a_polynomial() const
  {
    for (const auto term: _denominator) {
      if (term->type() != CONSTANT) {
        return false;
      }
    }

    return true;
  }

  /**
   * @brief Establish whether is in rational polynomial form
   *
   * @return `true` if and only if the expression is either
   *        polynomial or the ratio between two polynomials
   */
  bool is_in_rational_form() const
  {
    for (const auto term: _numerator) {
      if (!term->is_a_polynomial()) {
        return false;
      }
    }

    for (const auto term: _denominator) {
      if (!term->is_a_polynomial()) {
        return false;
      }
    }

    return true;
  }

  /**
   * @brief Complement the expression
   *
   * @return the complement expression
   */
  base_expression_type<C> *complement()
  {
    _constant = -_constant;

    return this;
  }

  /**
   * @brief Turn the expression into a sum of products, a constant, or a
   * symbol.
   *
   * This method returns an expression that is sum of products, a constant,
   * or a symbol and that is algebraic equivalent to the current object.
   *
   * @return A sum of products, a constant, or a symbol
   */
  base_expression_type<C> *expand() const
  {
    base_expression_type<C> *base = new finite_prod_type<C>(this->_constant);

    for (const auto &div: _denominator) {
      base = base->divided_by(div->clone());
    }

    std::list<base_expression_type<C> *> result;

    result.push_back(base);

    for (auto it = std::begin(_numerator); it != std::end(_numerator); ++it) {
      base_expression_type<C> *new_it = (*it)->expand();

      result = new_it->multiply(result);

      delete new_it;
    }
    return new finite_sum_type<C>(result);
  }

  /**
   * @brief Get the numerator of the product
   *
   * @return the numerator of the product
   */
  base_expression_type<C> *get_numerator() const
  {
    base_expression_type<C> *prod = new finite_prod_type<C>(_constant);

    for (const auto term: _numerator) {
      prod = prod->multiply(term->clone());
    }

    return prod;
  }

  /**
   * @brief Get the denominator of the product
   *
   * @return the denominator of the product
   */
  base_expression_type<C> *get_denominator() const
  {
    if (_denominator.size() == 0) {
      return new constant_type<C>(1);
    }

    base_expression_type<C> *prod = new finite_prod_type<C>();

    for (const auto term: _denominator) {
      prod = prod->multiply(term->clone());
    }

    return prod;
  }

  /**
   * @brief Clone the current object
   *
   * @return A clone of the current object
   */
  base_expression_type<C> *clone() const
  {
    return new finite_prod_type<C>(*this);
  }

  /**
   * @brief Get the ids of the symbols in the expression
   *
   * @return The set of the symbol ids in the expression
   */
  std::set<SymbolIdType> get_symbol_ids() const
  {
    std::set<SymbolIdType> ids;

    for (auto &_list: {_numerator, _denominator}) {
      for (auto it = std::begin(_list); it != std::end(_list); ++it) {
        std::set<SymbolIdType> it_ids = (*it)->get_symbol_ids();

        ids.insert(std::begin(it_ids), std::end(it_ids));
      }
    }

    return ids;
  }

  /**
   * @brief Get the ids of the symbols in the expression
   *
   * @param max_cardinality is the maximum number of symbol ids to be returned
   * @return A set, whose cardinality is `max_cardinality` at most, of the
   * symbol identifiers in the expression
   */
  std::set<SymbolIdType>
  get_symbol_ids(const unsigned int max_cardinality) const
  {
    std::set<SymbolIdType> ids;

    unsigned int missing = max_cardinality;
    for (auto &_list: {_numerator, _denominator}) {
      for (auto it = std::begin(_list); it != std::end(_list); ++it) {
        std::set<SymbolIdType> it_ids = (*it)->get_symbol_ids(missing);

        ids.insert(std::begin(it_ids), std::end(it_ids));

        missing = max_cardinality - ids.size();
        if (missing == 0) {
          return ids;
        }
      }
    }

    return ids;
  }

  /**
   * @brief Test the presence of any symbol
   *
   * @return `true` if and only if the expression as one symbol at least
   */
  bool has_symbols() const
  {
    for (auto &_list: {_numerator, _denominator}) {
      for (auto it = std::begin(_list); it != std::end(_list); ++it) {
        if ((*it)->has_symbols()) {
          return true;
        }
      }
    }

    return false;
  }

  /**
   * @brief Replace symbol occurrences by using expressions
   *
   * This method replaces any occurrence of the symbols whose ids are in
   * the domain of `replacement` with a clone of the corresponding expressions.
   * This method updates this object.
   *
   * @param replacements associates symbol ids to their replacements
   * @return A pointer to the updated object
   */
  base_expression_type<C> *replace(
      const std::map<SymbolIdType, base_expression_type<C> *> &replacements)
  {
    base_expression_type<C> *result = new constant_type<C>(this->_constant);
    for (auto it = std::begin(_numerator); it != std::end(_numerator); ++it) {
      if (!(*it)->is_zero()) {
        result = result->multiply((*it)->replace(replacements));
      } else {
        delete *it;
      }
    }
    _numerator.clear();

    for (auto it = std::begin(_denominator); it != std::end(_denominator);
         ++it) {
      if (!(*it)->is_zero()) {
        result = result->divided_by((*it)->replace(replacements));
      } else {
        delete *it;
      }
    }
    _denominator.clear();

    delete this;

    return result;
  }

  /**
   * @brief Numerically evaluate the expression
   *
   * @return The numeric evaluation of the expression
   */
  C evaluate() const
  {
    C prod = this->_constant;
    for (auto it = std::begin(_numerator); it != std::end(_numerator); ++it) {
      prod *= (*it)->evaluate();
    }
    for (auto it = std::begin(_denominator); it != std::end(_denominator);
         ++it) {
      prod /= (*it)->evaluate();
    }

    return prod;
  }

  /**
   * @brief Apply a symbolic interpretation to an expression
   *
   * This method evaluates the expression using the symbol interpretation
   * provided as parameter.
   *
   * @param interpretation is an interpretation for the symbols
   * @return the expression obtained from the current object
   *         and `interpretation`
   */
  base_expression_type<C> *
  apply(const std::map<SymbolIdType, C> &interpretation) const
  {
    base_expression_type<C> *prod = new constant_type<C>(_constant);

    for (auto it = std::begin(_numerator); it != std::end(_numerator); ++it) {
      prod = prod->multiply((*it)->apply(interpretation));
    }
    for (auto it = std::begin(_denominator); it != std::end(_denominator);
         ++it) {
      prod = prod->divided_by((*it)->apply(interpretation));
    }

    return prod;
  }

  /**
   * @brief Get the coefficient of a term having a given degree
   *
   * @param symbol_id is the symbol id whose term coefficient is aimed
   * @param degree is the degree of the aimed term
   * @return The coefficient of the term in which the symbol having id
   *         `symbol_id` has degree `degree`
   */
  base_expression_type<C> *get_coeff(const SymbolIdType &symbol_id,
                                     const int &degree) const
  {
    C constant = this->_constant;

    try {
      for (auto it = std::begin(_denominator); it != std::begin(_denominator);
           ++it) {
        constant /= (*it)->evaluate();
      }
    } catch (std::runtime_error &) {
      SAPO_ERROR("non-constant denominator not supported", std::runtime_error);
    }

    base_expression_type<C> *res = new constant_type<C>(constant);

    int degree_counter = 0;
    for (auto it = std::begin(_numerator); it != std::end(_numerator); ++it) {
      switch ((*it)->type()) {
      case SYMBOL:
        if (((symbol_type<C> *)(*it))->get_id() == symbol_id) {
          ++degree_counter;
          if (degree_counter > degree) {
            delete res;

            return new constant_type<C>(0);
          }
        } else {
          res = res->multiply((*it)->clone());
        }
        break;
      case CONSTANT:
        res = res->multiply((*it)->clone());
        break;
      case FINITE_PROD:
      case FINITE_SUM:
      default:
        SAPO_ERROR("'get_coeff' call only admitted after 'expand' call",
                   std::runtime_error);
      }
    }

    if (degree_counter == degree) {
      return res;
    }

    delete res;

    return new constant_type<C>(0);
  }

  /**
   * @brief Get the coefficients of the polynomial expression
   *
   * @param symbol_id is the symbol id whose term coefficient is aimed
   * @return The coefficients of the polynomial expression on the
   *         symbol whose id is `symbol_id`
   */
  std::map<int, base_expression_type<C> *>
  get_coeffs(const SymbolIdType &symbol_id) const
  {
    if (_denominator.size() > 0) {
      SAPO_ERROR("non-constant denominator is not supported",
                 std::runtime_error);
    }

    std::map<int, base_expression_type<C> *> res;

    res[0] = new constant_type<C>(this->_constant);

    for (auto it = std::begin(_numerator); it != std::end(_numerator); ++it) {
      auto map_it = (*it)->get_coeffs(symbol_id);

      // if 0 is among the returned coefficients
      auto constant = map_it.find(0);
      if (constant != std::end(map_it) && !constant->second->has_symbols()
          && constant->second->evaluate() == 0) {

        // delete it
        map_it.erase(constant);
      }

      std::map<int, base_expression_type<C> *> next_res;
      for (auto res_coeff = std::begin(res); res_coeff != std::end(res);
           ++res_coeff) {
        for (auto it_coeff = std::begin(map_it); it_coeff != std::end(map_it);
             ++it_coeff) {
          const int degree = res_coeff->first + it_coeff->first;
          auto res_coeff_clone = res_coeff->second->clone();
          res_coeff_clone
              = res_coeff_clone->multiply(it_coeff->second->clone());
          auto coeff_degree = next_res.find(degree);
          if (coeff_degree == std::end(next_res)) {
            next_res[degree] = res_coeff_clone;
          } else {
            (coeff_degree->second)
                = (coeff_degree->second)->add(res_coeff_clone);
          }
        }
        delete res_coeff->second;
      }
      for (auto it_coeff = std::begin(map_it); it_coeff != std::end(map_it);
           ++it_coeff) {
        delete it_coeff->second;
      }
      res.clear();
      std::swap(res, next_res);
    }

    if (res.size() == 0) {
      res[0] = 0;
    }
    return res;
  }

  /**
   * @brief Get the degree of a symbol
   *
   * @param symbol_id is the id of the symbol whose degree is aimed
   * @return the degree of the parameter in the expression, i.e,
   *         the sum degree of the expression degree in the
   *          numerator minus the sum degree of the expression
   *          degree in the denominator.
   */
  int degree(const SymbolIdType &symbol_id) const
  {
    int total_degree = 0;

    for (auto it = std::begin(_numerator); it != std::end(_numerator); ++it) {
      total_degree += (*it)->degree(symbol_id);
    }

    for (auto it = std::begin(_denominator); it != std::end(_denominator);
         ++it) {
      total_degree -= (*it)->degree(symbol_id);
    }

    return total_degree;
  }

  /**
   * @brief Compute the derivative of an expression
   *
   * This method computes the derivative of the current expression
   * with respect to a symbol provided as the parameter.
   *
   * @param symbol_id is the id of the symbol with respect to the
   *               derivative must be computed
   * @return the derivative of the current expression with respect
   *         to the symbol associated to `symbol_id`
   */
  base_expression_type<C> *
  get_derivative_wrt(const SymbolIdType &symbol_id) const
  {
    // let be c*f(symbol)/g(symbol) the expression where c is a constant
    // and f and g are products of expressions in symbol

    base_expression_type<C> *res = new constant_type<C>(_constant);

    // evaluate the derivative of f(symbol)
    base_expression_type<C> *num_dev
        = _prod_list_derivative(_numerator, symbol_id);

    // if the denominator is empty, i.e., g(symbol) == 1
    if (_denominator.size() == 0) {

      // then the derivative is c*f'(symbol)
      return res->multiply(num_dev);
    }

    // otherwise, compute the derivative of f(symbol)
    base_expression_type<C> *denom_dev
        = _prod_list_derivative(_denominator, symbol_id);

    // evaluate g'*f
    for (const auto &term: _numerator) {
      denom_dev = denom_dev->multiply(term->clone());
    }

    // evaluate f'*g
    for (const auto &term: _denominator) {
      num_dev = num_dev->multiply(term->clone());
    }

    // compute c*( g'*f + f'*g )
    res = res->multiply(num_dev->subtract(denom_dev));

    // divide everything by g^2
    for (const auto &term: _denominator) {
      res = res->divided_by(term->clone());
      res = res->divided_by(term->clone());
    }

    return res;
  }

  /**
   * @brief Print the expression in an output stream
   *
   * @param os is the output stream in which the expression must be printed
   */
  void print(std::ostream &os) const
  {
    C constant = _constant;
    if (constant < 0) {
      os << "-";
      constant = -constant;
    }

    if (_numerator.size() > 0 || constant != 1) {
      print_list(os, _numerator, constant);
    }
    if (_denominator.size() > 0) {
      os << "/";
      print_list(os, _denominator);
    }
  }

  ~finite_prod_type()
  {
    for (auto &_list: {_numerator, _denominator}) {
      for (auto it = std::begin(_list); it != std::end(_list); ++it) {
        delete *it;
      }
    }

    _numerator.clear();
    _denominator.clear();
  }

  template<typename T>
  friend class base_expression_type;

  template<typename T>
  friend class finite_sum_type;

  template<typename T>
  friend struct less;
};

/**
 * @brief A class to represent symbols
 *
 * @tparam C is the type of expression numeric constants
 */
template<typename C>
class symbol_type;

} // end of low_level namespace

/**
 * @brief Exception raised when a symbol is evaluated
 */
class symbol_evaluation_error : public std::exception
{
private:
  std::string _symbol_name; //!< The name of the evaluated symbol
public:
  /**
   * @brief Constructor
   *
   * @tparam C is the type of expression numeric constants
   * @param symbol is the symbol that has been evaluated
   */
  template<typename C>
  explicit symbol_evaluation_error(const Symbol<C> &symbol):
      _symbol_name(symbol.get_name())
  {
  }

  /**
   * @brief Constructor
   *
   * @tparam C is the type of expression numeric constants
   * @param symbol is the symbol that has been evaluated
   */
  template<typename C>
  explicit symbol_evaluation_error(const low_level::symbol_type<C> &symbol):
      _symbol_name(Symbol<C>::get_symbol_name(symbol.get_id()))
  {
  }

  /**
   * @brief Get the name of the evaluated symbol
   *
   * @return the name of the evaluated symbol
   */
  inline const std::string &get_symbol_name() const noexcept
  {
    return _symbol_name;
  }

  /**
   * @brief Destroy the symbol evaluation error
   */
  virtual ~symbol_evaluation_error() noexcept {}

  /**
   * @brief Return the error message
   *
   * @return the error message
   */
  virtual const char *what() const noexcept
  {
    return "Symbols cannot be evaluated";
  }
};

namespace low_level
{

/**
 * @brief A class to represent symbols
 *
 * @tparam C is the type of numeric constants
 */
template<typename C>
class symbol_type : public base_expression_type<C>
{
public:
  typedef typename Symbol<C>::SymbolIdType SymbolIdType;

protected:
  SymbolIdType _id; //!< The symbol id.

public:
  /**
   * @brief Build a new symbol
   *
   * @param id is the id of the new symbol
   */
  symbol_type(const SymbolIdType id): base_expression_type<C>(), _id(id) {}

  /**
   * @brief Copy constructor
   *
   * @param orig is the model for the new symbol
   */
  symbol_type(const symbol_type<C> &orig):
      base_expression_type<C>(), _id(orig._id)
  {
  }

  /**
   * @brief Get the id of the symbol
   *
   * @return The id of the symbol
   */
  const SymbolIdType &get_id() const
  {
    return _id;
  }

  /**
   * @brief Get the expression type
   *
   * @return The value SYMBOL
   */
  ExpressionType type() const
  {
    return SYMBOL;
  }

  /**
   * @brief Get the rational polynomial representation of an expression
   *
   * @return the rational polynomial representation of this expression
   */
  base_expression_type<C> *get_rational_form() const
  {
    return this->clone();
  }

  /**
   * @brief Establish whether is a polynomial expression
   *
   * @return `true` if and only if the expression is a polynomial
   */
  bool is_a_polynomial() const
  {
    return true;
  }

  /**
   * @brief Establish whether is in rational polynomial form
   *
   * @return `true` if and only if the expression is either
   *        polynomial or the ratio between two polynomials
   */
  bool is_in_rational_form() const
  {
    return true;
  }

  /**
   * @brief Complement the expression
   *
   * @return the complement expression
   */
  base_expression_type<C> *complement()
  {
    return this->multiply(new constant_type<C>(static_cast<C>(-1)));
  }

  /**
   * @brief Clone the current object
   *
   * @return A clone of the current object
   */
  base_expression_type<C> *clone() const
  {
    return new symbol_type(*this);
  }

  /**
   * @brief Get the ids of the symbols in the expression
   *
   * @return The set of the symbol ids in the expression
   */
  std::set<SymbolIdType> get_symbol_ids() const
  {
    std::set<SymbolIdType> ids;

    ids.insert(_id);

    return ids;
  }

  /**
   * @brief Get the ids of the symbols in the expression
   *
   * @param max_cardinality is the maximum number of symbol ids to be returned
   * @return A set, whose cardinality is `max_cardinality` at most, of the
   * symbol identifiers in the expression
   */
  std::set<SymbolIdType>
  get_symbol_ids(const unsigned int max_cardinality) const
  {
    if (max_cardinality == 0) {
      return std::set<SymbolIdType>();
    }

    return std::set<SymbolIdType>({_id});
  }

  /**
   * @brief Test the presence of any symbol
   *
   * @return `true` if and only if the expression as one symbol at least
   */
  bool has_symbols() const
  {
    return true;
  }

  /**
   * @brief Replace symbol occurrences by using expressions
   *
   * This method replaces any occurrence of the symbols whose ids are in
   * the domain of `replacement` with a clone of the corresponding expressions.
   * This method updates this object.
   *
   * @param replacements associates symbol ids to their replacements
   * @return A pointer to the updated object
   */
  base_expression_type<C> *replace(
      const std::map<SymbolIdType, base_expression_type<C> *> &replacements)
  {
    auto it = replacements.find(_id);

    if (it == std::end(replacements)) {
      return this;
    }

    delete this;

    return it->second->clone();
  }

  /**
   * @brief Numerically evaluate the expression
   *
   * @return The numeric evaluation of the expression
   */
  C evaluate() const
  {
    throw symbol_evaluation_error(*this);
  }

  /**
   * @brief Apply a symbolic interpretation to an expression
   *
   * This method evaluates the expression using the symbol interpretation
   * provided as parameter.
   *
   * @param interpretation is an interpretation for the symbols
   * @return the expression obtained from the current object
   *         and `interpretation`
   */
  base_expression_type<C> *
  apply(const std::map<SymbolIdType, C> &interpretation) const
  {
    auto found = interpretation.find(_id);

    if (found == std::end(interpretation)) {
      return this->clone();
    }

    return new constant_type<C>(found->second);
  }

  /**
   * @brief Get the coefficient of a term having a given degree
   *
   * @param symbol_id is the symbol id whose term coefficient is aimed
   * @param degree is the degree of the aimed term
   * @return The coefficient of the term in which the symbol having id
   *         `symbol_id` has degree `degree`
   */
  base_expression_type<C> *get_coeff(const SymbolIdType &symbol_id,
                                     const int &degree) const
  {
    if (symbol_id == _id) {
      return new constant_type<C>(degree == 1 ? 1 : 0);
    } else {
      if (degree == 0) {
        return this->clone();
      } else {
        return new constant_type<C>(0);
      }
    }
  }

  /**
   * @brief Get the coefficients of the polynomial expression
   *
   * @param symbol_id is the symbol id whose term coefficient is aimed
   * @return The coefficients of the polynomial expression on the
   *         symbol whose id is `symbol_id`
   */
  std::map<int, base_expression_type<C> *>
  get_coeffs(const SymbolIdType &symbol_id) const
  {
    std::map<int, base_expression_type<C> *> res;

    if (symbol_id == _id) {
      res[1] = new constant_type<C>(1);
    } else {
      res[0] = this->clone();
    }

    return res;
  }

  /**
   * @brief Get the degree of a symbol
   *
   * @param symbol_id is the id of the symbol whose degree is aimed
   * @return the degree of the parameter in the expression, i.e,
   *         1, if the current symbol is the parameter; 0, otherwise
   */
  int degree(const SymbolIdType &symbol_id) const
  {
    return ((((const symbol_type<C> *)this)->_id == symbol_id) ? 1 : 0);
  }

  /**
   * @brief Compute the derivative of an expression
   *
   * This method computes the derivative of the current expression
   * with respect to a symbol provided as the parameter.
   *
   * @param symbol_id is the id of the symbol with respect to the
   *               derivative must be computed
   * @return the derivative of the current expression with respect
   *         to the symbol associated to `symbol_id`
   */
  base_expression_type<C> *
  get_derivative_wrt(const SymbolIdType &symbol_id) const
  {
    if (symbol_id == this->_id) {
      return new constant_type<C>(1);
    }

    return new constant_type<C>(0);
  }

  /**
   * @brief Print the expression in an output stream
   *
   * @param os is the output stream in which the expression must be printed
   */
  void print(std::ostream &os) const
  {
    os << Symbol<C>::get_symbol_name(_id);
  }
};

}

template<typename C>
Expression<C>::Expression(low_level::base_expression_type<C> *ex): _ex(ex)
{
}

template<typename C>
Expression<C>::Expression():
    _ex(new low_level::constant_type<C>(static_cast<C>(0)))
{
}

template<typename C>
Expression<C>::Expression(const int value):
    _ex(new low_level::constant_type<C>(static_cast<C>(value)))
{
}

template<typename C>
Expression<C>::Expression(const C value):
    _ex(new low_level::constant_type<C>(value))
{
}

template<typename C>
Expression<C>::Expression(const Expression<C> &orig):
    _ex((orig._ex == nullptr ? nullptr : orig._ex->clone()))
{
}

template<typename C>
Expression<C> &
Expression<C>::replace(const Expression<C>::replacement_type &replacement)
{
  if (this->_ex == nullptr) {
    return *this;
  }
  std::map<typename Symbol<C>::SymbolIdType,
           low_level::base_expression_type<C> *>
      base_repl;

  for (auto it = std::begin(replacement); it != std::end(replacement); ++it) {
    base_repl[it->first.get_id()] = it->second._ex;
  }

  this->_ex = this->_ex->replace(base_repl);

  return *this;
}

template<typename C>
Expression<C> Expression<C>::get_rational_form() const
{
  return Expression<C>(_ex->get_rational_form());
}

template<typename C>
bool Expression<C>::is_a_polynomial() const
{
  return _ex->is_a_polynomial();
}

template<typename C>
bool Expression<C>::is_in_rational_form() const
{
  return _ex->is_in_rational_form();
}

template<typename C>
bool Expression<C>::is_a_constant() const
{
  return _ex->type() == CONSTANT;
}

template<typename C>
Expression<C> Expression<C>::get_numerator() const
{
  if (_ex->type() == FINITE_PROD) {
    const auto prod = static_cast<const low_level::finite_prod_type<C> *>(_ex);
    return Expression<C>(prod->get_numerator());
  }

  return Expression<C>(_ex->clone());
}

template<typename C>
Expression<C> Expression<C>::get_derivative_wrt(const Symbol<C> &symbol) const
{
  return Expression<C>(_ex->get_derivative_wrt(symbol.get_id()));
}

template<typename C>
Expression<C> Expression<C>::get_denominator() const
{
  if (_ex->type() == FINITE_PROD) {
    const auto prod = static_cast<const low_level::finite_prod_type<C> *>(_ex);
    return Expression<C>(prod->get_denominator());
  }

  return Expression<C>(1);
}

template<typename C>
Expression<C> Expression<C>::apply(
    const Expression<C>::interpretation_type &interpretation) const
{
  if (this->_ex == nullptr) {
    return 0;
  }
  std::map<typename Symbol<C>::SymbolIdType, C> base_interpretation;

  for (auto it = std::begin(interpretation); it != std::end(interpretation);
       ++it) {
    base_interpretation[it->first.get_id()] = it->second;
  }

  return Expression<C>(_ex->apply(base_interpretation));
}

template<typename C>
Expression<C> &Expression<C>::expand()
{
  low_level::base_expression_type<C> *new_ex = _ex->expand();

  delete _ex;

  _ex = new_ex;

  return *this;
}

template<typename C>
Expression<C>::~Expression()
{
  if (_ex != nullptr) {
    delete _ex;
    _ex = nullptr;
  }
}

template<typename C>
constexpr Expression<C> pow(const Expression<C> &base, unsigned int exp)
{
  if (exp == 0) {
    return 1;
  }
  Expression<C> res = pow(base, exp / 2);
  res *= res;
  if (exp % 2 == 0) {
    return res;
  }

  return res * base;
}

template<typename C>
Expression<C> &Expression<C>::operator=(const Expression<C> &rhs)
{
  delete _ex;

  _ex = rhs._ex->clone();

  return *this;
}

template<typename C>
Expression<C> &Expression<C>::operator=(Expression<C> &&rhs)
{
  std::swap(_ex, rhs._ex);

  return *this;
}

template<typename C>
const Expression<C> &Expression<C>::operator+=(const Expression<C> &rhs)
{
  low_level::base_expression_type<C> *rhs_ex = rhs._ex->clone();

  _ex = _ex->add(rhs_ex);

  return *this;
}

template<typename C>
const Expression<C> &Expression<C>::operator-=(const Expression<C> &rhs)
{
  low_level::base_expression_type<C> *rhs_ex = rhs._ex->clone();

  _ex = _ex->subtract(rhs_ex);

  return *this;
}

template<typename C>
const Expression<C> &Expression<C>::operator*=(const Expression<C> &rhs)
{
  low_level::base_expression_type<C> *rhs_ex = rhs._ex->clone();

  _ex = _ex->multiply(rhs_ex);

  return *this;
}

template<typename C>
const Expression<C> &Expression<C>::operator/=(const Expression<C> &rhs)
{
  low_level::base_expression_type<C> *rhs_ex = rhs._ex->clone();

  _ex = _ex->divided_by(rhs_ex);

  return *this;
}

template<typename C>
const Expression<C> &Expression<C>::operator+=(Expression<C> &&rhs)
{
  low_level::base_expression_type<C> *rhs_ex = rhs._ex;
  rhs._ex = nullptr;

  _ex = _ex->add(rhs_ex);

  return *this;
}

template<typename C>
const Expression<C> &Expression<C>::operator-=(Expression<C> &&rhs)
{
  low_level::base_expression_type<C> *rhs_ex = rhs._ex;
  rhs._ex = nullptr;

  _ex = _ex->subtract(rhs_ex);

  return *this;
}

template<typename C>
const Expression<C> &Expression<C>::operator*=(Expression<C> &&rhs)
{
  low_level::base_expression_type<C> *rhs_ex = rhs._ex;
  rhs._ex = nullptr;

  _ex = _ex->multiply(rhs_ex);

  return *this;
}

template<typename C>
const Expression<C> &Expression<C>::operator/=(Expression<C> &&rhs)
{
  low_level::base_expression_type<C> *rhs_ex = rhs._ex;
  rhs._ex = nullptr;

  _ex = _ex->divided_by(rhs_ex);

  return *this;
}

/**
 * @brief Simplify an expression
 *
 * This function produces an expression which is a constant
 * or in which all the coefficients of the polynomial
 * representation on one of its symbol have been simplified.
 *
 * @tparam C is the constant value type
 * @param exp is the expression to be simplified
 * @return a simplified expression
 */
template<typename C>
Expression<C> simplify(const Expression<C> &exp)
{
  // get at most one of the symbols in the expression
  std::set<Symbol<C>> symbols = exp.get_symbols(1);

  // if the expression has no symbols
  if (symbols.size() == 0) {

    // it is already simplified
    return exp;
  }

  // otherwise, get the polynomial representation of the
  // expression on the found symbol
  int p = 0;
  Expression<C> simpl_ex(0), power(1);

  Symbol<C> x = *std::begin(symbols);
  auto x_coeffs = exp.get_coeffs(x);

  // simplify all the coefficients
  for (auto c_it = std::begin(x_coeffs); c_it != std::end(x_coeffs); ++c_it) {
    c_it->second = simplify(c_it->second);

    if (c_it->second.has_symbols() || c_it->second.evaluate() != 0) {
      while (p < c_it->first) {
        ++p;
        power *= x;
      }
      simpl_ex += power * c_it->second;
    }
  }

  return simpl_ex;
}

/**
 * @brief Check whether the two expressions are equivalent
 *
 * @tparam C is the numeric type of constants
 * @param lhs is an expression
 * @param rhs is an expression
 * @return `true` if and only if the two expressions are equivalent
 */
template<typename C>
bool are_equivalent(const Expression<C> &lhs, const Expression<C> &rhs)
{
  Expression<C> diff = simplify(lhs - rhs);

  return !(diff.has_symbols() || diff.evaluate() != 0);
}

/**
 * @brief Check whether an expression is semantically equivalent to a constant
 *
 * @tparam C is the numeric type of constants
 * @param lhs is an expression
 * @param rhs is a constant value
 * @return `true` if and only if the expression is equivalent to the
 *         constant value
 */
template<typename C>
inline bool are_equivalent(const Expression<C> &lhs, const C &rhs)
{
  return are_equivalent(lhs, Expression<C>(rhs));
}

/**
 * @brief Check whether an expression is semantically equivalent to a constant
 *
 * @tparam C is the numeric type of constants
 * @param lhs is a constant value
 * @param rhs is an expression
 * @return `true` if and only if the constant value is equivalent to the
 *         expression
 */
template<typename C>
inline bool are_equivalent(const C &lhs, const Expression<C> &rhs)
{
  return are_equivalent(rhs, lhs);
}

/**
 * @brief Check whether an expression is semantically equivalent to a constant
 *
 * @tparam C is the numeric type of constants
 * @param lhs is an expression
 * @param rhs is a constant value
 * @return `true` if and only if the expression is equivalent to the
 *         constant value
 */
template<typename C>
inline bool operator==(const Expression<C> &lhs, const C rhs)
{
  return are_equivalent(lhs, rhs);
}

/**
 * @brief Check whether an expression is semantically equivalent to a constant
 *
 * @tparam C is the numeric type of constants
 * @param lhs is an expression
 * @param rhs is a constant value
 * @return `true` if and only if the first parameter is equivalent to the
 *         second one
 */
template<typename C>
inline bool operator==(const Expression<C> &lhs, const int rhs)
{
  return lhs == static_cast<C>(rhs);
}

/**
 * @brief Check whether an expression is semantically equivalent to a constant
 *
 * @tparam C is the numeric type of constants
 * @param lhs is a constant value
 * @param rhs is an expression
 * @return `true` if and only if the first parameter is equivalent to the
 *         second one
 */
template<typename C>
inline bool operator==(const C lhs, const Expression<C> &rhs)
{
  return rhs == lhs;
}

/**
 * @brief Check whether an expression is semantically equivalent to a constant
 *
 * @tparam C is the numeric type of constants
 * @param lhs is a constant value
 * @param rhs is an expression
 * @return `true` if and only if the first parameter is equivalent to the
 *         second one
 */
template<typename C>
inline bool operator==(const int lhs, const Expression<C> &rhs)
{
  return rhs == static_cast<C>(lhs);
}

/**
 * @brief Check whether an expression is not equivalent to a constant value
 *
 * @tparam C is the numeric type of constants
 * @param lhs is an expression
 * @param rhs is a constant value
 * @return `true` if and only if the first parameter is not equivalent to the
 *         second one
 */
template<typename C>
inline bool operator!=(const Expression<C> &lhs, const C rhs)
{
  return !(lhs == rhs);
}

/**
 * @brief Check whether an expression is not equivalent to a constant value
 *
 * @tparam C is the numeric type of constants
 * @param lhs is an expression
 * @param rhs is a constant value
 * @return `true` if and only if the first parameter is not equivalent to the
 *         second one
 */
template<typename C>
inline bool operator!=(const Expression<C> &lhs, const int rhs)
{
  return !(lhs == rhs);
}

/**
 * @brief Check whether an expression is not equivalent to a constant value
 *
 * @tparam C is the numeric type of constants
 * @param lhs is a constant value
 * @param rhs is an expression
 * @return `true` if and only if the first parameter is not equivalent to the
 *         second one
 */
template<typename C>
inline bool operator!=(const C lhs, const Expression<C> &rhs)
{
  return !(lhs == rhs);
}

/**
 * @brief Check whether an expression is not equivalent to a constant value
 *
 * @tparam C is the numeric type of constants
 * @param lhs is a constant value
 * @param rhs is an expression
 * @return `true` if and only if the first parameter is not equivalent to the
 *         second one
 */
template<typename C>
inline bool operator!=(const int lhs, const Expression<C> &rhs)
{
  return !(lhs == rhs);
}

/**
 * @brief Sum two expressions
 *
 * This method builds an expression that represents the sum of the two
 * expressions passed as parameters.
 *
 * @tparam C is the type of numeric constants
 * @param lhs is an expression
 * @param rhs is an expression
 * @return An expression that represents the sum `lhs + rhs`
 */
template<typename C>
Expression<C> operator+(const Expression<C> &lhs, const Expression<C> &rhs)
{
  low_level::base_expression_type<C> *lhs_ex = lhs._ex->clone();

  low_level::base_expression_type<C> *rhs_ex = rhs._ex->clone();

  return Expression<C>(lhs_ex->add(rhs_ex));
}

/**
 * @brief Sum two expressions
 *
 * This method builds an expression that represents the sum of the two
 * expressions passed as parameters.
 *
 * @tparam C is the type of numeric constants
 * @param lhs is an expression
 * @param rhs is an expression
 * @return An expression that represents the sum `lhs + rhs`
 */
template<typename C>
Expression<C> operator+(Expression<C> &&lhs, const Expression<C> &rhs)
{
  low_level::base_expression_type<C> *lhs_ex = lhs._ex;
  lhs._ex = nullptr;

  low_level::base_expression_type<C> *rhs_ex = rhs._ex->clone();

  return Expression<C>(lhs_ex->add(rhs_ex));
}

/**
 * @brief Sum two expressions
 *
 * This method builds an expression that represents the sum of the two
 * expressions passed as parameters.
 *
 * @tparam C is the type of numeric constants
 * @param lhs is an expression
 * @param rhs is an expression
 * @return An expression that represents the sum `lhs + rhs`
 */
template<typename C>
Expression<C> operator+(const Expression<C> &lhs, Expression<C> &&rhs)
{
  low_level::base_expression_type<C> *lhs_ex = lhs._ex->clone();

  low_level::base_expression_type<C> *rhs_ex = rhs._ex;
  rhs._ex = nullptr;

  return Expression<C>(lhs_ex->add(rhs_ex));
}

/**
 * @brief Sum two expressions
 *
 * This method builds an expression that represents the sum of the two
 * expressions passed as parameters.
 *
 * @tparam C is the type of numeric constants
 * @param lhs is an expression
 * @param rhs is an expression
 * @return An expression that represents the sum `lhs + rhs`
 */
template<typename C>
Expression<C> operator+(Expression<C> &&lhs, Expression<C> &&rhs)
{
  low_level::base_expression_type<C> *lhs_ex = lhs._ex;
  lhs._ex = nullptr;

  low_level::base_expression_type<C> *rhs_ex = rhs._ex;
  rhs._ex = nullptr;

  return Expression<C>(lhs_ex->add(rhs_ex));
}

/**
 * @brief Subtract two expressions
 *
 * This method builds an expression that represents the subtraction between the
 * two expressions passed as parameters.
 *
 * @tparam C is the type of numeric constants
 * @param lhs is an expression
 * @param rhs is an expression
 * @return An expression that represents the subtraction `lhs - rhs`
 */
template<typename C>
Expression<C> operator-(const Expression<C> &lhs, const Expression<C> &rhs)
{
  low_level::base_expression_type<C> *lhs_ex = lhs._ex->clone();

  low_level::base_expression_type<C> *rhs_ex = rhs._ex->clone();

  return Expression<C>(lhs_ex->subtract(rhs_ex));
}

/**
 * @brief Subtract two expressions
 *
 * This method builds an expression that represents the subtraction between the
 * two expressions passed as parameters.
 *
 * @tparam C is the type of numeric constants
 * @param lhs is an expression
 * @param rhs is an expression
 * @return An expression that represents the subtraction `lhs - rhs`
 */
template<typename C>
Expression<C> operator-(Expression<C> &&lhs, const Expression<C> &rhs)
{
  low_level::base_expression_type<C> *lhs_ex = lhs._ex;
  lhs._ex = nullptr;

  low_level::base_expression_type<C> *rhs_ex = rhs._ex->clone();

  return Expression<C>(lhs_ex->subtract(rhs_ex));
}

/**
 * @brief Subtract two expressions
 *
 * This method builds an expression that represents the subtraction between the
 * two expressions passed as parameters.
 *
 * @tparam C is the type of numeric constants
 * @param lhs is an expression
 * @param rhs is an expression
 * @return An expression that represents the subtraction `lhs - rhs`
 */
template<typename C>
Expression<C> operator-(const Expression<C> &lhs, Expression<C> &&rhs)
{
  low_level::base_expression_type<C> *lhs_ex = lhs._ex->clone();

  low_level::base_expression_type<C> *rhs_ex = rhs._ex;
  rhs._ex = nullptr;

  return Expression<C>(lhs_ex->subtract(rhs_ex));
}

/**
 * @brief Subtract two expressions
 *
 * This method builds an expression that represents the subtraction between the
 * two expressions passed as parameters.
 *
 * @tparam C is the type of numeric constants
 * @param lhs is an expression
 * @param rhs is an expression
 * @return An expression that represents the subtraction `lhs - rhs`
 */
template<typename C>
Expression<C> operator-(Expression<C> &&lhs, Expression<C> &&rhs)
{
  low_level::base_expression_type<C> *lhs_ex = lhs._ex;
  lhs._ex = nullptr;

  low_level::base_expression_type<C> *rhs_ex = rhs._ex;
  rhs._ex = nullptr;

  return Expression<C>(lhs_ex->subtract(rhs_ex));
}

/**
 * @brief Multiply two expressions
 *
 * This method builds an expression that represents the multiplication of the
 * two expressions passed as parameters.
 *
 * @tparam C is the type of numeric constants
 * @param lhs is an expression
 * @param rhs is an expression
 * @return An expression that represents the multiplication `lhs * rhs`
 */
template<typename C>
Expression<C> operator*(const Expression<C> &lhs, const Expression<C> &rhs)
{
  low_level::base_expression_type<C> *lhs_ex = lhs._ex->clone();

  low_level::base_expression_type<C> *rhs_ex = rhs._ex->clone();

  return Expression<C>(lhs_ex->multiply(rhs_ex));
}

/**
 * @brief Multiply two expressions
 *
 * This method builds an expression that represents the multiplication of the
 * two expressions passed as parameters.
 *
 * @tparam C is the type of numeric constants
 * @param lhs is an expression
 * @param rhs is an expression
 * @return An expression that represents the multiplication `lhs * rhs`
 */
template<typename C>
Expression<C> operator*(Expression<C> &&lhs, const Expression<C> &rhs)
{
  low_level::base_expression_type<C> *lhs_ex = lhs._ex;
  lhs._ex = nullptr;

  low_level::base_expression_type<C> *rhs_ex = rhs._ex->clone();

  return Expression<C>(lhs_ex->multiply(rhs_ex));
}

/**
 * @brief Multiply two expressions
 *
 * This method builds an expression that represents the multiplication of the
 * two expressions passed as parameters.
 *
 * @tparam C is the type of numeric constants
 * @param lhs is an expression
 * @param rhs is an expression
 * @return An expression that represents the multiplication `lhs * rhs`
 */
template<typename C>
Expression<C> operator*(const Expression<C> &lhs, Expression<C> &&rhs)
{
  low_level::base_expression_type<C> *lhs_ex = lhs._ex->clone();

  low_level::base_expression_type<C> *rhs_ex = rhs._ex;
  rhs._ex = nullptr;

  return Expression<C>(lhs_ex->multiply(rhs_ex));
}

/**
 * @brief Multiply two expressions
 *
 * This method builds an expression that represents the multiplication of the
 * two expressions passed as parameters.
 *
 * @tparam C is the type of numeric constants
 * @param lhs is an expression
 * @param rhs is an expression
 * @return An expression that represents the multiplication `lhs * rhs`
 */
template<typename C>
Expression<C> operator*(Expression<C> &&lhs, Expression<C> &&rhs)
{
  low_level::base_expression_type<C> *lhs_ex = lhs._ex;
  lhs._ex = nullptr;

  low_level::base_expression_type<C> *rhs_ex = rhs._ex;
  rhs._ex = nullptr;

  return Expression<C>(lhs_ex->multiply(rhs_ex));
}

/**
 * @brief Divide two expressions
 *
 * This method builds an expression that represents the division between the
 * two expressions passed as parameters.
 *
 * @tparam C is the type of numeric constants
 * @param lhs is an expression
 * @param rhs is an expression
 * @return An expression that represents the division `lhs / rhs`
 */
template<typename C>
Expression<C> operator/(const Expression<C> &lhs, const Expression<C> &rhs)
{
  low_level::base_expression_type<C> *lhs_ex = lhs._ex->clone();

  low_level::base_expression_type<C> *rhs_ex = rhs._ex->clone();

  return Expression<C>(lhs_ex->divided_by(rhs_ex));
}

/**
 * @brief Divide two expressions
 *
 * This method builds an expression that represents the division between the
 * two expressions passed as parameters.
 *
 * @tparam C is the type of numeric constants
 * @param lhs is an expression
 * @param rhs is an expression
 * @return An expression that represents the division `lhs / rhs`
 */
template<typename C>
Expression<C> operator/(Expression<C> &&lhs, const Expression<C> &rhs)
{
  low_level::base_expression_type<C> *lhs_ex = lhs._ex;
  lhs._ex = nullptr;

  low_level::base_expression_type<C> *rhs_ex = rhs._ex->clone();

  return Expression<C>(lhs_ex->divided_by(rhs_ex));
}

/**
 * @brief Divide two expressions
 *
 * This method builds an expression that represents the division between the
 * two expressions passed as parameters.
 *
 * @tparam C is the type of numeric constants
 * @param lhs is an expression
 * @param rhs is an expression
 * @return An expression that represents the division `lhs / rhs`
 */
template<typename C>
Expression<C> operator/(const Expression<C> &lhs, Expression<C> &&rhs)
{
  low_level::base_expression_type<C> *lhs_ex = lhs._ex->clone();

  low_level::base_expression_type<C> *rhs_ex = rhs._ex;
  rhs._ex = nullptr;

  return Expression<C>(lhs_ex->divided_by(rhs_ex));
}

/**
 * @brief Divide two expressions
 *
 * This method builds an expression that represents the division between the
 * two expressions passed as parameters.
 *
 * @tparam C is the type of numeric constants
 * @param lhs is an expression
 * @param rhs is an expression
 * @return An expression that represents the division `lhs / rhs`
 */
template<typename C>
Expression<C> operator/(Expression<C> &&lhs, Expression<C> &&rhs)
{
  low_level::base_expression_type<C> *lhs_ex = lhs._ex;
  lhs._ex = nullptr;

  low_level::base_expression_type<C> *rhs_ex = rhs._ex;
  rhs._ex = nullptr;

  return Expression<C>(lhs_ex->divided_by(rhs_ex));
}

template<typename C>
Symbol<C>::Symbol(): Expression<C>(nullptr)
{
}

template<typename C>
Symbol<C>::Symbol(const char *name): Symbol<C>(std::string(name))
{
}

template<typename C>
Symbol<C>::Symbol(const std::string &name): Expression<C>(nullptr)
{
#ifdef WITH_THREADS
  std::unique_lock<std::mutex> lock(_mutex);
#endif // WITH_THREADS

  if (name.size() == 0) {
    SAPO_ERROR("the name must be non-empty", std::domain_error);
  }

  auto found_symbol = _declared_symbols.find(name);

  if (found_symbol == std::end(_declared_symbols)) {
    SymbolIdType symbol_id = _symbol_names.size();
    _symbol_names.push_back(name);
    _declared_symbols[name] = symbol_id;
    this->_ex = new low_level::symbol_type<C>(symbol_id);
  } else {
    this->_ex = new low_level::symbol_type<C>(found_symbol->second);
  }
}

template<typename C>
Symbol<C>::Symbol(const Symbol<C> &orig):
    Expression<C>((orig._ex == nullptr ? nullptr : orig._ex->clone()))
{
}

template<typename C>
Symbol<C> &Symbol<C>::operator=(const Symbol<C> &orig)
{
  Expression<C>::operator=(orig);

  return *this;
}

template<typename C>
Symbol<C> &Symbol<C>::operator=(Symbol<C> &&orig)
{
  Expression<C>::operator=(std::move(orig));

  return *this;
}

}

/*! @} End of SymbolicAlgebra group */

namespace std
{
/**
 * @brief Swap the content of two Expression objects
 *
 * @tparam C is the type of numeric constants
 * @param a is an expression
 * @param b is an expression
 */
template<typename C>
void swap(SymbolicAlgebra::Expression<C> &a, SymbolicAlgebra::Expression<C> &b)
{
  std::swap(a._ex, b._ex);
}

/**
 * @brief Write an expression in an output stream
 *
 * @tparam C is the type of numeric constants
 * @param os is the output stream
 * @param ex is the expression to be printed
 * @return the output stream
 */
template<typename C>
std::ostream &operator<<(std::ostream &os,
                         const SymbolicAlgebra::Expression<C> &ex)
{
  if (ex._ex != nullptr) {
    ex._ex->print(os);
  }

  return os;
}
} // namespace std

#endif // SYMBOLIC_ALGEBRA_H_
