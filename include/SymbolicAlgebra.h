#ifndef SYMBOLIC_ALGEBRA_H_
#define SYMBOLIC_ALGEBRA_H_

#include <iostream>
#include <list>
#include <map>
#include <set>
#include <vector>
#include <type_traits>

#include <gmpxx.h>


/*!
 *  \addtogroup SymbolicAlgebra
 *  @{
 */

//! Symbolic algebra namespace
namespace SymbolicAlgebra
{

typedef enum { CONSTANT, SYMBOL, FINITE_SUM, FINITE_PROD } ExpressionType;

template<typename C>
class _constant_type;

template<typename C>
class _symbol_type;

template<typename C>
class _finite_sum_type;

template<typename C>
class _finite_prod_type;

template<typename C>
class _base_expression_type;

template<typename C>
class Expression;

/**
 * @brief A class to represent algebraic symbols
 *
 * This class represents algebraic symbols. It univocally associates any symbol
 * name to an id by using two static members: the `_declared_symbols` map and
 * the `_symbol_names` vector. The former maps any symbol name to an id; the
 * latter connects an id to the corresponding symbol name.
 *
 * @tparam C is the type of numeric constants.
 */
template<typename C = double>
class Symbol : public Expression<C>
{
public:
  typedef unsigned int SymbolIdType; //!< The type of symbol identificators

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
   * @brief Build a new Symbol object.
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
   * @brief Build a new Symbol object.
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
   * @param orig is the Symbol object to be copied.
   */
  Symbol(const Symbol<C> &orig);

  /**
   * @brief Get the Symbol id.
   *
   * @return The id of the represented symbol.
   */
  const SymbolIdType &get_id() const
  {
    return ((_symbol_type<C> *)(this->_ex))->get_id();
  }

  /**
   * @brief Get the name associated to a symbol id
   *
   * @param id is the query symbol id
   * @return the name corresponding to the query symbol id.
   */
  static const std::string &get_symbol_name(const SymbolIdType &id)
  {
    return _symbol_names[id];
  }
};

template<typename C>
std::map<std::string, typename Symbol<C>::SymbolIdType>
    Symbol<C>::_declared_symbols;
template<typename C>
std::vector<std::string> Symbol<C>::_symbol_names;

/**
 * @brief A class to represent algebraic symbolic expressions.
 *
 * This class is a wrapper class: the expression representation is achieved by
 * the `_base_expression_type` class and its hierarchy.
 *
 * @tparam C is the type of numeric constants.
 */
template<typename C = double>
class Expression
{
protected:
  _base_expression_type<C> *_ex; //!< A pointer to a base expression object.

  /**
   * @brief Build a new Expression object.
   *
   * @param _ex is the expression representation.
   */
  Expression(_base_expression_type<C> *_ex);

public:
  typedef std::map<Symbol<C>, Expression<C>> replacement_type;

  /**
   * @brief Build an empty Expression object.
   */
  Expression();

  /**
   * @brief Build a constant expression.
   *
   * @param value is the value of the expression.
   */
  Expression(const int value);

  /**
   * @brief Build a constant expression.
   *
   * @param value is the value of the expression.
   */
  template<typename T = C>
  Expression(const double value,
             typename std::enable_if<!std::is_same<T, double>::value,
                                     int>::type * = nullptr):
      _ex(new _constant_type<C>(static_cast<C>(value)))
  {
  }

  /**
   * @brief Build a constant expression.
   *
   * @param value is the value of the expression.
   */
  Expression(const C value);

  /**
   * @brief Copy constructor.
   *
   * @param orig is the model for the new object.
   */
  Expression(const Expression<C> &orig);

  /**
   * @brief Replace symbol occurences by using expressions.
   *
   * This method replaces any occurence of the symbols in the domain
   * of `replacement` with the corresponding expressions. This method
   * updates this object.
   *
   * @param replacements associates symbols to their replacements.
   * @return A reference to the updated object.
   */
  Expression<C> &replace(const replacement_type &replacements);

  /**
   * @brief Turn the expression into a sum of products.
   *
   * This method turns this expression into an algebraic equivalent
   * sum of products.
   *
   * @return A reference to to the updated object.
   */
  Expression<C> &expand();

  /**
   * @brief Numerically evaluate the expression.
   *
   * @return The numeric evaluation of the expression.
   */
  template<typename T,
           typename
           = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
  inline T evaluate() const
  {
    return static_cast<T>(_ex->evaluate());
  }

  /**
   * @brief Get the coefficient of a term having a given degree.
   *
   * @param symb is the symbol whose term coefficient is aimed.
   * @param degree is the degree of the aimed term.
   * @return The coefficient of the term `symb^{degree}`.
   */
  inline Expression<C> get_coeff(const Symbol<C> &symb, const int degree) const
  {
    return Expression<C>(_ex->get_coeff(symb.get_id(), degree));
  }

  /**
   * @brief Get the coefficients of the polynomial expression.
   *
   * @param symbol_id is the symbol id whose term coefficient is aimed.
   * @return The coefficients of the polynomial expression on the
   *          symbol whose id is `symbol_id`.
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
   * @brief Get the degree of a symbol.
   *
   * @param symb is the symbol whose degree is aimed.
   * @return the degree of the parameter in the expression.
   */
  inline int degree(const Symbol<C> &symb) const
  {
    return _ex->degree(((const _symbol_type<C> *)(symb._ex))->get_id());
  }

  /**
   * @brief Destroy the Expression object.
   */
  ~Expression();

  /**
   * @brief Assign an expression.
   *
   * @param rhs is the expression to be assigned.
   * @return A reference to to the updated object.
   */
  const Expression<C> &operator=(const Expression<C> &rhs);

  /**
   * @brief Assign an expression.
   *
   * @param rhs is the expression to be assigned.
   * @return A reference to to the updated object.
   */
  const Expression<C> &operator=(Expression<C> &&rhs);

  /**
   * @brief Add an expression to the current one.
   *
   * @param rhs is the expression to be added.
   * @return A reference to to the updated object.
   */
  const Expression<C> &operator+=(const Expression<C> &rhs);

  /**
   * @brief Subtract an expression to the current one.
   *
   * @param rhs is the expression to be subtracted.
   * @return A reference to to the updated object.
   */
  const Expression<C> &operator-=(const Expression<C> &rhs);

  /**
   * @brief Multiply an expression to the current one.
   *
   * @param rhs is the expression to be multiplied.
   * @return A reference to to the updated object.
   */
  const Expression<C> &operator*=(const Expression<C> &rhs);

  /**
   * @brief Divide an expression to the current one.
   *
   * @param rhs is the expression to be divided.
   * @return A reference to to the updated object.
   */
  const Expression<C> &operator/=(const Expression<C> &rhs);

  /**
   * @brief Add an expression to the current one.
   *
   * @param rhs is the expression to be added.
   * @return A reference to to the updated object.
   */
  const Expression<C> &operator+=(Expression<C> &&rhs);

  /**
   * @brief Subtract an expression to the current one.
   *
   * @param rhs is the expression to be subtracted.
   * @return A reference to to the updated object.
   */
  const Expression<C> &operator-=(Expression<C> &&rhs);

  /**
   * @brief Multiply an expression to the current one.
   *
   * @param rhs is the expression to be multiplied.
   * @return A reference to to the updated object.
   */
  const Expression<C> &operator*=(Expression<C> &&rhs);

  /**
   * @brief Divide an expression to the current one.
   *
   * @param rhs is the expression to be divided.
   * @return A reference to to the updated object.
   */
  const Expression<C> &operator/=(Expression<C> &&rhs);

  /**
   * @brief Complement the expression.
   *
   * @return the complementar expression.
   */
  Expression<C> operator-() const
  {
    return Expression<C>((_ex->clone())->complement());
  }

  template<typename T>
  friend bool operator==(const Expression<T> &lhs, const T rhs);
  template<typename T>
  friend bool operator>(const Expression<T> &lhs, const T rhs);
  template<typename T>
  friend bool operator>(const T lhs, const Expression<T> &rhs);

  template<typename T>
  friend Expression<T> operator+(const Expression<T> &lhs,
                                 const Expression<T> &rhs);
  template<typename T>
  friend Expression<T> operator-(const Expression<T> &lhs,
                                 const Expression<T> &rhs);
  template<typename T>
  friend Expression<T> operator*(const Expression<T> &lhs,
                                 const Expression<T> &rhs);
  template<typename T>
  friend Expression<T> operator/(const Expression<T> &lhs,
                                 const Expression<T> &rhs);

  template<typename T>
  friend Expression<T> operator+(Expression<T> &&lhs,
                                 const Expression<T> &rhs);
  template<typename T>
  friend Expression<T> operator-(Expression<T> &&lhs,
                                 const Expression<T> &rhs);
  template<typename T>
  friend Expression<T> operator*(Expression<T> &&lhs,
                                 const Expression<T> &rhs);
  template<typename T>
  friend Expression<T> operator/(Expression<T> &&lhs,
                                 const Expression<T> &rhs);

  template<typename T>
  friend Expression<T> operator+(const Expression<T> &lhs,
                                 Expression<T> &&rhs);
  template<typename T>
  friend Expression<T> operator-(const Expression<T> &lhs,
                                 Expression<T> &&rhs);
  template<typename T>
  friend Expression<T> operator*(const Expression<T> &lhs,
                                 Expression<T> &&rhs);
  template<typename T>
  friend Expression<T> operator/(const Expression<T> &lhs,
                                 Expression<T> &&rhs);

  template<typename T>
  friend Expression<T> operator+(Expression<T> &&lhs, Expression<T> &&rhs);
  template<typename T>
  friend Expression<T> operator-(Expression<T> &&lhs, Expression<T> &&rhs);
  template<typename T>
  friend Expression<T> operator*(Expression<T> &&lhs, Expression<T> &&rhs);
  template<typename T>
  friend Expression<T> operator/(Expression<T> &&lhs, Expression<T> &&rhs);

  template<typename T>
  friend void std::swap(Expression<T> &a, Expression<T> &b);

  template<typename T>
  friend std::ostream &std::operator<<(std::ostream &os,
                                       const Expression<T> &ex);
};

/**
 * @brief A base class to represent algebraic expression.
 *
 * @tparam C is the type of numeric constants.
 */
template<typename C>
class _base_expression_type
{
public:
  typedef typename Symbol<C>::SymbolIdType SymbolIdType;

  /**
   * @brief Build an empty expression.
   */
  _base_expression_type() {}

  /**
   * @brief Get the expression type.
   *
   * @return The type of the expression.
   */
  virtual ExpressionType type() const = 0;

  /**
   * @brief Add an expression to the current one.
   *
   * This method builds an expression that represents the sum among the
   * current object and the parameter. This is done by merging the two
   * representations and deallocating the non-necessary components. The
   * original expressions depicted by the current object and by the
   * parameter will not be available anymore after the execution.
   *
   * @param op is the expression to be added.
   * @return a pointer to an expression that represents the sum among
   *              the current object and the parameter.
   */
  virtual _base_expression_type<C> *add(_base_expression_type<C> *op)
  {
    switch (op->type()) {
    case CONSTANT: {
      _finite_sum_type<C> *sum = new _finite_sum_type<C>();

      sum->_constant = ((_constant_type<C> *)op)->get_value();

      delete op;

      return sum->add(this);
    }
    case SYMBOL:
    case FINITE_PROD: {
      _finite_sum_type<C> *sum = new _finite_sum_type<C>();

      sum->_sum.push_back(op);
      return sum->add(this);
    }
    case FINITE_SUM:
    default:
      return op->add(this);
    }
  }

  /**
   * @brief Subtract an expression from the current one.
   *
   * This method builds an expression that represents the subtraction
   * of the parameter from the current object. This is done by merging
   * the two representations and deallocating the non-necessary components.
   * The original expressions depicted by the current object and by the
   * parameter will not be available anymore after the execution.
   *
   * @param op is the expression to be subtracted.
   * @return a pointer to an expression that represents the subtraction
   *          of the parameter from the current object.
   */
  virtual _base_expression_type<C> *subtract(_base_expression_type<C> *op)
  {
    if (op->type() == CONSTANT
        && ((_constant_type<C> *)op)->get_value() == 0) {
      delete op;

      return this;
    }
    return this->add(op->complement());
  }

  /**
   * @brief Multiply an expression to the current one.
   *
   * This method builds an expression that represents the multiplication
   * of the parameter and the current object. This is done by merging
   * the two representations and deallocating the non-necessary components.
   * The original expressions depicted by the current object and by the
   * parameter will not be available anymore after the execution.
   *
   * @param op is the expression that multiplies the current object.
   * @return a pointer to an expression that represents the multiplication
   *          of the parameter and the current object.
   */
  virtual _base_expression_type<C> *multiply(_base_expression_type<C> *op)
  {
    switch (op->type()) {
    case CONSTANT: {
      _finite_prod_type<C> *prod = new _finite_prod_type<C>();

      prod->_constant = ((_constant_type<C> *)op)->get_value();

      delete op;

      return prod->multiply(this);
    }
    case SYMBOL:
    case FINITE_SUM: {
      _finite_prod_type<C> *prod = new _finite_prod_type<C>();

      prod->_numerator.push_back(op);
      return prod->multiply(this);
    }
    case FINITE_PROD:
    default:
      return op->multiply(this);
    }
  }

  /**
   * @brief Divide the current expression by another one.
   *
   * This method builds an expression that represents the division
   * of the current object by the parameter. This is done by merging
   * the two representations and deallocating the non-necessary components.
   * The original expressions depicted by the current object and by the
   * parameter will not be available anymore after the execution.
   *
   * @param op is the expression that divides the current object.
   * @return a pointer to an expression that represents the division
   *          of the current object by the parameter.
   */
  virtual _base_expression_type<C> *divide(_base_expression_type<C> *op)
  {
    switch (op->type()) {
    case CONSTANT: {
      _finite_prod_type<C> *prod = new _finite_prod_type<C>();

      prod->_constant = 1 / ((_constant_type<C> *)op)->get_value();

      delete op;

      return prod->multiply(this);
    }
    case SYMBOL:
    case FINITE_SUM: {
      _finite_prod_type<C> *prod = new _finite_prod_type<C>();

      prod->_denominator.push_back(op);
      return prod->multiply(this);
    }
    case FINITE_PROD:
    default:
      _finite_prod_type<C> *prod = new _finite_prod_type<C>();
      _finite_prod_type<C> *op_prod = (_finite_prod_type<C> *)op;

      std::swap(prod->_denominator, op_prod->_numerator);
      std::swap(prod->_numerator, op_prod->_denominator);
      prod->_constant = 1 / op_prod->_constant;
      return prod->multiply(this);
    }
  }

  /**
   * @brief Complement the expression.
   *
   * @return the complementar expression.
   */
  virtual _base_expression_type<C> *complement() = 0;

  /**
   * @brief Multiply all the products in a list by the current object.
   *
   * This method multiplies all the product in a list by the current object and
   * returns the list of the products. The list of products will be not
   * available anymore after the execution.
   *
   * @param prods is a list of products.
   * @return the list of the products between the current object and the
   * products in the `prods`.
   */
  virtual std::list<_base_expression_type<C> *>
  multiply(std::list<_base_expression_type<C> *> &prods)
  {
    std::list<_base_expression_type<C> *> result;

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
   * @return A sum of products, a constant, or a symbol.
   */
  virtual _base_expression_type<C> *expand() const
  {
    return this->clone();
  }

  /**
   * @brief Clone the current object.
   *
   * @return A clone of the current object.
   */
  virtual _base_expression_type<C> *clone() const = 0;

  /**
   * @brief Get the ids of the symbols in the expression.
   *
   * @return The set of the symbol ids in the expression.
   */
  virtual std::set<SymbolIdType> get_symbol_ids() const = 0;

  /**
   * @brief Replace symbol occurences by using expressions.
   *
   * This method replaces any occurence of the symbols whose ids are in
   * the domain of `replacement` with a clone of the corresponding expressions.
   * This method updates this object.
   *
   * @param replacements associates symbol ids to their replacements.
   * @return A pointer to the updated object.
   */
  virtual _base_expression_type<C> *replace(
      const std::map<SymbolIdType, _base_expression_type<C> *> &replacements)
      = 0;

  /**
   * @brief Numerically evaluate the expression.
   *
   * @return The numeric evaluation of the expression.
   */
  virtual C evaluate() const = 0;

  /**
   * @brief Get the coefficient of a term having a given degree.
   *
   * @param symbol_id is the symbol id whose term coefficient is aimed.
   * @param degree is the degree of the aimed term.
   * @return The coefficient of the term in which the symbol having id
   * `symbol_id` has degree `degree`.
   */
  virtual _base_expression_type<C> *get_coeff(const SymbolIdType &symbol_id,
                                              const int &degree) const = 0;

  /**
   * @brief Get the coefficients of the polynomial expression.
   *
   * @param symbol_id is the symbol id whose term coefficient is aimed.
   * @return The coefficients of the polynomial expression on the
   *          symbol whose id is `symbol_id`.
   */
  virtual std::map<int, _base_expression_type<C> *>
  get_coeffs(const SymbolIdType &symbol_id) const = 0;

  /**
   * @brief Get the degree of a symbol.
   *
   * @param symbol_id is the id of the symbol whose degree is aimed.
   * @return the degree of the parameter in the expression.
   */
  virtual int degree(const SymbolIdType &symbol_id) const = 0;

  /**
   * @brief Print the expression in an output stream.
   *
   * @param os is the output stream in which the expression must be printed.
   */
  virtual void print(std::ostream &os) const = 0;

  /**
   * @brief Check whether this is the constant 0.
   *
   * @return true if and only if this is the constant 0.
   */
  bool is_zero() const
  {
    return (type() == CONSTANT
            && ((_constant_type<C> *)this)->get_value() == 0);
  }

  /**
   * @brief Destroy the base expression type object
   */
  virtual ~_base_expression_type() {}
};

/**
 * @brief A class to represent constants.
 *
 * @tparam C is the type of numeric constants.
 */
template<typename C>
class _constant_type : public _base_expression_type<C>
{
  C _value; //!< The numerical value of the constant.
public:
  typedef typename Symbol<C>::SymbolIdType SymbolIdType;

  /**
   * @brief Build a new constant expression.
   *
   * @param value is the value of the constant.
   */
  _constant_type(const C value): _base_expression_type<C>(), _value(value) {}

  /**
   * @brief Copy constructor.
   *
   * @param orig is the model for the new constant.
   */
  _constant_type(const _constant_type<C> &orig):
      _base_expression_type<C>(), _value(orig._value)
  {
  }

  /**
   * @brief Get the constant value.
   *
   * @return the value of the constant.
   */
  const C &get_value() const
  {
    return _value;
  }

  /**
   * @brief Get the expression type.
   *
   * @return The value CONSTANT.
   */
  ExpressionType type() const
  {
    return CONSTANT;
  }

  /**
   * @brief Add an expression to the current one.
   *
   * This method builds an expression that represents the sum among the
   * current object and the parameter. This is done by merging the two
   * representations and deallocating the non-necessary components. The
   * original expressions depicted by the current object and by the
   * parameter will not be available anymore after the execution.
   *
   * @param op is the expression to be added.
   * @return a pointer to an expression that represents the sum among
   *              the current object and the parameter.
   */
  _base_expression_type<C> *add(_base_expression_type<C> *op)
  {
    switch (op->type()) {
    case CONSTANT:
      _value += ((_constant_type<C> *)op)->_value;

      delete op;

      return this;

    default:
      return op->add(this);
    }
  }

  /**
   * @brief Subtract an expression from the current one.
   *
   * This method builds an expression that represents the subtraction
   * of the parameter from the current object. This is done by merging
   * the two representations and deallocating the non-necessary components.
   * The original expressions depicted by the current object and by the
   * parameter will not be available anymore after the execution.
   *
   * @param op is the expression to be subtracted.
   * @return a pointer to an expression that represents the subtraction
   *          of the parameter from the current object.
   */
  _base_expression_type<C> *subtract(_base_expression_type<C> *op)
  {
    switch (op->type()) {
    case CONSTANT:
      _value -= ((_constant_type<C> *)op)->_value;

      delete op;

      return this;

    default:
      _constant_type<C> *minus1 = new _constant_type<C>(-1);
      return this->add(minus1->multiply(op));
    }
  }

  /**
   * @brief Multiply an expression to the current one.
   *
   * This method builds an expression that represents the multiplication
   * of the parameter and the current object. This is done by merging
   * the two representations and deallocating the non-necessary components.
   * The original expressions depicted by the current object and by the
   * parameter will not be available anymore after the execution.
   *
   * @param op is the expression that multiplies the current object.
   * @return a pointer to an expression that represents the multiplication
   *          of the parameter and the current object.
   */
  _base_expression_type<C> *multiply(_base_expression_type<C> *op)
  {
    switch (op->type()) {
    case CONSTANT:
      _value *= ((_constant_type<C> *)op)->_value;

      delete op;

      return this;

    default:
      return op->multiply(this);
    }
  }

  /**
   * @brief Divide the current expression by another one.
   *
   * This method builds an expression that represents the division
   * of the current object by the parameter. This is done by merging
   * the two representations and deallocating the non-necessary components.
   * The original expressions depicted by the current object and by the
   * parameter will not be available anymore after the execution.
   *
   * @param op is the expression that divides the current object.
   * @return a pointer to an expression that represents the division
   *          of the current object by the parameter.
   */
  _base_expression_type<C> *divide(_base_expression_type<C> *op)
  {
    switch (op->type()) {
    case CONSTANT:
      _value /= ((_constant_type<C> *)op)->_value;

      delete op;

      return this;

    default:
      _finite_prod_type<C> *prod = new _finite_prod_type<C>();

      prod->multiply(this);
      return prod->divide(op);
    }
  }

  /**
   * @brief Complement the expression.
   *
   * @return the complementar expression.
   */
  _base_expression_type<C> *complement()
  {
    _value = -_value;

    return this;
  }

  /**
   * @brief Replace symbol occurences by using expressions.
   *
   * This method replaces any occurence of the symbols whose ids are in
   * the domain of `replacement` with a clone of the corresponding expressions.
   * This method updates this object.
   *
   * @param replacements associates symbol ids to their replacements.
   * @return A pointer to the updated object.
   */
  _base_expression_type<C> *replace(
      const std::map<SymbolIdType, _base_expression_type<C> *> &replacements)
  {
    (void)replacements;

    return this;
  }

  /**
   * @brief Clone the current object.
   *
   * @return A clone of the current object.
   */
  _base_expression_type<C> *clone() const
  {
    return new _constant_type<C>(*this);
  }

  /**
   * @brief Get the ids of the symbols in the expression.
   *
   * @return The set of the symbol ids in the expression.
   */
  std::set<SymbolIdType> get_symbol_ids() const
  {
    return std::set<SymbolIdType>();
  }

  /**
   * @brief Numerically evaluate the expression.
   *
   * @return The numeric evaluation of the expression.
   */
  C evaluate() const
  {
    return _value;
  }

  /**
   * @brief Get the coefficient of a term having a given degree.
   *
   * @param symbol_id is the symbol id whose term coefficient is aimed.
   * @param degree is the degree of the aimed term.
   * @return The coefficient of the term in which the symbol having id
   * `symbol_id` has degree `degree`.
   */
  _base_expression_type<C> *get_coeff(const SymbolIdType &symbol_id,
                                      const int &degree) const
  {
    (void)symbol_id;

    return new _constant_type(degree == 0 ? _value : 0);
  }

  /**
   * @brief Get the coefficients of the polynomial expression.
   *
   * @param symbol_id is the symbol id whose term coefficient is aimed.
   * @return The coefficients of the polynomial expression on the
   *          symbol whose id is `symbol_id`.
   */
  std::map<int, _base_expression_type<C> *>
  get_coeffs(const SymbolIdType &symbol_id) const
  {
    (void)symbol_id;

    std::map<int, _base_expression_type<C> *> res;

    res[0] = new _constant_type<C>(_value);

    return res;
  }

  /**
   * @brief Get the degree of a symbol.
   *
   * @param symbol_id is the id of the symbol whose degree is aimed.
   * @return the degree of the parameter in the expression, i.e., 0.
   */
  int degree(const SymbolIdType &symbol_id) const
  {
    (void)symbol_id;

    return 0;
  }

  /**
   * @brief Print the expression in an output stream.
   *
   * @param os is the output stream in which the expression must be printed.
   */
  void print(std::ostream &os) const
  {
    os << _value;
  }
};

/**
 * @brief A class to represent finite sums.
 *
 * @tparam C is the type of numeric constants.
 */
template<typename C>
class _finite_sum_type : public _base_expression_type<C>
{
  C _constant;
  std::list<_base_expression_type<C> *>
      _sum; //!< The finite list of non-constant expressions.
public:
  typedef typename Symbol<C>::SymbolIdType SymbolIdType;

  /**
   * @brief Build an empty finite sum.
   */
  _finite_sum_type(): _base_expression_type<C>(), _constant(0), _sum() {}

  /**
   * @brief Copy constructor.
   *
   * @param orig is the model for the new finite sum.
   */
  _finite_sum_type(const _finite_sum_type<C> &orig):
      _base_expression_type<C>(), _constant(orig._constant), _sum()
  {
    for (auto it = std::begin(orig._sum); it != std::end(orig._sum); ++it) {
      _sum.push_back((*it)->clone());
    }
  }

  /**
   * @brief Build a finite sum of products.
   *
   * @param orig is a list of products that must be added.
   */
  _finite_sum_type(std::list<_base_expression_type<C> *> &model):
      _base_expression_type<C>(), _constant(0), _sum()
  {
    for (auto it = std::begin(model); it != std::end(model); ++it) {
      _sum.push_back(*it);
    }

    model.clear();
  }

  /**
   * @brief Get the expression type.
   *
   * @return The value FINITE_SUM.
   */
  ExpressionType type() const
  {
    return FINITE_SUM;
  }

  /**
   * @brief Add an expression to the current one.
   *
   * This method builds an expression that represents the sum among the
   * current object and the parameter. This is done by merging the two
   * representations and deallocating the non-necessary components. The
   * original expressions depicted by the current object and by the
   * parameter will not be available anymore after the execution.
   *
   * @param op is the expression to be added.
   * @return a pointer to an expression that represents the sum among
   *              the current object and the parameter.
   */
  _base_expression_type<C> *add(_base_expression_type<C> *op)
  {
    switch (op->type()) {
    case FINITE_SUM: {
      _finite_sum_type<C> *op_sum = (_finite_sum_type<C> *)op;

      _constant += op_sum->_constant;
      _sum.splice(std::end(_sum), op_sum->_sum);

      delete op;

      return this;
    }
    case CONSTANT:
      _constant += ((_constant_type<C> *)op)->get_value();

      delete op;

      if (_sum.size() == 0) {
        _constant_type<C> *res = new _constant_type<C>(_constant);

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
   * @brief Multiply all the products in a list by the current object.
   *
   * This method multiplies all the product in a list by the current object and
   * returns the list of the products. The list of products will be not
   * available anymore after the execution.
   *
   * @param prods is a list of products.
   * @return the list of the products between the current object and the
   * products in the `prods`.
   */
  std::list<_base_expression_type<C> *>
  multiply(std::list<_base_expression_type<C> *> &prods)
  {
    std::list<_base_expression_type<C> *> result;

    for (auto s_it = std::begin(_sum); s_it != std::end(_sum); ++s_it) {
      for (auto p_it = std::begin(prods); p_it != std::end(prods); ++p_it) {
        _base_expression_type<C> *p_clone = (*p_it)->clone();

        result.push_back(p_clone->multiply((*s_it)->clone()));
      }
      delete *s_it;
    }
    _sum.clear();

    for (auto p_it = std::begin(prods); p_it != std::end(prods); ++p_it) {
      result.push_back((*p_it)->multiply(new _constant_type<C>(_constant)));
    }
    prods.clear();

    return result;
  }

  /**
   * @brief Complement the expression.
   *
   * @return the complementar expression.
   */
  _base_expression_type<C> *complement()
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
   * @return A sum of products, a constant, or a symbol.
   */
  _base_expression_type<C> *expand() const
  {
    _finite_sum_type<C> *new_obj = new _finite_sum_type<C>();
    std::list<_base_expression_type<C> *> &new_sum = new_obj->_sum;
    C &constant = new_obj->_constant;

    constant = _constant;

    for (auto it = std::begin(_sum); it != std::end(_sum); ++it) {
      _base_expression_type<C> *ex_it = (*it)->expand();

      switch (ex_it->type()) {
      case FINITE_SUM: {
        _finite_sum_type<C> *ex_it_sum = (_finite_sum_type<C> *)ex_it;
        for (auto e_it = std::begin(ex_it_sum->_sum);
             e_it != std::end(ex_it_sum->_sum); ++e_it) {
          if ((*e_it)->type() == CONSTANT) {
            constant += ((_constant_type<C> *)(*e_it))->get_value();

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
        constant += ((_constant_type<C> *)ex_it)->get_value();

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
   * @brief Clone the current object.
   *
   * @return A clone of the current object.
   */
  _base_expression_type<C> *clone() const
  {
    return new _finite_sum_type<C>(*this);
  }

  /**
   * @brief Get the ids of the symbols in the expression.
   *
   * @return The set of the symbol ids in the expression.
   */
  std::set<SymbolIdType> get_symbol_ids() const
  {
    std::set<SymbolIdType> ids;

    for (auto it = std::begin(_sum); it != std::end(_sum); ++it) {
      std::set<SymbolIdType> it_ids = (*it)->get_symbol_ids();

      for (auto it_it = std::begin(it_ids); it_it != std::end(it_ids);
           ++it_it) {
        ids.insert(*it_it);
      }
    }

    return ids;
  }

  /**
   * @brief Replace symbol occurences by using expressions.
   *
   * This method replaces any occurence of the symbols whose ids are in
   * the domain of `replacement` with a clone of the corresponding expressions.
   * This method updates this object.
   *
   * @param replacements associates symbol ids to their replacements.
   * @return A pointer to the updated object.
   */
  _base_expression_type<C> *replace(
      const std::map<SymbolIdType, _base_expression_type<C> *> &replacements)
  {
    std::list<_base_expression_type<C> *> new_sum;
    for (auto it = std::begin(_sum); it != std::end(_sum); ++it) {
      _base_expression_type<C> *r_it = (*it)->replace(replacements);
      switch (r_it->type()) {
      case FINITE_SUM: {
        _finite_sum_type<C> *r_it_sum = (_finite_sum_type<C> *)r_it;
        _constant += r_it_sum->_constant;

        new_sum.splice(std::end(new_sum), r_it_sum->_sum);

        delete r_it;
        break;
      }
      case CONSTANT:
        _constant += ((_constant_type<C> *)r_it)->get_value();

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

    _constant_type<C> *result = new _constant_type<C>(_constant);
    delete this;

    return result;
  }

  /**
   * @brief Numerically evaluate the expression.
   *
   * @return The numeric evaluation of the expression.
   */
  C evaluate() const
  {
    C total = 0;
    for (auto it = std::begin(_sum); it != std::end(_sum); ++it) {
      total += (*it)->evaluate();
    }

    return total;
  }

  /**
   * @brief Get the coefficient of a term having a given degree.
   *
   * @param symbol_id is the symbol id whose term coefficient is aimed.
   * @param degree is the degree of the aimed term.
   * @return The coefficient of the term in which the symbol having id
   * `symbol_id` has degree `degree`.
   */
  _base_expression_type<C> *get_coeff(const SymbolIdType &symbol_id,
                                      const int &degree) const
  {
    _base_expression_type<C> *total_coeff = new _finite_sum_type<C>();

    if (degree == 0) {
      ((_finite_sum_type<C> *)total_coeff)->_constant = this->_constant;
    }

    for (auto it = std::begin(_sum); it != std::end(_sum); ++it) {
      total_coeff = total_coeff->add((*it)->get_coeff(symbol_id, degree));
    }

    return total_coeff;
  }

  /**
   * @brief Get the coefficients of the polynomial expression.
   *
   * @param symbol_id is the symbol id whose term coefficient is aimed.
   * @return The coefficients of the polynomial expression on the
   *          symbol whose id is `symbol_id`.
   */
  std::map<int, _base_expression_type<C> *>
  get_coeffs(const SymbolIdType &symbol_id) const
  {
    std::map<int, _base_expression_type<C> *> res;

    res[0] = new _constant_type<C>(this->_constant);

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
   * @brief Get the degree of a symbol.
   *
   * @param symbol_id is the id of the symbol whose degree is aimed.
   * @return the degree of the parameter in the expression, i.e.,
   *          the maximum degree among the addends.
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
   * @brief Print the expression in an output stream.
   *
   * @param os is the output stream in which the expression must be printed.
   */
  void print(std::ostream &os) const
  {
    if (_constant != 0 || _sum.size() == 0) {
      os << _constant;
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

  ~_finite_sum_type()
  {
    for (auto it = std::begin(_sum); it != std::end(_sum); ++it) {
      delete *it;
    }
  }

  template<typename T>
  friend class _base_expression_type;
};

/**
 * @brief A class to represent finite products.
 *
 * @tparam C is the type of numeric constants.
 */
template<typename C>
class _finite_prod_type : public _base_expression_type<C>
{
  C _constant; //!< The numeric constant that multiplies this product.
  std::list<_base_expression_type<C> *>
      _numerator; //!< The numerator list of this product.
  std::list<_base_expression_type<C> *>
      _denominator; //!< The denominator list of this product.

  static void print_list(std::ostream &os,
                         const std::list<_base_expression_type<C> *> &_list,
                         const C constant = 1)
  {
    if (constant != 1 || _list.size() == 0) {
      os << constant;
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
   * @brief Constructor.
   *
   * @param constant is a constant to initialize the new product.
   */
  _finite_prod_type(const C &constant):
      _base_expression_type<C>(), _constant(constant), _numerator(),
      _denominator()
  {
  }

public:
  typedef typename Symbol<C>::SymbolIdType SymbolIdType;

  /**
   * @brief Build an empty finite product.
   */
  _finite_prod_type():
      _base_expression_type<C>(), _constant(1), _numerator(), _denominator()
  {
  }

  /**
   * @brief Copy constructor.
   *
   * @param orig is a model for the new finite product.
   */
  _finite_prod_type(const _finite_prod_type<C> &orig):
      _base_expression_type<C>(), _constant(orig._constant), _numerator(),
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
   * @brief Get the expression type.
   *
   * @return The value FINITE_PROD.
   */
  ExpressionType type() const
  {
    return FINITE_PROD;
  }

  /**
   * @brief Multiply an expression to the current one.
   *
   * This method builds an expression that represents the multiplication
   * of the parameter and the current object. This is done by merging
   * the two representations and deallocating the non-necessary components.
   * The original expressions depicted by the current object and by the
   * parameter will not be available anymore after the execution.
   *
   * @param op is the expression that multiplies the current object.
   * @return a pointer to an expression that represents the multiplication
   *          of the parameter and the current object.
   */
  _base_expression_type<C> *multiply(_base_expression_type<C> *op)
  {
    switch (op->type()) {
    case FINITE_PROD: {
      _finite_prod_type<C> *op_prod = (_finite_prod_type<C> *)op;

      _constant *= op_prod->_constant;

      _numerator.splice(std::end(_numerator), op_prod->_numerator);
      _denominator.splice(std::end(_denominator), op_prod->_denominator);

      delete op_prod;

      break;
    }
    case CONSTANT: {
      _constant_type<C> *op_const = (_constant_type<C> *)op;

      _constant *= op_const->get_value();

      delete op_const;

      break;
    }
    default:
      _numerator.push_back(op);
    }

    if (_denominator.size() == 0) {
      if (_constant == 0 || _numerator.size() == 0) {
        _constant_type<C> *result = new _constant_type<C>(_constant);

        delete this;

        return result;
      }
      if (_constant == 1 && _numerator.size() == 1) {

        _base_expression_type<C> *result = _numerator.front();

        _numerator.clear();

        delete this;

        return result;
      }
    }

    return this;
  }

  /**
   * @brief Divide the current expression by another one.
   *
   * This method builds an expression that represents the division
   * of the current object by the parameter. This is done by merging
   * the two representations and deallocating the non-necessary components.
   * The original expressions depicted by the current object and by the
   * parameter will not be available anymore after the execution.
   *
   * @param op is the expression that divides the current object.
   * @return a pointer to an expression that represents the division
   *          of the current object by the parameter.
   */
  _base_expression_type<C> *divide(_base_expression_type<C> *op)
  {
    switch (op->type()) {
    case FINITE_PROD: {
      _finite_prod_type<C> *op_prod = (_finite_prod_type<C> *)op;

      _constant /= op_prod->_constant;

      _numerator.splice(std::end(_numerator), op_prod->_denominator);
      _denominator.splice(std::end(_denominator), op_prod->_numerator);

      delete op;

      break;
    }
    case CONSTANT: {
      _constant_type<C> *op_const = (_constant_type<C> *)op;

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
   * @brief Complement the expression.
   *
   * @return the complementar expression.
   */
  _base_expression_type<C> *complement()
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
   * @return A sum of products, a constant, or a symbol.
   */
  _base_expression_type<C> *expand() const
  {
    C base_const = this->_constant;

    try {
      for (auto it = std::begin(_denominator); it != std::begin(_denominator);
           ++it) {
        base_const /= (*it)->evaluate();
      }
    } catch (std::runtime_error &) {
      throw std::runtime_error("Non-constant denominator not supported.");
    }

    std::list<_base_expression_type<C> *> result;

    result.push_back(new _constant_type<C>(base_const));

    for (auto it = std::begin(_numerator); it != std::end(_numerator); ++it) {
      _base_expression_type<C> *new_it = (*it)->expand();

      result = new_it->multiply(result);

      delete new_it;
    }
    return new _finite_sum_type<C>(result);
  }

  /**
   * @brief Clone the current object.
   *
   * @return A clone of the current object.
   */
  _base_expression_type<C> *clone() const
  {
    return new _finite_prod_type<C>(*this);
  }

  /**
   * @brief Get the ids of the symbols in the expression.
   *
   * @return The set of the symbol ids in the expression.
   */
  std::set<SymbolIdType> get_symbol_ids() const
  {
    std::set<SymbolIdType> ids;

    for (auto &_list: {_numerator, _denominator}) {
      for (auto it = std::begin(_list); it != std::end(_list); ++it) {
        std::set<SymbolIdType> it_ids = (*it)->get_symbol_ids();

        for (auto it_it = std::begin(it_ids); it_it != std::end(it_ids);
             ++it_it) {
          ids.insert(*it_it);
        }
      }
    }

    return ids;
  }

  /**
   * @brief Replace symbol occurences by using expressions.
   *
   * This method replaces any occurence of the symbols whose ids are in
   * the domain of `replacement` with a clone of the corresponding expressions.
   * This method updates this object.
   *
   * @param replacements associates symbol ids to their replacements.
   * @return A pointer to the updated object.
   */
  _base_expression_type<C> *replace(
      const std::map<SymbolIdType, _base_expression_type<C> *> &replacements)
  {
    _base_expression_type<C> *result = new _constant_type<C>(this->_constant);
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
        result->divide((*it)->replace(replacements));
      } else {
        delete *it;
      }
    }
    _denominator.clear();

    return result;
  }

  /**
   * @brief Numerically evaluate the expression.
   *
   * @return The numeric evaluation of the expression.
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
   * @brief Get the coefficient of a term having a given degree.
   *
   * @param symbol_id is the symbol id whose term coefficient is aimed.
   * @param degree is the degree of the aimed term.
   * @return The coefficient of the term in which the symbol having id
   * `symbol_id` has degree `degree`.
   */
  _base_expression_type<C> *get_coeff(const SymbolIdType &symbol_id,
                                      const int &degree) const
  {
    C constant = this->_constant;

    try {
      for (auto it = std::begin(_denominator); it != std::begin(_denominator);
           ++it) {
        constant /= (*it)->evaluate();
      }
    } catch (std::runtime_error &) {
      throw std::runtime_error("Non-constant denominator not supported.");
    }

    _base_expression_type<C> *res = new _constant_type<C>(constant);

    int degree_counter = 0;
    for (auto it = std::begin(_numerator); it != std::end(_numerator); ++it) {
      switch ((*it)->type()) {
      case SYMBOL:
        if (((_symbol_type<C> *)(*it))->get_id() == symbol_id) {
          ++degree_counter;
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
        throw std::runtime_error(
            "'coeff' call only admitted after 'expand' call.");
      }
    }

    if (degree_counter == degree) {
      return res;
    }

    delete res;

    return new _constant_type<C>(0);
  }

  /**
   * @brief Get the coefficients of the polynomial expression.
   *
   * @param symbol_id is the symbol id whose term coefficient is aimed.
   * @return The coefficients of the polynomial expression on the
   *          symbol whose id is `symbol_id`.
   */
  std::map<int, _base_expression_type<C> *>
  get_coeffs(const SymbolIdType &symbol_id) const
  {
    if (_denominator.size() > 0) {
      throw std::runtime_error(
          "coeff does not support non-constant denominator");
    }

    std::map<int, _base_expression_type<C> *> res;

    res[0] = new _constant_type<C>(this->_constant);

    for (auto it = std::begin(_numerator); it != std::end(_numerator); ++it) {
      auto map_it = (*it)->get_coeffs(symbol_id);
      std::map<int, _base_expression_type<C> *> next_res;
      for (auto res_coeff = std::begin(res); res_coeff != std::end(res);
           ++res_coeff) {
        for (auto it_coeff = std::begin(map_it); it_coeff != std::end(map_it);
             ++it_coeff) {
          const int degree = res_coeff->first + it_coeff->first;
          auto res_coeff_clone = res_coeff->second->clone();
          res_coeff_clone = res_coeff_clone->multiply(it_coeff->second);
          if (next_res.find(degree) == std::end(next_res)) {
            next_res[degree] = res_coeff_clone;
          } else {
            next_res[degree] = next_res[degree]->add(res_coeff_clone);
          }
        }
        delete res_coeff->second;
      }

      std::swap(res, next_res);
    }

    return res;
  }

  /**
   * @brief Get the degree of a symbol.
   *
   * @param symbol_id is the id of the symbol whose degree is aimed.
   * @return the degree of the parameter in the expression, i.e.,
   *          the sum degree of the expression degree in the
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
   * @brief Print the expression in an output stream.
   *
   * @param os is the output stream in which the expression must be printed.
   */
  void print(std::ostream &os) const
  {
    if (_numerator.size() > 0 || _constant != 1) {
      if (_denominator.size() > 0) {
        os << "(";
        print_list(os, _numerator, _constant);
        os << ")";
      } else {
        print_list(os, _numerator, _constant);
      }
    }
    if (_denominator.size() > 0) {
      os << "/";
      if (_denominator.size() > 1) {
        os << "(";
        print_list(os, _denominator);
        os << ")";
      } else {
        print_list(os, _denominator);
      }
    }
  }

  ~_finite_prod_type()
  {
    for (auto &_list: {_numerator, _denominator}) {
      for (auto it = std::begin(_list); it != std::end(_list); ++it) {
        delete *it;
      }
    }
  }

  template<typename T>
  friend class _base_expression_type;
};

/**
 * @brief A class to represent symbols.
 *
 * @tparam C is the type of numeric constants.
 */
template<typename C>
class _symbol_type : public _base_expression_type<C>
{
public:
  typedef typename Symbol<C>::SymbolIdType SymbolIdType;

protected:
  SymbolIdType _id; //!< The symbol id.

public:
  /**
   * @brief Build a new symbol
   *
   * @param id is the id of the new symbol.
   */
  _symbol_type(const SymbolIdType id): _base_expression_type<C>(), _id(id) {}

  /**
   * @brief Copy constructor.
   *
   * @param orig is the model for the new symbol.
   */
  _symbol_type(const _symbol_type<C> &orig):
      _base_expression_type<C>(), _id(orig._id)
  {
  }

  /**
   * @brief Get the id of the symbol.
   *
   * @return The id of the symbol.
   */
  const SymbolIdType &get_id() const
  {
    return _id;
  }

  /**
   * @brief Get the expression type.
   *
   * @return The value SYMBOL.
   */
  ExpressionType type() const
  {
    return SYMBOL;
  }

  /**
   * @brief Complement the expression.
   *
   * @return the complementar expression.
   */
  _base_expression_type<C> *complement()
  {
    return this->multiply(new _constant_type<C>(static_cast<C>(-1)));
  }

  /**
   * @brief Clone the current object.
   *
   * @return A clone of the current object.
   */
  _base_expression_type<C> *clone() const
  {
    return new _symbol_type(*this);
  }

  /**
   * @brief Get the ids of the symbols in the expression.
   *
   * @return The set of the symbol ids in the expression.
   */
  std::set<SymbolIdType> get_symbol_ids() const
  {
    std::set<SymbolIdType> ids;

    ids.insert(_id);

    return ids;
  }

  /**
   * @brief Replace symbol occurences by using expressions.
   *
   * This method replaces any occurence of the symbols whose ids are in
   * the domain of `replacement` with a clone of the corresponding expressions.
   * This method updates this object.
   *
   * @param replacements associates symbol ids to their replacements.
   * @return A pointer to the updated object.
   */
  _base_expression_type<C> *replace(
      const std::map<SymbolIdType, _base_expression_type<C> *> &replacements)
  {
    auto it = replacements.find(_id);

    if (it == std::end(replacements)) {
      return this;
    }

    delete this;

    return it->second->clone();
  }

  /**
   * @brief Numerically evaluate the expression.
   *
   * @return The numeric evaluation of the expression.
   */
  C evaluate() const
  {
    throw std::runtime_error("Symbols cannot be evaluated");
  }

  /**
   * @brief Get the coefficient of a term having a given degree.
   *
   * @param symbol_id is the symbol id whose term coefficient is aimed.
   * @param degree is the degree of the aimed term.
   * @return The coefficient of the term in which the symbol having id
   * `symbol_id` has degree `degree`.
   */
  _base_expression_type<C> *get_coeff(const SymbolIdType &symbol_id,
                                      const int &degree) const
  {
    return new _constant_type<C>((degree == 1 && symbol_id == _id) ? 1 : 0);
  }

  /**
   * @brief Get the coefficients of the polynomial expression.
   *
   * @param symbol_id is the symbol id whose term coefficient is aimed.
   * @return The coefficients of the polynomial expression on the
   *          symbol whose id is `symbol_id`.
   */
  std::map<int, _base_expression_type<C> *>
  get_coeffs(const SymbolIdType &symbol_id) const
  {
    std::map<int, _base_expression_type<C> *> res;

    if (symbol_id == _id) {
      res[1] = new _constant_type<C>(1);
    } else {
      res[0] = this->clone();
    }

    return res;
  }

  /**
   * @brief Get the degree of a symbol.
   *
   * @param symbol_id is the id of the symbol whose degree is aimed.
   * @return the degree of the parameter in the expression, i.e.,
   *          1, if the current symbol is the parameter; 0, otherwise
   */
  int degree(const SymbolIdType &symbol_id) const
  {
    return ((((const _symbol_type<C> *)this)->_id == symbol_id) ? 1 : 0);
  }
  /**
   * @brief Print the expression in an output stream.
   *
   * @param os is the output stream in which the expression must be printed.
   */
  void print(std::ostream &os) const
  {
    os << Symbol<C>::get_symbol_name(_id);
  }
};

template<>
template<>
inline double Expression<mpq_class>::evaluate<double>() const
{
  mpq_class value = _ex->evaluate();

  return value.get_d();
}

template<typename C>
Expression<C>::Expression(_base_expression_type<C> *_ex): _ex(_ex)
{
}

template<typename C>
Expression<C>::Expression(): _ex(NULL)
{
}

template<typename C>
Expression<C>::Expression(const int value):
    _ex(new _constant_type<C>(static_cast<C>(value)))
{
}

/*
template<typename C>
Expression<C>::Expression(const double value):
    _ex(new _constant_type<C>(static_cast<C>(value)))
{
}
*/

template<typename C>
Expression<C>::Expression(const C value): _ex(new _constant_type<C>(value))
{
}

template<typename C>
Expression<C>::Expression(const Expression<C> &orig):
    _ex((orig._ex == NULL ? NULL : orig._ex->clone()))
{
}

template<typename C>
Expression<C> &
Expression<C>::replace(const Expression<C>::replacement_type &replacement)
{
  if (this->_ex == NULL) {
    return *this;
  }
  std::map<typename Symbol<C>::SymbolIdType, _base_expression_type<C> *>
      base_repl;

  for (auto it = std::begin(replacement); it != std::end(replacement); ++it) {
    base_repl[it->first.get_id()] = it->second._ex;
  }

  this->_ex = this->_ex->replace(base_repl);

  return *this;
}

template<typename C>
Expression<C> &Expression<C>::expand()
{
  _base_expression_type<C> *new_ex = _ex->expand();

  delete _ex;

  _ex = new_ex;

  return *this;
}

template<typename C>
Expression<C>::~Expression()
{
  if (_ex != NULL) {
    delete _ex;
    _ex = NULL;
  }
}

template<typename C>
const Expression<C> &Expression<C>::operator=(const Expression<C> &rhs)
{
  delete _ex;

  _ex = rhs._ex->clone();

  return *this;
}

template<typename C>
const Expression<C> &Expression<C>::operator=(Expression<C> &&rhs)
{
  std::swap(_ex, rhs._ex);

  return *this;
}

template<typename C>
const Expression<C> &Expression<C>::operator+=(const Expression<C> &rhs)
{
  _base_expression_type<C> *rhs_ex = rhs._ex->clone();

  _ex = _ex->add(rhs_ex);

  return *this;
}

template<typename C>
const Expression<C> &Expression<C>::operator-=(const Expression<C> &rhs)
{
  _base_expression_type<C> *rhs_ex = rhs._ex->clone();

  _ex = _ex->subtract(rhs_ex);

  return *this;
}

template<typename C>
const Expression<C> &Expression<C>::operator*=(const Expression<C> &rhs)
{
  _base_expression_type<C> *rhs_ex = rhs._ex->clone();

  _ex = _ex->multiply(rhs_ex);

  return *this;
}

template<typename C>
const Expression<C> &Expression<C>::operator/=(const Expression<C> &rhs)
{
  _base_expression_type<C> *rhs_ex = rhs._ex->clone();

  _ex = _ex->divide(rhs_ex);

  return *this;
}

template<typename C>
const Expression<C> &Expression<C>::operator+=(Expression<C> &&rhs)
{
  _base_expression_type<C> *rhs_ex = rhs._ex;
  rhs._ex = NULL;

  _ex = _ex->add(rhs_ex);

  return *this;
}

template<typename C>
const Expression<C> &Expression<C>::operator-=(Expression<C> &&rhs)
{
  _base_expression_type<C> *rhs_ex = rhs._ex;
  rhs._ex = NULL;

  _ex = _ex->subtract(rhs_ex);

  return *this;
}

template<typename C>
const Expression<C> &Expression<C>::operator*=(Expression<C> &&rhs)
{
  _base_expression_type<C> *rhs_ex = rhs._ex;
  rhs._ex = NULL;

  _ex = _ex->multiply(rhs_ex);

  return *this;
}

template<typename C>
const Expression<C> &Expression<C>::operator/=(Expression<C> &&rhs)
{
  _base_expression_type<C> *rhs_ex = rhs._ex;
  rhs._ex = NULL;

  _ex = _ex->divide(rhs_ex);

  return *this;
}

/**
 * @brief Check whether an expression is equivalent to a constant value.
 *
 * @tparam C is the numeric type of constants.
 * @param lhs is an expression.
 * @param rhs is a constant value.
 * @return true if and only if the first parameter is equivalent to the
 *         second one.
 */
template<typename C>
inline bool operator==(const Expression<C> &lhs, const C rhs)
{
  return (lhs._ex->get_symbol_ids().size() == 0 && lhs._ex->evaluate() == rhs);
}

/**
 * @brief Check whether an expression is equivalent to a constant value.
 *
 * @tparam C is the numeric type of constants.
 * @param lhs is an expression.
 * @param rhs is a constant value.
 * @return true if and only if the first parameter is equivalent to the
 *         second one.
 */
template<typename C>
inline bool operator==(const Expression<C> &lhs, const int rhs)
{
  return lhs == static_cast<C>(rhs);
}

/**
 * @brief Check whether an expression is equivalent to a constant value.
 *
 * @tparam C is the numeric type of constants.
 * @param lhs is a constant value.
 * @param rhs is an expression.
 * @return true if and only if the first parameter is equivalent to the
 *         second one.
 */
template<typename C>
inline bool operator==(const C lhs, const Expression<C> &rhs)
{
  return rhs == lhs;
}

/**
 * @brief Check whether an expression is equivalent to a constant value.
 *
 * @tparam C is the numeric type of constants.
 * @param lhs is a constant value.
 * @param rhs is an expression.
 * @return true if and only if the first parameter is equivalent to the
 *         second one.
 */
template<typename C>
inline bool operator==(const int lhs, const Expression<C> &rhs)
{
  return rhs == static_cast<C>(lhs);
}

/**
 * @brief Check whether an expression is not equivalent to a constant value.
 *
 * @tparam C is the numeric type of constants.
 * @param lhs is an expression.
 * @param rhs is a constant value.
 * @return true if and only if the first parameter is not equivalent to the
 *         second one.
 */
template<typename C>
inline bool operator!=(const Expression<C> &lhs, const C rhs)
{
  return !(lhs == rhs);
}

/**
 * @brief Check whether an expression is not equivalent to a constant value.
 *
 * @tparam C is the numeric type of constants.
 * @param lhs is an expression.
 * @param rhs is a constant value.
 * @return true if and only if the first parameter is not equivalent to the
 *         second one.
 */
template<typename C>
inline bool operator!=(const Expression<C> &lhs, const int rhs)
{
  return !(lhs == rhs);
}

/**
 * @brief Check whether an expression is not equivalent to a constant value.
 *
 * @tparam C is the numeric type of constants.
 * @param lhs is a constant value.
 * @param rhs is an expression.
 * @return true if and only if the first parameter is not equivalent to the
 *         second one.
 */
template<typename C>
inline bool operator!=(const C lhs, const Expression<C> &rhs)
{
  return !(lhs == rhs);
}

/**
 * @brief Check whether an expression is not equivalent to a constant value.
 *
 * @tparam C is the numeric type of constants.
 * @param lhs is a constant value.
 * @param rhs is an expression.
 * @return true if and only if the first parameter is not equivalent to the
 *         second one.
 */
template<typename C>
inline bool operator!=(const int lhs, const Expression<C> &rhs)
{
  return !(lhs == rhs);
}

/**
 * @brief Check whether an expression is greater than a constant value.
 *
 * @tparam C is the numeric type of constants.
 * @param lhs is an expression.
 * @param rhs is a constant value.
 * @return true if and only if the first parameter is constant and it is
 *         greater than the second one.
 */
template<typename C>
inline bool operator>(const Expression<C> &lhs, const C rhs)
{
  return (lhs._ex->get_symbol_ids().size() == 0 && lhs._ex->evaluate() > rhs);
}

/**
 * @brief Check whether an expression is greater than a constant value.
 *
 * @tparam C is the numeric type of constants.
 * @param lhs is an expression.
 * @param rhs is a constant value.
 * @return true if and only if the first parameter is constant and it is
 *         greater than the second one.
 */
template<typename C>
inline bool operator>(const Expression<C> &lhs, const int rhs)
{
  return lhs > static_cast<C>(rhs);
}

/**
 * @brief Check whether a constant value is greater than an expression.
 *
 * @tparam C is the numeric type of constants.
 * @param lhs is a constant value.
 * @param rhs is an expression.
 * @return true if and only if the second parameter is a constant
 *         expression and it is lesser than the first one.
 */
template<typename C>
inline bool operator>(const C lhs, const Expression<C> &rhs)
{
  return (rhs._ex.get_symbol_ids().size() == 0 && rhs._ex->evaluate() < lhs);
}

/**
 * @brief Check whether a constant value is greater than an expression.
 *
 * @tparam C is the numeric type of constants.
 * @param lhs is a constant value.
 * @param rhs is an expression.
 * @return true if and only if the second parameter is a constant
 *         expression and it is lesser than the first one.
 */
template<typename C>
inline bool operator>(const int lhs, const Expression<C> &rhs)
{
  return static_cast<int>(lhs) > rhs;
}

/**
 * @brief Check whether a constant value is greater than an expression.
 *
 * @tparam C is the numeric type of constants.
 * @param lhs is an expression.
 * @param rhs is a constant value.
 * @return true if and only if the first parameter is a constant
 *         expression and it is lesser than the second one.
 */
template<typename C>
inline bool operator<(const Expression<C> &lhs, const C rhs)
{
  return rhs > lhs;
}

/**
 * @brief Check whether a constant value is greater than an expression.
 *
 * @tparam C is the numeric type of constants.
 * @param lhs is an expression.
 * @param rhs is a constant value.
 * @return true if and only if the first parameter is a constant
 *         expression and it is lesser than the second one.
 */
template<typename C>
inline bool operator<(const Expression<C> &lhs, const int rhs)
{
  return rhs > lhs;
}

/**
 * @brief Check whether a constant value is lesser than an expression.
 *
 * @tparam C is the numeric type of constants.
 * @param lhs is a constant value.
 * @param rhs is an expression.
 * @return true if and only if the second parameter is a constant
 *         expression and it is greater than the first one.
 */
template<typename C>
inline bool operator<(const C lhs, const Expression<C> &rhs)
{
  return rhs > lhs;
}

/**
 * @brief Check whether a constant value is lesser than an expression.
 *
 * @tparam C is the numeric type of constants.
 * @param lhs is a constant value.
 * @param rhs is an expression.
 * @return true if and only if the second parameter is a constant
 *         expression and it is greater than the first one.
 */
template<typename C>
inline bool operator<(const int lhs, const Expression<C> &rhs)
{
  return rhs > lhs;
}

/**
 * @brief Sum two expressions.
 *
 * This method builds an expression that represents the sum of the two
 * expressions passed as parameters.
 *
 * @tparam C is the type of numeric constants.
 * @param lhs is an expression.
 * @param rhs is an expression.
 * @return An expression that represents the sum `lhs + rhs`.
 */
template<typename C>
Expression<C> operator+(const Expression<C> &lhs, const Expression<C> &rhs)
{
  _base_expression_type<C> *lhs_ex = lhs._ex->clone();

  _base_expression_type<C> *rhs_ex = rhs._ex->clone();

  return Expression<C>(lhs_ex->add(rhs_ex));
}

/**
 * @brief Sum an expression and a constant value.
 *
 * This method builds an expression that represents the sum of an
 * an expression and a constant value.
 *
 * @tparam C1 is the type of expression numeric constants.
 * @tparam C2 is the type of constant value.
 * @param lhs is an expression.
 * @param rhs is a constant value.
 * @return An expression that represents the sum `lhs + rhs`.
 */
template<typename C1, typename C2,
         typename
         = typename std::enable_if<std::is_arithmetic<C2>::value, C2>::type>
inline Expression<C1> operator+(const Expression<C1> &lhs, const C2 rhs)
{
  return lhs + Expression<C1>(static_cast<C1>(rhs));
}

/**
 * @brief Sum a constant value and an expression.
 *
 * This method builds an expression that represents the sum of
 * a constant value and an expression.
 *
 * @tparam C1 is the type of expression numeric constants.
 * @tparam C2 is the type of constant value.
 * @param lhs is a constant value.
 * @param rhs is an expression.
 * @return An expression that represents the sum `lhs + rhs`.
 */
template<typename C1, typename C2,
         typename
         = typename std::enable_if<std::is_arithmetic<C2>::value, C2>::type>
inline Expression<C1> operator+(const C2 lhs, const Expression<C1> &rhs)
{
  return rhs + lhs;
}

/**
 * @brief Sum two expressions.
 *
 * This method builds an expression that represents the sum of the two
 * expressions passed as parameters.
 *
 * @tparam C is the type of numeric constants.
 * @param lhs is an expression.
 * @param rhs is an expression.
 * @return An expression that represents the sum `lhs + rhs`.
 */
template<typename C>
Expression<C> operator+(Expression<C> &&lhs, const Expression<C> &rhs)
{
  _base_expression_type<C> *lhs_ex = lhs._ex;
  lhs._ex = NULL;

  _base_expression_type<C> *rhs_ex = rhs._ex->clone();

  return Expression<C>(lhs_ex->add(rhs_ex));
}

/**
 * @brief Sum two expressions.
 *
 * This method builds an expression that represents the sum of the two
 * expressions passed as parameters.
 *
 * @tparam C is the type of numeric constants.
 * @param lhs is an expression.
 * @param rhs is an expression.
 * @return An expression that represents the sum `lhs + rhs`.
 */
template<typename C>
Expression<C> operator+(const Expression<C> &lhs, Expression<C> &&rhs)
{
  _base_expression_type<C> *lhs_ex = lhs._ex->clone();

  _base_expression_type<C> *rhs_ex = rhs._ex;
  rhs._ex = NULL;

  return Expression<C>(lhs_ex->add(rhs_ex));
}

/**
 * @brief Sum two expressions.
 *
 * This method builds an expression that represents the sum of the two
 * expressions passed as parameters.
 *
 * @tparam C is the type of numeric constants.
 * @param lhs is an expression.
 * @param rhs is an expression.
 * @return An expression that represents the sum `lhs + rhs`.
 */
template<typename C>
Expression<C> operator+(Expression<C> &&lhs, Expression<C> &&rhs)
{
  _base_expression_type<C> *lhs_ex = lhs._ex;
  lhs._ex = NULL;

  _base_expression_type<C> *rhs_ex = rhs._ex;
  rhs._ex = NULL;

  return Expression<C>(lhs_ex->add(rhs_ex));
}

/**
 * @brief Subtract two expressions.
 *
 * This method builds an expression that represents the subtraction between the
 * two expressions passed as parameters.
 *
 * @tparam C is the type of numeric constants.
 * @param lhs is an expression.
 * @param rhs is an expression.
 * @return An expression that represents the subtraction `lhs - rhs`.
 */
template<typename C>
Expression<C> operator-(const Expression<C> &lhs, const Expression<C> &rhs)
{
  _base_expression_type<C> *lhs_ex = lhs._ex->clone();

  _base_expression_type<C> *rhs_ex = rhs._ex->clone();

  return Expression<C>(lhs_ex->subtract(rhs_ex));
}

/**
 * @brief Subtract a constant value from an expression.
 *
 * This method builds an expression that represents the subtraction of a
 * constant value from an expression.
 *
 * @tparam C1 is the type of expression numeric constants.
 * @tparam C2 is the type of constant value.
 * @param lhs is an expression.
 * @param rhs is a constant value.
 * @return An expression that represents the subtraction `lhs - rhs`.
 */
template<typename C1, typename C2,
         typename
         = typename std::enable_if<std::is_arithmetic<C2>::value, C2>::type>
inline Expression<C1> operator-(const Expression<C1> &lhs, const C2 rhs)
{
  return lhs - Expression<C1>(static_cast<C1>(rhs));
}

/**
 * @brief Subtract a constant value from an expression.
 *
 * This method builds an expression that represents the subtraction of a
 * constant value from an expression.
 *
 * @tparam C1 is the type of expression numeric constants.
 * @tparam C2 is the type of constant value.
 * @param lhs is a constant value.
 * @param rhs is an expression.
 * @return An expression that represents the subtraction `lhs - rhs`.
 */
template<typename C1, typename C2,
         typename
         = typename std::enable_if<std::is_arithmetic<C2>::value, C2>::type>
inline Expression<C1> operator-(const C2 lhs, const Expression<C1> &rhs)
{
  return Expression<C1>(static_cast<C1>(lhs)) - rhs;
}

/**
 * @brief Subtract two expressions.
 *
 * This method builds an expression that represents the subtraction between the
 * two expressions passed as parameters.
 *
 * @tparam C is the type of numeric constants.
 * @param lhs is an expression.
 * @param rhs is an expression.
 * @return An expression that represents the subtraction `lhs - rhs`.
 */
template<typename C>
Expression<C> operator-(Expression<C> &&lhs, const Expression<C> &rhs)
{
  _base_expression_type<C> *lhs_ex = lhs._ex;
  lhs._ex = NULL;

  _base_expression_type<C> *rhs_ex = rhs._ex->clone();

  return Expression<C>(lhs_ex->subtract(rhs_ex));
}

/**
 * @brief Subtract two expressions.
 *
 * This method builds an expression that represents the subtraction between the
 * two expressions passed as parameters.
 *
 * @tparam C is the type of numeric constants.
 * @param lhs is an expression.
 * @param rhs is an expression.
 * @return An expression that represents the subtraction `lhs - rhs`.
 */
template<typename C>
Expression<C> operator-(const Expression<C> &lhs, Expression<C> &&rhs)
{
  _base_expression_type<C> *lhs_ex = lhs._ex->clone();

  _base_expression_type<C> *rhs_ex = rhs._ex;
  rhs._ex = NULL;

  return Expression<C>(lhs_ex->subtract(rhs_ex));
}

/**
 * @brief Subtract two expressions.
 *
 * This method builds an expression that represents the subtraction between the
 * two expressions passed as parameters.
 *
 * @tparam C is the type of numeric constants.
 * @param lhs is an expression.
 * @param rhs is an expression.
 * @return An expression that represents the subtraction `lhs - rhs`.
 */
template<typename C>
Expression<C> operator-(Expression<C> &&lhs, Expression<C> &&rhs)
{
  _base_expression_type<C> *lhs_ex = lhs._ex;
  lhs._ex = NULL;

  _base_expression_type<C> *rhs_ex = rhs._ex;
  rhs._ex = NULL;

  return Expression<C>(lhs_ex->subtract(rhs_ex));
}

/**
 * @brief Mutiply two expressions.
 *
 * This method builds an expression that represents the multiplication of the
 * two expressions passed as parameters.
 *
 * @tparam C is the type of numeric constants.
 * @param lhs is an expression.
 * @param rhs is an expression.
 * @return An expression that represents the multiplication `lhs * rhs`.
 */
template<typename C>
Expression<C> operator*(const Expression<C> &lhs, const Expression<C> &rhs)
{
  _base_expression_type<C> *lhs_ex = lhs._ex->clone();

  _base_expression_type<C> *rhs_ex = rhs._ex->clone();

  return Expression<C>(lhs_ex->multiply(rhs_ex));
}

/**
 * @brief Multiply an expression and a constant value.
 *
 * This method builds an expression that represents the mutiplication
 * of an an expression and a constant value.
 *
 * @tparam C1 is the type of expression numeric constants.
 * @tparam C2 is the type of constant value.
 * @param lhs is an expression.
 * @param rhs is a constant value.
 * @return An expression that represents the multiplication `lhs * rhs`.
 */
template<typename C1, typename C2,
         typename
         = typename std::enable_if<std::is_arithmetic<C2>::value, C2>::type>
inline Expression<C1> operator*(const Expression<C1> &lhs, const C2 rhs)
{
  return lhs * Expression<C1>(static_cast<C1>(rhs));
}

/**
 * @brief Multiply a constant value and an expression.
 *
 * This method builds an expression that represents the mutiplication
 * of an an expression and a constant value.
 *
 * @tparam C1 is the type of expression numeric constants.
 * @tparam C2 is the type of constant value.
 * @param lhs is a constant value.
 * @param rhs is an expression.
 * @return An expression that represents the multiplication `lhs * rhs`.
 */
template<typename C1, typename C2,
         typename
         = typename std::enable_if<std::is_arithmetic<C2>::value, C2>::type>
inline Expression<C1> operator*(const C2 lhs, const Expression<C1> &rhs)
{
  return rhs * lhs;
}

/**
 * @brief Mutiply two expressions.
 *
 * This method builds an expression that represents the multiplication of the
 * two expressions passed as parameters.
 *
 * @tparam C is the type of numeric constants.
 * @param lhs is an expression.
 * @param rhs is an expression.
 * @return An expression that represents the multiplication `lhs * rhs`.
 */
template<typename C>
Expression<C> operator*(Expression<C> &&lhs, const Expression<C> &rhs)
{
  _base_expression_type<C> *lhs_ex = lhs._ex;
  lhs._ex = NULL;

  _base_expression_type<C> *rhs_ex = rhs._ex->clone();

  return Expression<C>(lhs_ex->multiply(rhs_ex));
}

/**
 * @brief Mutiply two expressions.
 *
 * This method builds an expression that represents the multiplication of the
 * two expressions passed as parameters.
 *
 * @tparam C is the type of numeric constants.
 * @param lhs is an expression.
 * @param rhs is an expression.
 * @return An expression that represents the multiplication `lhs * rhs`.
 */
template<typename C>
Expression<C> operator*(const Expression<C> &lhs, Expression<C> &&rhs)
{
  _base_expression_type<C> *lhs_ex = lhs._ex->clone();

  _base_expression_type<C> *rhs_ex = rhs._ex;
  rhs._ex = NULL;

  return Expression<C>(lhs_ex->multiply(rhs_ex));
}

/**
 * @brief Mutiply two expressions.
 *
 * This method builds an expression that represents the multiplication of the
 * two expressions passed as parameters.
 *
 * @tparam C is the type of numeric constants.
 * @param lhs is an expression.
 * @param rhs is an expression.
 * @return An expression that represents the multiplication `lhs * rhs`.
 */
template<typename C>
Expression<C> operator*(Expression<C> &&lhs, Expression<C> &&rhs)
{
  _base_expression_type<C> *lhs_ex = lhs._ex;
  lhs._ex = NULL;

  _base_expression_type<C> *rhs_ex = rhs._ex;
  rhs._ex = NULL;

  return Expression<C>(lhs_ex->multiply(rhs_ex));
}

/**
 * @brief Divide two expressions.
 *
 * This method builds an expression that represents the division between the
 * two expressions passed as parameters.
 *
 * @tparam C is the type of numeric constants.
 * @param lhs is an expression.
 * @param rhs is an expression.
 * @return An expression that represents the division `lhs / rhs`.
 */
template<typename C>
Expression<C> operator/(const Expression<C> &lhs, const Expression<C> &rhs)
{
  _base_expression_type<C> *lhs_ex = lhs._ex->clone();

  _base_expression_type<C> *rhs_ex = rhs._ex->clone();

  return Expression<C>(lhs_ex->divide(rhs_ex));
}

/**
 * @brief Divide an expression by a constant value.
 *
 * This method builds an expression that represents the division among
 * an expression and a constat value.
 *
 * @tparam C1 is the type of expression numeric constants.
 * @tparam C2 is the type of constant value.
 * @param lhs is an expression.
 * @param rhs is a constant value.
 * @return An expression that represents the division `lhs / rhs`.
 */
template<typename C1, typename C2,
         typename
         = typename std::enable_if<std::is_arithmetic<C2>::value, C2>::type>
inline Expression<C1> operator/(const Expression<C1> &lhs, const C2 rhs)
{
  return lhs / Expression<C1>(static_cast<C1>(rhs));
}

/**
 * @brief Divide a constant value by an expression.
 *
 * This method builds an expression that the division among a
 * constant value and an expression.
 *
 * @tparam C1 is the type of expression numeric constants.
 * @tparam C2 is the type of constant value.
 * @param lhs is a constant value.
 * @param rhs is an expression.
 * @return An expression that represents the subtraction `lhs - rhs`.
 */
template<typename C1, typename C2,
         typename
         = typename std::enable_if<std::is_arithmetic<C2>::value, C2>::type>
inline Expression<C1> operator/(const C2 lhs, const Expression<C1> &rhs)
{
  return Expression<C1>(static_cast<C1>(lhs)) / rhs;
}

/**
 * @brief Divide two expressions.
 *
 * This method builds an expression that represents the division between the
 * two expressions passed as parameters.
 *
 * @tparam C is the type of numeric constants.
 * @param lhs is an expression.
 * @param rhs is an expression.
 * @return An expression that represents the division `lhs / rhs`.
 */
template<typename C>
Expression<C> operator/(Expression<C> &&lhs, const Expression<C> &rhs)
{
  _base_expression_type<C> *lhs_ex = lhs._ex;
  lhs._ex = NULL;

  _base_expression_type<C> *rhs_ex = rhs._ex->clone();

  return Expression<C>(lhs_ex->divide(rhs_ex));
}

/**
 * @brief Divide two expressions.
 *
 * This method builds an expression that represents the division between the
 * two expressions passed as parameters.
 *
 * @tparam C is the type of numeric constants.
 * @param lhs is an expression.
 * @param rhs is an expression.
 * @return An expression that represents the division `lhs / rhs`.
 */
template<typename C>
Expression<C> operator/(const Expression<C> &lhs, Expression<C> &&rhs)
{
  _base_expression_type<C> *lhs_ex = lhs._ex->clone();

  _base_expression_type<C> *rhs_ex = rhs._ex;
  rhs._ex = NULL;

  return Expression<C>(lhs_ex->divide(rhs_ex));
}

/**
 * @brief Divide two expressions.
 *
 * This method builds an expression that represents the division between the
 * two expressions passed as parameters.
 *
 * @tparam C is the type of numeric constants.
 * @param lhs is an expression.
 * @param rhs is an expression.
 * @return An expression that represents the division `lhs / rhs`.
 */
template<typename C>
Expression<C> operator/(Expression<C> &&lhs, Expression<C> &&rhs)
{
  _base_expression_type<C> *lhs_ex = lhs._ex;
  lhs._ex = NULL;

  _base_expression_type<C> *rhs_ex = rhs._ex;
  rhs._ex = NULL;

  return Expression<C>(lhs_ex->divide(rhs_ex));
}

/**
 * @brief Sum an expression and a constant.
 *
 * This method builds an expression that represents the sum of
 * an expression and a constant.
 *
 * @tparam C is the type of numeric constants.
 * @param lhs is an expression.
 * @param rhs is a constant value.
 * @return An expression that represents the sum `lhs + rhs`.
 */
template<typename T>
Expression<T> operator+(const Expression<T> &lhs, const T rhs)
{
  if (rhs == 0) {
    return lhs;
  }

  return lhs + Expression<T>(rhs);
}

/**
 * @brief Subtract a constant from an expression.
 *
 * This method builds an expression that represents the subtraction of
 * a constant from an expression
 *
 * @tparam C is the type of numeric constants.
 * @param lhs is an expression.
 * @param rhs is a constant value.
 * @return An expression that represents the sum `lhs - rhs`.
 */
template<typename T>
Expression<T> operator-(const Expression<T> &lhs, const T rhs)
{
  if (rhs == 0) {
    return lhs;
  }

  return lhs - Expression<T>(rhs);
}

/**
 * @brief Multiply an expression and a constant.
 *
 * This method builds an expression that represents the multiplication between
 * an expression and a constant.
 *
 * @tparam C is the type of numeric constants.
 * @param lhs is an expression.
 * @param rhs is a constant value.
 * @return An expression that represents the multiplication `lhs * rhs`.
 */
template<typename T>
Expression<T> operator*(const Expression<T> &lhs, const T rhs)
{
  if (rhs == 1) {
    return lhs;
  }

  if (rhs == 0) {
    return Expression<T>(rhs);
  }

  return lhs * Expression<T>(rhs);
}
/**
 * @brief Divide an expression by a constant.
 *
 * This method builds an expression that represents the division between an
 * expression and a constant value.
 *
 * @tparam C is the type of numeric constants.
 * @param lhs is an expression.
 * @param rhs is a constant value.
 * @return An expression that represents the division `lhs / rhs`.
 */
template<typename T>
Expression<T> operator/(const Expression<T> &lhs, const T rhs)
{
  if (rhs == 1) {
    return lhs;
  }

  return lhs / Expression<T>(rhs);
}

/**
 * @brief Sum a constant and an expression.
 *
 * This method builds an expression that represents the sum of
 * a constant and an expression.
 *
 * @tparam C is the type of numeric constants.
 * @param lhs is a constant value.
 * @param rhs is an expression.
 * @return An expression that represents the sum `lhs + rhs`.
 */
template<typename T>
Expression<T> operator+(const T lhs, const Expression<T> &rhs)
{
  if (lhs == 0) {
    return rhs;
  }

  return Expression<T>(lhs) + rhs;
}

/**
 * @brief Subtract an expression from a constant.
 *
 * This method builds an expression that represents the subtraction of
 * an expression from a constant.
 *
 * @tparam C is the type of numeric constants.
 * @param lhs is a constant value.
 * @param rhs is an expression.
 * @return An expression that represents the subtraction `lhs - rhs`.
 */
template<typename T>
Expression<T> operator-(const T lhs, const Expression<T> &rhs)
{
  if (lhs == 0) {
    return rhs;
  }

  return Expression<T>(lhs) - rhs;
}
/**
 * @brief Multiply a constant and an expression.
 *
 * This method builds an expression that represents the multiplication of a
 * constant value and an expression.
 *
 * @tparam C is the type of numeric constants.
 * @param lhs is a constant value.
 * @param rhs is an expression.
 * @return An expression that represents the multiplication `lhs * rhs`.
 */
template<typename T>
Expression<T> operator*(const T lhs, const Expression<T> &rhs)
{
  if (lhs == 1) {
    return rhs;
  }

  if (lhs == 0) {
    return Expression<T>(lhs);
  }

  return Expression<T>(lhs) * rhs;
}

/**
 * @brief Divide a constant by an expression.
 *
 * This method builds an expression that represents the division between a
 * constant value and an expression.
 *
 * @tparam C is the type of numeric constants.
 * @param lhs is a constant value.
 * @param rhs is an expression.
 * @return An expression that represents the division `lhs / rhs`.
 */
template<typename T>
Expression<T> operator/(const T lhs, const Expression<T> &rhs)
{
  if (lhs == 0) {
    return Expression<T>(lhs);
  }

  return Expression<T>(lhs) / rhs;
}

template<typename C>
Symbol<C>::Symbol(): Expression<C>()
{
}

template<typename C>
Symbol<C>::Symbol(const char *name): Symbol<C>(std::string(name))
{
}

template<typename C>
Symbol<C>::Symbol(const std::string &name): Expression<C>()
{
  auto found_symb = _declared_symbols.find(name);

  if (found_symb == std::end(_declared_symbols)) {
    SymbolIdType symb_id = _symbol_names.size();
    _symbol_names.push_back(name);
    _declared_symbols[name] = symb_id;
    this->_ex = new _symbol_type<C>(symb_id);
  } else {
    this->_ex = new _symbol_type<C>(found_symb->second);
  }
}

template<typename C>
Symbol<C>::Symbol(const Symbol<C> &orig):
    Expression<C>((orig._ex == NULL ? NULL : orig._ex->clone()))
{
}

template<typename C>
constexpr bool operator<(const Symbol<C> &a, const Symbol<C> &b)
{
  return a.get_id() < b.get_id();
}

}

/*! @} End of SymbolicAlgebra group */

namespace std
{
/**
 * @brief Swap the content of two Expression objects.
 *
 * @tparam C is the type of numeric constants.
 * @param a is an expression.
 * @param b is an expression.
 */
template<typename C>
void swap(SymbolicAlgebra::Expression<C> &a, SymbolicAlgebra::Expression<C> &b)
{
  std::swap(a._ex, b._ex);
}

/**
 * @brief Write an expression in an output stream.
 *
 * @tparam C is the type of numeric constants.
 * @param os is the output stream.
 * @param ex is the expression to be printed.
 * @return the output stream.
 */
template<typename C>
std::ostream &operator<<(std::ostream &os,
                         const SymbolicAlgebra::Expression<C> &ex)
{
  if (ex._ex != NULL) {
    ex._ex->print(os);
  }

  return os;
}
} // namespace std

#endif // SYMBOLIC_ALGEBRA_H_