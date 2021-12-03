#ifndef SYMBOLIC_ALGEBRA_H_
#define SYMBOLIC_ALGEBRA_H_

#include <iostream>
#include <list>
#include <map>
#include <vector>

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
template<typename C>
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
template<typename C>
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
  Expression(const C &value);

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
  inline C evaluate() const
  {
    return _ex->evaluate();
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
    case FINITE_SUM:
      return op->add(this);
    case CONSTANT:
      if (((_constant_type<C> *)op)->get_value() == 0) {
        delete op;
        return this;
      }
    default:
      _finite_sum_type<C> *sum = new _finite_sum_type<C>();

      sum->add(this);
      return sum->add(op);
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
    _constant_type<C> *minus1 = new _constant_type<C>(-1);
    return this->add(minus1->multiply(op));
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
    case FINITE_PROD:
      return op->multiply(this);

    case CONSTANT:
      if (((_constant_type<C> *)op)->get_value() == 1) {
        delete op;

        return this;
      }
    default:
      _finite_prod_type<C> *prod = new _finite_prod_type<C>();

      prod->multiply(this);
      return prod->multiply(op);
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
    if (op->type() == CONSTANT
        && ((_constant_type<C> *)op)->get_value() == 1) {
      delete op;

      return this;
    }

    _finite_prod_type<C> *prod = new _finite_prod_type<C>();

    prod->multiply(this);
    return prod->divide(op);
  }

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
  virtual std::list<_finite_prod_type<C> *>
  multiply(std::list<_finite_prod_type<C> *> &prods)
  {
    std::list<_finite_prod_type<C> *> result;

    for (auto p_it = std::begin(prods); p_it != std::end(prods); ++p_it) {
      result.push_back(
          (_finite_prod_type<C> *)((*p_it)->multiply(this->clone())));
    }

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
   * @brief Print the expression in an output stream.
   *
   * @param os is the output stream in which the expression must be printed.
   */
  virtual void print(std::ostream &os) const = 0;

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
   * @brief Numerically evaluate the expression.
   *
   * @return The numeric evaluation of the expression.
   */
  C evaluate() const
  {
    return _value;
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
  std::list<_base_expression_type<C> *>
      _sum; //!< The finite list of summed expressions.
public:
  typedef typename Symbol<C>::SymbolIdType SymbolIdType;

  /**
   * @brief Build an empty finite sum.
   */
  _finite_sum_type(): _base_expression_type<C>(), _sum() {}

  /**
   * @brief Copy constructor.
   *
   * @param orig is the model for the new finite sum.
   */
  _finite_sum_type(const _finite_sum_type<C> &orig):
      _base_expression_type<C>(), _sum()
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
  _finite_sum_type(std::list<_finite_prod_type<C> *> &model):
      _base_expression_type<C>(), _sum()
  {
    for (auto it = std::begin(model); it != std::end(model); ++it) {
      _sum.push_back((_base_expression_type<C> *)*it);
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
      _sum.splice(std::end(_sum), op_sum->_sum);

      delete op;

      break;
    }
    case CONSTANT: {
      _constant_type<C> *op_const = (_constant_type<C> *)op;
      if (op_const->get_value() == 0) {
        delete op;

        return this;
      }
    }
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
  std::list<_finite_prod_type<C> *>
  multiply(std::list<_finite_prod_type<C> *> &prods)
  {
    std::list<_finite_prod_type<C> *> result;

    for (auto s_it = std::begin(_sum); s_it != std::end(_sum); ++s_it) {
      for (auto p_it = std::begin(prods); p_it != std::end(prods); ++p_it) {
        _finite_prod_type<C> *new_prod = new _finite_prod_type<C>(*(*p_it));

        _base_expression_type<C> *new_s = new_prod->multiply((*s_it)->clone());

        result.push_back((_finite_prod_type<C> *)new_s);
      }
      delete *s_it;
    }
    _sum.clear();

    for (auto p_it = std::begin(prods); p_it != std::end(prods); ++p_it) {
      delete *p_it;
    }

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
  _base_expression_type<C> *expand() const
  {
    std::list<_base_expression_type<C> *> new_list;
    for (auto it = std::begin(_sum); it != std::end(_sum); ++it) {
      _base_expression_type<C> *new_it = (*it)->expand();

      if (new_it->type() == FINITE_SUM) {
        new_list.splice(std::end(new_list),
                        ((_finite_sum_type<C> *)new_it)->_sum);

        delete new_it;
      } else {
        new_list.push_back(new_it);
      }
    }

    _finite_sum_type<C> *new_sum = new _finite_sum_type<C>();

    std::swap(new_sum->_sum, new_list);

    return new_sum;
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
    for (auto it = std::begin(_sum); it != std::end(_sum); ++it) {
      *it = (*it)->replace(replacements);
    }

    return this;
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
   * @brief Print the expression in an output stream.
   *
   * @param os is the output stream in which the expression must be printed.
   */
  void print(std::ostream &os) const
  {
    for (auto it = std::begin(_sum); it != std::end(_sum); ++it) {
      if (it != std::begin(_sum)) {
        os << " + ";
      }
      switch ((*it)->type()) {
      case FINITE_SUM:
        os << "(";
        (*it)->print(os);
        os << ")";
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
      throw std::runtime_error("Non-constant denominator not supported");
    }

    std::list<_finite_prod_type<C> *> result;

    _finite_prod_type<C> *base = new _finite_prod_type<C>();

    base->_constant = base_const;

    result.push_back(base);

    for (auto it = std::begin(_numerator); it != std::end(_numerator); ++it) {
      _base_expression_type<C> *new_it = (*it)->expand();

      result = new_it->multiply(result);

      delete new_it;
    }
    return new _finite_sum_type(result);
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
    _finite_prod_type<C> *new_prod = new _finite_prod_type<C>();
    new_prod->_constant = this->_constant;

    for (auto it = std::begin(_numerator); it != std::end(_numerator); ++it) {
      new_prod->multiply((*it)->replace(replacements));
    }
    _numerator.clear();

    for (auto it = std::begin(_denominator); it != std::end(_denominator);
         ++it) {
      new_prod->divide((*it)->replace(replacements));
    }
    _denominator.clear();

    std::swap(_numerator, new_prod->_numerator);
    std::swap(_denominator, new_prod->_denominator);
    std::swap(_constant, new_prod->_constant);

    delete new_prod;

    return this;
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
   * @brief Clone the current object.
   *
   * @return A clone of the current object.
   */
  _base_expression_type<C> *clone() const
  {
    return new _symbol_type(*this);
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
   * @brief Print the expression in an output stream.
   *
   * @param os is the output stream in which the expression must be printed.
   */
  void print(std::ostream &os) const
  {
    os << Symbol<C>::get_symbol_name(_id);
  }
};

template<typename C>
Expression<C>::Expression(_base_expression_type<C> *_ex): _ex(_ex)
{
}

template<typename C>
Expression<C>::Expression(): _ex(NULL)
{
}

template<typename C>
Expression<C>::Expression(const C &value): _ex(new _constant_type(value))
{
}

template<typename C>
Expression<C>::Expression(const Expression<C> &orig): _ex(orig._ex->clone())
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
Symbol<C>::Symbol(const std::string &name): Expression<C>()
{
  auto found_symb = _declared_symbols.find(name);

  if (found_symb == std::end(_declared_symbols)) {
    SymbolIdType symb_id = _symbol_names.size();
    _symbol_names.push_back(name);
    _declared_symbols[name] = symb_id;
    this->_ex = new _symbol_type<C>(symb_id);
  }
}

template<typename C>
Symbol<C>::Symbol(const Symbol<C> &orig): Expression<C>(orig._ex->clone())
{
}

template<typename C>
constexpr bool operator<(const Symbol<C> &a, const Symbol<C> &b)
{
  return a.get_id() < b.get_id();
}

} // namespace SymbolicAlgebra

namespace std
{
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