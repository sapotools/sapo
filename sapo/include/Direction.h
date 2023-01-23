#ifndef __DIRECTION_H__
#define __DIRECTION_H__

#include <cmath>
#include <limits>

#include "SymbolicAlgebra.h"
#include "Expr.h"
#include "AbsSynIO.h"

namespace AbsSyn
{

template<typename T=double>
class Direction
{
public:
  enum Type {
    LT, // <
    LE, // <=
    GT, // >
    GE, // >=
    EQ, // =
    IN  // lhs in [a,b]
  };

protected:
  SymbolicAlgebra::Expression<T> _linear_expression;
  const SymbolicAlgebra::Symbol<T> *_symbol;
  T _lower_bound;
  T _upper_bound;
  Type _type;

private:
  Direction(const SymbolicAlgebra::Expression<T>& linear_expression,
            const SymbolicAlgebra::Symbol<T> *symbol,
            const T& lower_bound, const T& upper_bound,
            const Type& type):
    _linear_expression(linear_expression), _symbol(symbol),
    _lower_bound(lower_bound), _upper_bound(upper_bound),
    _type(type)
  {}

  static Type get_complementary(const Type type)
  {
    switch (type) {
      case Type::LT:
        return Type::GT;
      case Type::LE:
        return Type::GE;
      case Type::GT:
        return Type::LT;
      case Type::GE:
        return Type::LE;
      case Type::EQ:
      case Type::IN:
        return type;
      default:
        throw std::domain_error("Unknown Type");
    }
  }
public:

  Direction(const Direction<T>& orig):
    _linear_expression(orig._linear_expression),
    _symbol(orig._symbol), _lower_bound(orig._lower_bound), 
    _upper_bound(orig._upper_bound), _type(orig._type)
  {}

  Direction(Direction<T>&& orig):
    _linear_expression(std::move(orig._linear_expression)),
    _symbol(orig._symbol), _lower_bound(orig._lower_bound), 
    _upper_bound(orig._upper_bound), _type(orig._type)
  {}

  Direction(const SymbolicAlgebra::Expression<T>& linear_expression, 
            const T& lower_bound, const T& upper_bound,
            const SymbolicAlgebra::Symbol<T> *symbol = nullptr):
      _linear_expression(linear_expression), _symbol(symbol),
      _lower_bound(lower_bound), _upper_bound(upper_bound),
      _type(Type::IN)
  {
    auto const_term = linear_expression.get_constant_term();
    _linear_expression -= const_term;

    const T value_const_term = const_term.evaluate();
    _lower_bound -= value_const_term;
    _upper_bound -= value_const_term;
  }

  /*
  Direction(SymbolicAlgebra::Expression<> e1, SymbolicAlgebra::Expression<> e2,
            Type t, double lb, double ub, SymbolicAlgebra::Symbol<> *sym):
      lhs(e1),
      rhs(e2), type(t), LB(lb), UB(ub), s(sym)
  {
  }
  */

  Direction(const SymbolicAlgebra::Expression<T>& linear_expression1, 
            const SymbolicAlgebra::Expression<T>& linear_expression2,
            const Type type):
      _linear_expression(linear_expression1 - linear_expression2),
      _symbol(nullptr), 
      _lower_bound(-std::numeric_limits<T>::infinity()),
      _upper_bound(std::numeric_limits<T>::infinity()),
      _type(type)
  {
    auto const_term = _linear_expression.get_constant_term();
    _linear_expression -= const_term;  // e1-e2 == _linear_expression+const_term

    const T value_const_term = const_term.evaluate();

    switch (type) {
    case Type::LE:     // e1 <= e2: _linear_expression <= -const_term  
    case Type::LT:     // e1 < e2: _linear_expression < -const_term  
      _upper_bound = -value_const_term;
      break;
    case Type::GE:    // e1 >= e2: _linear_expression >= -const_term    
    case Type::GT:    // e1 > e2: _linear_expression > -const_term  
      _lower_bound = -value_const_term;
      break;
    case Type::EQ:    // e1 == e2
      _upper_bound = -value_const_term;
      _lower_bound = _upper_bound;
      break;
    case Type::IN:
      throw std::domain_error("Unhandled direction type (IN)");
      break;
    default:
      throw std::domain_error("Unknown direction type");
    }
  }

  inline
  Type get_type() const
  {
    return _type;
  }

  inline
  const SymbolicAlgebra::Expression<T> &
  get_linear_expression() const
  {
    return _linear_expression;
  }

  ~Direction() {}

  inline
  std::set<SymbolicAlgebra::Symbol<T>> get_variables() const
  {
    return _linear_expression.get_symbols();
  }

  std::vector<T> 
  get_variable_coefficients(const std::vector<SymbolicAlgebra::Symbol<T>>& variables) const
  {
    std::vector<T> res;

    for (const auto& var: variables) {
      res.push_back(getCoefficient(_linear_expression, var));
    }

    return res;
  }

  std::string getName() const
  {
    if (_symbol == nullptr) {
      throw std::runtime_error("Direction has no name");
    }
    return SymbolicAlgebra::Symbol<T>::get_symbol_name(_symbol->get_id());
  }

  inline
  const SymbolicAlgebra::Symbol<T> *getSymbol() const
  {
    return _symbol;
  }

  inline
  void setSymbol(const SymbolicAlgebra::Symbol<T> *symbol)
  {
    _symbol = symbol;
  }

  inline T get_lower_bound() const 
  {
    return _lower_bound;
  }

  inline T get_upper_bound() const 
  {
    return _upper_bound;
  }
  
  inline bool has_lower_bound() const
  {
    return _lower_bound != -std::numeric_limits<T>::infinity();
  }

  inline bool has_upper_bound() const
  {
    return _upper_bound != std::numeric_limits<T>::infinity();
  }

  void set_lower_bound(const T& lower_bound)
  {
    if (this->has_upper_bound()) {
      _type = Type::IN;
    }

    _lower_bound = (lower_bound == 0 ? 0 : lower_bound);
  }

  void set_upper_bound(const T& upper_bound)
  {
    if (this->has_lower_bound()) {
      _type = Type::IN;
    }

    _upper_bound = (upper_bound == 0 ? 0 : upper_bound);
  }

  inline
  Type getType() const
  {
    return _type;
  }

  inline
  bool contains(const std::vector<SymbolicAlgebra::Symbol<T>>& symbols) const
  {
    const auto dir_symbols = _linear_expression.get_symbols();

    std::less<SymbolicAlgebra::Symbol<T>> cmp;
    return std::includes(std::begin(dir_symbols), std::end(dir_symbols),
                         std::begin(symbols), std::end(symbols), cmp);
  }

  //Direction<T> *copy() const; // deep copy of direction

  inline
  Direction<T> get_complementary() const
  {
    return Direction<T>(-_linear_expression, _symbol, -_upper_bound, -_lower_bound, 
                        Direction<T>::get_complementary(_type));
  }

  // checks if the symbol named "name" is present in the direction
  inline bool covers(const Symbol<T> &s) const
  {
    return _linear_expression.degree(s) != 0;
  }
};


template<typename T>
std::ostream &operator<<(std::ostream &os, const Direction<T> &direction)
{
  if (direction.getSymbol() != nullptr) {
    os << *(direction.getSymbol()) << ": ";
  }

  os << direction.get_linear_expression();

  switch (direction.get_type()) {
  case Direction<T>::Type::LT:
    return os << " < "<< direction.get_upper_bound();
  case Direction<T>::Type::LE:
    return os << " <= "<< direction.get_upper_bound();
  case Direction<T>::Type::GT:
    return os << " > " << direction.get_lower_bound();
  case Direction<T>::Type::GE:
    return os << " >= " << direction.get_lower_bound();
  case Direction<T>::Type::EQ:
    return os << " == " << direction.get_lower_bound();
  case Direction<T>::Type::IN:
    return os << " in [" << direction.get_lower_bound() << "," 
                         << direction.get_upper_bound() << "]";
  default:
    throw std::logic_error("unsupported direction type");
  }
  return os;
}
}

#endif
