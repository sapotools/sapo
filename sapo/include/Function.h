#ifndef __FUNCTION_H__
#define __FUNCTION_H__

#include <functional>

#include <SymbolicAlgebra.h>
#include <ErrorHandling.h>

namespace AbsSyn
{

template<typename T>
class Function
{
  std::string _name; //!< function name
  std::vector<SymbolicAlgebra::Symbol<T>>
      _formal_parameters;                     //!< formal parameters
  SymbolicAlgebra::Expression<T> _expression; //!< expression

public:
  Function(const std::string &name,
           const std::vector<SymbolicAlgebra::Symbol<T>> &formal_parameters,
           const SymbolicAlgebra::Expression<T> expression):
      _name(name),
      _formal_parameters(formal_parameters), _expression(expression)
  {
  }

  ~Function() {}

  inline const std::string &name() const
  {
    return _name;
  }

  inline const std::vector<SymbolicAlgebra::Symbol<T>> &
  formal_parameters() const
  {
    return _formal_parameters;
  }

  inline const SymbolicAlgebra::Expression<T> &expression() const
  {
    return _expression;
  }

  SymbolicAlgebra::Expression<T> operator()(
      const std::vector<SymbolicAlgebra::Expression<T>> &actual_parameters)
      const
  {
    using namespace SymbolicAlgebra;

    if (actual_parameters.size() != _formal_parameters.size()) {
      SAPO_ERROR("Wrong number of parameters", std::domain_error);
    }

    typename Expression<T>::replacement_type replacement;

    for (size_t i = 0; i < _formal_parameters.size(); ++i) {
      replacement[_formal_parameters[i]] = actual_parameters[i];
    }

    return SymbolicAlgebra::Expression<T>(_expression).replace(replacement);
  }
};

template<typename T>
std::ostream &operator<<(std::ostream &os, const Function<T> &function)
{
  os << function.name() << "(";
  auto sep = "";
  for (const auto &parameter: function.formal_parameters()) {
    os << sep << parameter;
    sep = ",";
  }
  os << ") = " << function.expression();
}

typedef std::pair<std::string, size_t> FunctionSignature;

struct SignatureOrder {
  bool operator()(const FunctionSignature &a, const FunctionSignature &b) const
  {
    return ((a.second < b.second)
            || ((a.second == b.second) && (a.first < b.first)));
  }
};

}

#endif // __FUNCTION_H__
