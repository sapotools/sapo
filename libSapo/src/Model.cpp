/**
 * @file Model.cpp
 * Model implementation
 *
 * @author Alberto Casagrande <acasagrande@units.it>
 * @version 0.1
 */

#include "Model.h"

Model::Model(const std::vector<SymbolicAlgebra::Symbol<>> &variables,
             const std::vector<SymbolicAlgebra::Expression<>> &dynamics,
             const Bundle &init_set, const std::string name):
    Model(variables, {}, dynamics, init_set, PolytopesUnion(), name)
{
}

Model::Model(const std::vector<SymbolicAlgebra::Symbol<>> &variables,
             const std::vector<SymbolicAlgebra::Symbol<>> &parameters,
             const std::vector<SymbolicAlgebra::Expression<>> &dynamics,
             const Bundle &init_set, const PolytopesUnion &parameter_set,
             const std::string name):
    _vars(variables),
    _params(parameters), _dynamics(dynamics),
    _init_set(std::make_shared<Bundle>(init_set)), _param_set(parameter_set),
    _spec(nullptr), _assumptions(), _name(name)
{
  using namespace SymbolicAlgebra;
  std::set<Symbol<>> dyn_symbols;

  for (const auto &dyn: dynamics) {
    auto symbols = dyn.get_symbols();
    dyn_symbols.insert(std::begin(symbols), std::end(symbols));
  }

  for (const auto &var: variables) {
    dyn_symbols.erase(var);
  }

  for (const auto &param: parameters) {
    dyn_symbols.erase(param);
  }

  if (!dyn_symbols.empty()) {
    std::ostringstream oss;
    oss << "The dynamic law symbols " << dyn_symbols
        << "are neither variables nor parameters";
    throw std::domain_error(oss.str());
  }

  if (init_set.dim() != variables.size()) {
    std::ostringstream oss;
    oss << "The space dimension of the initial set (" << init_set.dim()
        << ") differs from the "
        << "variable number (" << variables.size() << "). "
        << "They must be the same.";
    throw std::domain_error(oss.str());
  }

  if (parameter_set.dim() != parameters.size()) {
    std::ostringstream oss;
    oss << "The space dimension of the parameter set (" << parameter_set.dim()
        << ") differs from the "
        << "parameter number (" << parameters.size() << "). "
        << "They must be the same.";
    throw std::domain_error(oss.str());
  }
}

Model &Model::set_specification(std::shared_ptr<STL::STL> specification)
{
  if (specification != nullptr) {
    auto spec_vars = specification->get_variables();

    for (const auto &var: _vars) {
      spec_vars.erase(var);
    }

    if (!spec_vars.empty()) {
      std::ostringstream oss;
      oss << "The specification variables " << spec_vars
          << "are not among the model variables";
      throw std::domain_error(oss.str());
    }

    _spec = specification;
  }

  return *this;
}

Model &Model::set_assumptions(const LinearSystem &assumptions)
{
  if (assumptions.size() > 0 && assumptions.dim() != _vars.size()) {
    throw std::domain_error("The number of dimensions of the "
                            "assumptions space differs from "
                            "that of the model variables.");
  }

  _assumptions = LinearSystem(assumptions);

  return *this;
}