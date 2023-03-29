/**
 * @file Model.cpp
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @author Alberto Casagrande <acasagrande@units.it>
 * @brief Discrete-time dynamical models
 * @version 0.1
 * @date 2015-10-15
 *
 * @copyright Copyright (c) 2015-2022
 */

#include "Model.h"

#include "ErrorHandling.h"

Model::Model(const Bundle<double> &init_set, const std::string name):
    Model(init_set, SetsUnion<Polytope<double>>(), name)
{
}

Model::Model(const Bundle<double> &init_set,
             const SetsUnion<Polytope<double>> &parameter_set,
             const std::string name):
    _init_set(std::make_shared<Bundle<double>>(init_set)),
    _param_set(parameter_set), _spec(nullptr), _assumptions(), _invariant(),
    _name(name)
{
}

void Model::validate_parameters(
    const std::vector<SymbolicAlgebra::Symbol<double>> &variables,
    const std::vector<SymbolicAlgebra::Symbol<double>> &parameters,
    const std::vector<SymbolicAlgebra::Expression<double>> &dynamics,
    const Bundle<double> &init_set,
    const SetsUnion<Polytope<double>> &parameter_set)
{
  using namespace SymbolicAlgebra;
  std::set<Symbol<double>> dyn_symbols;

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
    oss << "the dynamic law symbols " << dyn_symbols
        << "are neither variables nor parameters";
    SAPO_ERROR(oss.str(), std::domain_error);
  }

  if (init_set.dim() != variables.size()) {
    std::ostringstream oss;
    oss << "the space dimension of the initial set (" << init_set.dim()
        << ") differs from the variable number (" << variables.size() << ")";
    SAPO_ERROR(oss.str(), std::domain_error);
  }

  if (parameter_set.dim() != parameters.size()) {
    std::ostringstream oss;
    oss << "the space dimension of the parameter set (" << parameter_set.dim()
        << ") differs from the parameter number (" << parameter_set.dim()
        << ")";
    SAPO_ERROR(oss.str(), std::domain_error);
  }
}

Model &Model::set_specification(std::shared_ptr<STL::STL> specification)
{
  if (specification != nullptr) {
    auto spec_vars = specification->get_variables();

    for (const auto &var: variables()) {
      spec_vars.erase(var);
    }

    if (!spec_vars.empty()) {
      std::ostringstream oss;
      oss << "the specification variables " << spec_vars
          << "are not among the model variables";
      SAPO_ERROR(oss.str(), std::domain_error);
    }

    _spec = specification;
  }

  return *this;
}

Model &Model::set_assumptions(const LinearSystem<double> &assumptions)
{
  if (assumptions.size() > 0 && assumptions.dim() != dim()) {
    SAPO_ERROR("the assumptions space and the model variables differ "
               "in the number of dimensions",
               std::domain_error);
  }

  _assumptions = LinearSystem<double>(assumptions);

  return *this;
}

Model &Model::set_invariant(const LinearSystem<double> &invariant)
{
  if (invariant.size() > 0 && invariant.dim() != dim()) {
    SAPO_ERROR("the invariant space and the model variables differ "
               "in the number of dimensions",
               std::domain_error);
  }

  _invariant = LinearSystem<double>(invariant);

  return *this;
}

DiscreteModel::DiscreteModel(
    const std::vector<SymbolicAlgebra::Symbol<double>> &variables,
    const std::vector<SymbolicAlgebra::Expression<double>> &dynamics,
    const Bundle<double> &init_set, const std::string name):
    DiscreteModel(variables, {}, dynamics, init_set, {}, name)
{
}

DiscreteModel::DiscreteModel(
    const std::vector<SymbolicAlgebra::Symbol<double>> &variables,
    const std::vector<SymbolicAlgebra::Symbol<double>> &parameters,
    const std::vector<SymbolicAlgebra::Expression<double>> &dynamics,
    const Bundle<double> &init_set,
    const SetsUnion<Polytope<double>> &param_set, const std::string name):
    Model(init_set, param_set, name),
    _discrete_system(variables, parameters, dynamics)
{
  validate_parameters(variables, parameters, dynamics, init_set, param_set);
}