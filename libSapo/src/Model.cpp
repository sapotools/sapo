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

Model::Model(const std::vector<SymbolicAlgebra::Symbol<>> &variables,
             const std::vector<SymbolicAlgebra::Expression<>> &dynamics,
             const Bundle &init_set, const std::string name):
    Model(variables, {}, dynamics, init_set, SetsUnion<Polytope>(), name)
{
}

Model::Model(const std::vector<SymbolicAlgebra::Symbol<>> &variables,
             const std::vector<SymbolicAlgebra::Symbol<>> &parameters,
             const std::vector<SymbolicAlgebra::Expression<>> &dynamics,
             const Bundle &init_set, const SetsUnion<Polytope> &parameter_set,
             const std::string name):
    Model(DynamicalSystem<double>(variables, parameters, dynamics), init_set,
          parameter_set, name)
{
}

Model::Model(const DynamicalSystem<double> &dynamical_system,
             const Bundle &init_set, const SetsUnion<Polytope> &parameter_set,
             const std::string name):
    _dynamical_system(dynamical_system),
    _init_set(std::make_shared<Bundle>(init_set)), _param_set(parameter_set),
    _spec(nullptr), _assumptions(), _name(name)
{
  using namespace SymbolicAlgebra;
  std::set<Symbol<>> dyn_symbols;

  for (const auto &dyn: _dynamical_system.dynamics()) {
    auto symbols = dyn.get_symbols();
    dyn_symbols.insert(std::begin(symbols), std::end(symbols));
  }

  for (const auto &var: _dynamical_system.variables()) {
    dyn_symbols.erase(var);
  }

  for (const auto &param: _dynamical_system.parameters()) {
    dyn_symbols.erase(param);
  }

  if (!dyn_symbols.empty()) {
    std::ostringstream oss;
    oss << "The dynamic law symbols " << dyn_symbols
        << "are neither variables nor parameters";
    throw std::domain_error(oss.str());
  }

  if (init_set.dim() != _dynamical_system.dim()) {
    std::ostringstream oss;
    oss << "The space dimension of the initial set (" << init_set.dim()
        << ") differs from the "
        << "variable number (" << _dynamical_system.dim()
        << "). They must be the same.";
    throw std::domain_error(oss.str());
  }

  if (parameter_set.dim() != _dynamical_system.parameters().size()) {
    std::ostringstream oss;
    oss << "The space dimension of the parameter set (" << parameter_set.dim()
        << ") differs from the parameter number ("
        << _dynamical_system.parameters().size()
        << "). They must be the same.";
    throw std::domain_error(oss.str());
  }
}

Model &Model::set_specification(std::shared_ptr<STL::STL> specification)
{
  if (specification != nullptr) {
    auto spec_vars = specification->get_variables();

    for (const auto &var: _dynamical_system.variables()) {
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
  if (assumptions.size() > 0 && assumptions.dim() != _dynamical_system.dim()) {
    throw std::domain_error("The number of dimensions of the "
                            "assumptions space differs from "
                            "that of the model variables.");
  }

  _assumptions = LinearSystem(assumptions);

  return *this;
}

Model &Model::set_invariant(const LinearSystem &invariant)
{
  if (invariant.size() > 0 && invariant.dim() != _dynamical_system.dim()) {
    throw std::domain_error("The number of dimensions of the "
                            "invariant space differs from "
                            "that of the model variables.");
  }

  _invariant = LinearSystem(invariant);

  return *this;
}
