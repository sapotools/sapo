/**
 * @file DynamicalSystem.h
 * @author Alberto Casagrande <acasagrande@units.it>
 * @brief Represent dynamic system laws
 * @version 0.1
 * @date 2022-05-23
 *
 * @copyright Copyright (c) 2022
 */

#ifndef DYNAMICS_H_
#define DYNAMICS_H_

#include <map>
#include <vector>
#include <sstream>

#include "SymbolicAlgebra.h"

#include "LinearAlgebraIO.h"

/**
 * @brief This class represents dynamical systems
 *
 * This class represents dynamical systems. Each
 * variable of the system is updated according with
 * the corresponding expression.
 *
 * @tparam T is the type of expression value domain
 */
template<typename T>
class DynamicalSystem
{
protected:
  typedef typename SymbolicAlgebra::Symbol<T>::SymbolIdType
      SymbolIdType; //!< The type of symbol identificators

  std::vector<SymbolicAlgebra::Symbol<T>>
      _variables; //!< the vector of the variables
  std::vector<SymbolicAlgebra::Symbol<T>>
      _parameters; //!< the vector of the parameters
  std::vector<SymbolicAlgebra::Expression<T>>
      _dynamics; //!< the vector of the dynamic laws
  std::map<SymbolIdType, unsigned int>
      _var_index; //!< a map from the variable identifiers to the corresponding
                  //!< indices

  /**
   * @brief Validate constructor parameters and initialize objects
   */
  void validate_and_initialize()
  {
    using namespace SymbolicAlgebra;

    if (_variables.size() != _dynamics.size()) {
      throw std::domain_error("The vectors of the dynamics laws and that "
                              "of the variables must have the same size");
    }

    std::set<Symbol<T>> symbols;

    for (auto d_it = std::begin(_dynamics); d_it != std::end(_dynamics);
         ++d_it) {

      std::set<Symbol<T>> d_symbols = d_it->get_symbols();

      for (auto ds_it = std::begin(d_symbols); ds_it != std::end(d_symbols);
           ++ds_it) {
        symbols.insert(*ds_it);
      }
    }

    for (auto v_it = std::begin(_variables); v_it != std::end(_variables);
         ++v_it) {
      symbols.erase(*v_it);
      if (!_var_index.emplace(v_it->get_id(), _var_index.size()).second) {
        throw std::domain_error("The vector or the variables must not "
                                "contain duplicates");
      }
    }

    std::set<SymbolIdType> param_index;
    for (auto p_it = std::begin(_parameters); p_it != std::end(_parameters);
         ++p_it) {
      symbols.erase(*p_it);
      if (_var_index.find(p_it->get_id()) != std::end(_var_index)) {
        std::ostringstream oss;

        oss << "Symbol \"" << Symbol<>::get_symbol_name(p_it->get_id())
            << "\" is reported as both a variable and a parameter";

        throw std::domain_error(oss.str());
      }
      if (!param_index.emplace(p_it->get_id()).second) {
        throw std::domain_error("The vector of the parameters must not "
                                "contain duplicates");
      }
    }

    if (!symbols.empty()) {
      std::ostringstream oss;

      oss << "The symbols in " << symbols << " are neither variables "
          << "nor parameters";

      throw std::domain_error(oss.str());
    }
  }

public:
  /**
   * @brief A constructor
   *
   * @param[in] variables is the vector of the variables
   * @param[in] parameters is the vector of the parameters
   * @param[in] dynamics is the vector of the dynamic laws
   */
  DynamicalSystem(const std::vector<SymbolicAlgebra::Symbol<T>> &variables,
                  const std::vector<SymbolicAlgebra::Symbol<T>> &parameters,
                  const std::vector<SymbolicAlgebra::Expression<T>> &dynamics):
      _variables(variables),
      _parameters(parameters), _dynamics(dynamics), _var_index()
  {
    validate_and_initialize();
  }

  /**
   * @brief A swap constructor
   *
   * @param[in] variables is the vector of the variables
   * @param[in] parameters is the vector of the parameters
   * @param[in] dynamics is the vector of the dynamic laws
   */
  DynamicalSystem(std::vector<SymbolicAlgebra::Symbol<T>> &&variables,
                  std::vector<SymbolicAlgebra::Symbol<T>> &&parameters,
                  std::vector<SymbolicAlgebra::Expression<T>> &&dynamics):
      _variables(std::move(variables)),
      _parameters(std::move(parameters)), _dynamics(std::move(dynamics)),
      _var_index()
  {
    validate_and_initialize();
  }

  /**
   * @brief An empty constructor
   */
  DynamicalSystem(): DynamicalSystem({}, {}, {}) {}

  /**
   * @brief A constructor
   *
   * @param[in] dynamics is the map that relates variables and dynamics laws
   */
  DynamicalSystem(const std::map<SymbolicAlgebra::Symbol<T>,
                                 SymbolicAlgebra::Expression<T>> &dynamics):
      _variables(),
      _parameters(), _dynamics(), _var_index()
  {
    _variables.reserve(dynamics.size());
    _dynamics.reserve(dynamics.size());

    for (auto d_it = std::begin(dynamics); d_it != std::end(dynamics);
         ++d_it) {
      _variables.push_back(d_it->first);
      _dynamics.push_back(d_it->second);
    }

    validate_and_initialize();
  }

  /**
   * @brief A constructor
   *
   * @param[in] dynamics is the map that relates variables and dynamics laws
   * @param[in] parameters is the vector of the parameters
   */
  DynamicalSystem(const std::map<SymbolicAlgebra::Symbol<T>,
                                 SymbolicAlgebra::Expression<T>> &dynamics,
                  const std::vector<SymbolicAlgebra::Symbol<T>> &parameters):
      _variables(),
      _parameters(parameters), _dynamics(), _var_index()
  {
    _variables.reserve(dynamics.size());
    _dynamics.reserve(dynamics.size());

    for (auto d_it = std::begin(dynamics); d_it != std::end(dynamics);
         ++d_it) {
      _variables.push_back(d_it->first);
      _dynamics.push_back(d_it->second);
    }

    validate_and_initialize();
  }

  /**
   * @brief A constructor
   *
   * @param[in] dynamics is the map that relates variables and dynamic laws
   * @param[in] parameters is the vector of the parameters
   */
  DynamicalSystem(const std::map<SymbolicAlgebra::Symbol<T>,
                                 SymbolicAlgebra::Expression<T>> &&dynamics,
                  std::vector<SymbolicAlgebra::Symbol<T>> &&parameters):
      _variables(),
      _parameters(std::move(parameters)), _dynamics(), _var_index()
  {
    _variables.reserve(dynamics.size());
    _dynamics.reserve(dynamics.size());

    for (auto d_it = std::begin(dynamics); d_it != std::end(dynamics);
         ++d_it) {
      _variables.push_back(std::move(d_it->first));
      _dynamics.push_back(std::move(d_it->second));
    }

    validate_and_initialize();
  }

  /**
   * @brief A constructor
   *
   * @param[in] variables is the vector of the variables
   * @param[in] dynamics is the vector of the dynamic laws
   */
  DynamicalSystem(const std::vector<SymbolicAlgebra::Symbol<T>> &variables,
                  const std::vector<SymbolicAlgebra::Expression<T>> &dynamics):
      DynamicalSystem(variables, {}, dynamics)
  {
  }

  /**
   * @brief A swap constructor
   *
   * @param[in] variables is the vector of the variables
   * @param[in] dynamics is the vector of the dynamic laws
   */
  DynamicalSystem(std::vector<SymbolicAlgebra::Symbol<T>> &&variables,
                  std::vector<SymbolicAlgebra::Expression<T>> &&dynamics):
      DynamicalSystem(std::move(variables), {}, std::move(dynamics))
  {
  }

  /**
   * @brief A copy constructor
   *
   * @param[in] orig the original dynamic laws
   */
  DynamicalSystem(const DynamicalSystem<T> &orig):
      _variables(orig._variables), _parameters(orig._parameters),
      _dynamics(orig._dynamics), _var_index(orig._var_index)
  {
  }

  /**
   * @brief A swap constructor
   *
   * @param[in] orig the original dynamic laws
   */
  DynamicalSystem(DynamicalSystem<T> &&orig)
  {
    std::swap(orig._variables, _variables);
    std::swap(orig._parameters, _parameters);
    std::swap(orig._dynamics, _dynamics);
    std::swap(orig._var_index, _var_index);
  }

  /**
   * @brief Number of dynamic laws
   *
   * @return return the number of represented dynamic laws
   */
  size_t dim() const
  {
    return _variables.size();
  }

  /**
   * @brief A copy operator
   *
   * @param orig the original dynamic laws
   * @return a reference to the updated object
   */
  DynamicalSystem &operator=(const DynamicalSystem<T> &orig)
  {
    _variables = orig._variables;
    _parameters = orig._parameters;
    _dynamics = orig._dynamics;
    _var_index = orig._var_index;

    return *this;
  }

  /**
   * @brief A copy operator
   *
   * @param orig the original dynamic laws
   * @return a reference to the updated object
   */
  DynamicalSystem &operator=(DynamicalSystem<T> &&orig)
  {
    std::swap(orig._variables, _variables);
    std::swap(orig._parameters, _parameters);
    std::swap(orig._dynamics, _dynamics);
    std::swap(orig._var_index, _var_index);

    return *this;
  }

  /**
   * @brief Replace the dynamic law of a variable
   *
   * @param variable is the variable changed by the dynamic
   * @param dynamic is the new dynamic law for `variable`
   * @return a reference to the updated object
   */
  DynamicalSystem &replace(const SymbolicAlgebra::Symbol<T> &variable,
                           const SymbolicAlgebra::Expression<T> &dynamic)
  {
    _dynamics[_var_index.at[variable.get_id()]] = dynamic;

    return *this;
  }

  /**
   * @brief Replace the dynamic law of a variable
   *
   * @param variable is the variable changed by the dynamic
   * @param dynamic is the new dynamic law for `variable`
   * @return a reference to the updated object
   */
  DynamicalSystem &replace(const SymbolicAlgebra::Symbol<T> &variable,
                           SymbolicAlgebra::Expression<T> &&dynamic)
  {

    _dynamics[_var_index.at[variable.get_id()]] = std::move(dynamic);

    return *this;
  }

  /**
   * @brief Get the system variables
   *
   * @return a reference to the vector of the system variables
   */
  inline const std::vector<SymbolicAlgebra::Symbol<T>> &variables() const
  {
    return _variables;
  }

  /**
   * @brief Get the system parameters
   *
   * @return a reference to the vector of the system parameters
   */
  inline const std::vector<SymbolicAlgebra::Symbol<T>> &parameters() const
  {
    return _parameters;
  }

  /**
   * @brief Get the system dynamic laws
   *
   * @return a reference to the system dynamic law vector
   */
  inline const std::vector<SymbolicAlgebra::Expression<T>> &dynamics() const
  {
    return _dynamics;
  }

  /**
   * @brief Get the i-th variable
   *
   * @return a reference to the i-th variable
   */
  inline const SymbolicAlgebra::Symbol<T> &variable(const unsigned int i) const
  {
    return _variables[i];
  }

  /**
   * @brief Get the i-th parameter
   *
   * @return a reference to the i-th parameter
   */
  inline const SymbolicAlgebra::Symbol<T> &
  parameter(const unsigned int i) const
  {
    return _parameters[i];
  }

  /**
   * @brief Get the i-th dynamic law
   *
   * @return a reference to the i-th dynamic law
   */
  inline const SymbolicAlgebra::Expression<T> &
  dynamic(const unsigned int i) const
  {
    return _dynamics[i];
  }

  /**
   * @brief Get the dynamic law associated to a variable
   *
   * @return a reference to the dynamic law associated to a variable
   */
  inline const SymbolicAlgebra::Expression<T> &
  operator[](const SymbolicAlgebra::Symbol<T> &variable) const
  {
    return _dynamics[_var_index.at[variable.get_id()]];
  }
};

/**
 * @brief Integrate an ODE by using the 4th degree Runge-Kutta method
 *
 * @tparam T is the constant type in the dynamical system
 * @tparam DELTA_TYPE is the type of the delta step
 * @param variables is the variable vector
 * @param ODE is the ODE vector
 * @param time is the time variable
 * @param timestep is the integration time step
 * @return the integration of `ODE` according to the 4th degree
 *         Runge-Kutta method
 */
template<typename T, typename DELTA_TYPE>
std::vector<SymbolicAlgebra::Expression<T>>
runge_kutta4(const std::vector<SymbolicAlgebra::Symbol<T>> &variables,
             const std::vector<SymbolicAlgebra::Expression<T>> &ODE,
             const SymbolicAlgebra::Symbol<T> *time,
             const DELTA_TYPE &timestep)
{
  using namespace SymbolicAlgebra;
  using namespace LinearAlgebra;

  const auto &k1 = ODE;

  // get k2
  auto k2 = ODE;
  {
    Expression<>::replacement_type repl;
    for (size_t i = 0; i < variables.size(); ++i) {
      repl[variables[i]] = variables[i] + k1[i] * timestep / 2;
    }
    if (time != nullptr) {
      repl[*time] = *time + timestep / 2;
    }

    for (auto &k: k2) {
      k.replace(repl);
    }
  }

  // get k3
  auto k3 = ODE;
  {
    Expression<>::replacement_type repl;
    for (size_t i = 0; i < variables.size(); ++i) {
      repl[variables[i]] = variables[i] + k2[i] * timestep / 2;
    }
    if (time != nullptr) {
      repl[*time] = *time + timestep / 2;
    }

    for (auto &k: k3) {
      k.replace(repl);
    }
  }

  // get k4
  auto k4 = ODE;
  {
    Expression<>::replacement_type repl;
    for (size_t i = 0; i < variables.size(); ++i) {
      repl[variables[i]] = variables[i] + k3[i] * timestep;
    }
    if (time != nullptr) {
      repl[*time] = *time + timestep;
    }

    for (auto &k: k4) {
      k.replace(repl);
    }
  }

  std::vector<Expression<T>> rk_dyn(variables.size());

  for (size_t i = 0; i < variables.size(); ++i) {
    rk_dyn[i] = variables[i]
                + (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) * timestep / 6;
  }

  return rk_dyn;
}

/**
 * @brief Integrate an ODE by using the 4th degree Runge-Kutta method
 *
 * @tparam T is the constant type in the dynamical system
 * @tparam DELTA_TYPE is the type of the delta step
 * @param variables is the variable vector
 * @param ODE is the ODE vector
 * @param time is the time variable
 * @param timestep is the integration time step
 * @return the integration of `ODE` according to the 4th degree
 *         Runge-Kutta method
 */
template<typename T, typename DELTA_TYPE>
inline std::vector<SymbolicAlgebra::Expression<T>>
runge_kutta4(const std::vector<SymbolicAlgebra::Symbol<T>> &variables,
             const std::vector<SymbolicAlgebra::Expression<T>> &ODE,
             const SymbolicAlgebra::Symbol<T> &time,
             const DELTA_TYPE &timestep)
{
  return runge_kutta4<T, DELTA_TYPE>(variables, ODE, &time, timestep);
}

/**
 * @brief Integrate an time-invariant ODE by using the 4th degree Runge-Kutta
 * method
 *
 * @tparam T is the constant type in the dynamical system
 * @tparam DELTA_TYPE is the type of the delta step
 * @param variables is the variable vector
 * @param ODE is the ODE vector
 * @param timestep is the integration time step
 * @return the integration of `ODE` according to the 4th degree
 *         Runge-Kutta method
 */
template<typename T, typename DELTA_TYPE>
inline std::vector<SymbolicAlgebra::Expression<T>>
runge_kutta4(const std::vector<SymbolicAlgebra::Symbol<T>> &variables,
             const std::vector<SymbolicAlgebra::Expression<T>> &ODE,
             const DELTA_TYPE &timestep)
{
  return runge_kutta4<T, DELTA_TYPE>(variables, ODE, nullptr, timestep);
}

#endif // DYNAMICS_H_