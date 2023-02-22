/**
 * @file DynamicalSystem.h
 * @author Alberto Casagrande <acasagrande@units.it>
 * @brief Represent dynamic system laws
 * @version 0.1
 * @date 2022-05-23
 *
 * @copyright Copyright (c) 2022
 */

#ifndef _DYNAMICAL_SYSTEM_
#define _DYNAMICAL_SYSTEM_

#include <map>
#include <vector>
#include <sstream>

#include "ErrorHandling.h"
#include "SymbolicAlgebra.h"
#include "LinearAlgebraIO.h"

/**
 * @brief Dynamic type
 */
enum class DynamicType { DISCRETE, CONTINUOUS, UNDEFINED };

/**
 * @brief This class represents parametric polynomial dynamical systems
 *
 * This class represents parametric polynomial dynamical systems.
 * Each variable of the system is associated to an expression that
 * describes the variable dynamic law.
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
   * @brief Conclude initialization
   */
  void build_variable_index();

  /**
   * @brief Validate the object
   */
  void validate() const;

  /**
   * @brief A constructor
   *
   * @param[in] variables is the vector of the variables
   * @param[in] parameters is the vector of the parameters
   * @param[in] dynamics is the vector of the dynamic laws
   */
  DynamicalSystem(const std::vector<SymbolicAlgebra::Symbol<T>> &variables,
                  const std::vector<SymbolicAlgebra::Symbol<T>> &parameters,
                  const std::vector<SymbolicAlgebra::Expression<T>> &dynamics);

  /**
   * @brief A swap constructor
   *
   * @param[in] variables is the vector of the variables
   * @param[in] parameters is the vector of the parameters
   * @param[in] dynamics is the vector of the dynamic laws
   */
  DynamicalSystem(std::vector<SymbolicAlgebra::Symbol<T>> &&variables,
                  std::vector<SymbolicAlgebra::Symbol<T>> &&parameters,
                  std::vector<SymbolicAlgebra::Expression<T>> &&dynamics);

  /**
   * @brief A constructor
   *
   * @param[in] dynamics is the map that relates variables and dynamics laws
   */
  DynamicalSystem(const std::map<SymbolicAlgebra::Symbol<T>,
                                 SymbolicAlgebra::Expression<T>> &dynamics);

  /**
   * @brief A constructor
   *
   * @param[in] dynamics is the map that relates variables and dynamics laws
   * @param[in] parameters is the vector of the parameters
   */
  DynamicalSystem(const std::map<SymbolicAlgebra::Symbol<T>,
                                 SymbolicAlgebra::Expression<T>> &dynamics,
                  const std::vector<SymbolicAlgebra::Symbol<T>> &parameters);

  /**
   * @brief A constructor
   *
   * @param[in] dynamics is the map that relates variables and dynamic laws
   * @param[in] parameters is the vector of the parameters
   */
  DynamicalSystem(const std::map<SymbolicAlgebra::Symbol<T>,
                                 SymbolicAlgebra::Expression<T>> &&dynamics,
                  std::vector<SymbolicAlgebra::Symbol<T>> &&parameters);

  /**
   * @brief A constructor
   *
   * @param[in] variables is the vector of the variables
   * @param[in] dynamics is the vector of the dynamic laws
   */
  DynamicalSystem(const std::vector<SymbolicAlgebra::Symbol<T>> &variables,
                  const std::vector<SymbolicAlgebra::Expression<T>> &dynamics);

  /**
   * @brief A swap constructor
   *
   * @param[in] variables is the vector of the variables
   * @param[in] dynamics is the vector of the dynamic laws
   */
  DynamicalSystem(std::vector<SymbolicAlgebra::Symbol<T>> &&variables,
                  std::vector<SymbolicAlgebra::Expression<T>> &&dynamics);

  /**
   * @brief A copy constructor
   *
   * @param[in] orig the original dynamic laws
   */
  DynamicalSystem(const DynamicalSystem<T> &orig);

  /**
   * @brief A swap constructor
   *
   * @param[in] orig the original dynamic laws
   */
  DynamicalSystem(DynamicalSystem<T> &&orig);

public:
  /**
   * @brief An empty constructor
   */
  DynamicalSystem();

  /**
   * @brief Number of dynamic laws
   *
   * @return return the number of represented dynamic laws
   */
  size_t dim() const;

  /**
   * @brief A copy operator
   *
   * @param orig the original dynamic laws
   * @return a reference to the updated object
   */
  DynamicalSystem &operator=(const DynamicalSystem<T> &orig);

  /**
   * @brief A move operator
   *
   * @param orig the original dynamic laws
   * @return a reference to the updated object
   */
  DynamicalSystem &operator=(DynamicalSystem<T> &&orig);

  /**
   * @brief Get the system variables
   *
   * @return a reference to the vector of the system variables
   */
  const std::vector<SymbolicAlgebra::Symbol<T>> &variables() const;

  /**
   * @brief Get the system parameters
   *
   * @return a reference to the vector of the system parameters
   */
  const std::vector<SymbolicAlgebra::Symbol<T>> &parameters() const;

  /**
   * @brief Get the system dynamic laws
   *
   * @return a reference to the system dynamic law vector
   */
  const std::vector<SymbolicAlgebra::Expression<T>> &dynamics() const;

  /**
   * @brief Get the i-th variable
   *
   * @return a reference to the i-th variable
   */
  const SymbolicAlgebra::Symbol<T> &variable(const size_t i) const;

  /**
   * @brief Get the i-th parameter
   *
   * @return a reference to the i-th parameter
   */
  const SymbolicAlgebra::Symbol<T> &parameter(const size_t i) const;

  /**
   * @brief Get the i-th dynamic law
   *
   * @return a reference to the i-th dynamic law
   */
  const SymbolicAlgebra::Expression<T> &dynamic(const size_t i) const;

  /**
   * @brief Get the dynamic law associated to a variable
   *
   * @return a reference to the dynamic law associated to a variable
   */
  const SymbolicAlgebra::Expression<T> &
  operator[](const SymbolicAlgebra::Symbol<T> &variable) const;

  /**
   * @brief Obtain the type of the dynamical system
   *
   * @return the type of the dynamical system
   */
  virtual DynamicType type() const;
};

template<typename T>
void DynamicalSystem<T>::build_variable_index()
{
  using namespace SymbolicAlgebra;

  for (const auto &variable: _variables) {
    if (!_var_index.emplace(variable.get_id(), _var_index.size()).second) {
      SAPO_ERROR("The variable vector contains \""
                     << Symbol<T>::get_symbol_name(variable.get_id())
                     << "\" twice",
                 std::domain_error);
    }
  }

  std::set<SymbolIdType> param_index;
  for (const auto &parameter: _parameters) {
    if (_var_index.find(parameter.get_id()) != std::end(_var_index)) {
      SAPO_ERROR("Symbol \"" << Symbol<T>::get_symbol_name(parameter.get_id())
                             << "\" must be either a variable or a parameter",
                 std::domain_error);
    }
    if (!param_index.emplace(parameter.get_id()).second) {
      SAPO_ERROR("The parameter vector contains \""
                     << Symbol<T>::get_symbol_name(parameter.get_id())
                     << "\" twice",
                 std::domain_error);
    }
  }
}

template<typename T>
void DynamicalSystem<T>::validate() const
{
  if (_variables.size() != _dynamics.size()) {
    SAPO_ERROR("the vectors of the dynamics laws and that "
               "of the variables must have the same size",
               std::domain_error);
  }
}

template<typename T>
DynamicalSystem<T>::DynamicalSystem(
    const std::vector<SymbolicAlgebra::Symbol<T>> &variables,
    const std::vector<SymbolicAlgebra::Symbol<T>> &parameters,
    const std::vector<SymbolicAlgebra::Expression<T>> &dynamics):
    _variables(variables),
    _parameters(parameters), _dynamics(dynamics), _var_index()
{
  build_variable_index();
  validate();
}

template<typename T>
DynamicalSystem<T>::DynamicalSystem(
    std::vector<SymbolicAlgebra::Symbol<T>> &&variables,
    std::vector<SymbolicAlgebra::Symbol<T>> &&parameters,
    std::vector<SymbolicAlgebra::Expression<T>> &&dynamics):
    _variables(std::move(variables)),
    _parameters(std::move(parameters)), _dynamics(std::move(dynamics)),
    _var_index()
{
  build_variable_index();
  validate();
}

template<typename T>
DynamicalSystem<T>::DynamicalSystem(): DynamicalSystem({}, {}, {})
{
}

template<typename T>
DynamicalSystem<T>::DynamicalSystem(
    const std::map<SymbolicAlgebra::Symbol<T>, SymbolicAlgebra::Expression<T>>
        &dynamics):
    _variables(),
    _parameters(), _dynamics(), _var_index()
{
  _variables.reserve(dynamics.size());
  _dynamics.reserve(dynamics.size());

  for (auto d_it = std::begin(dynamics); d_it != std::end(dynamics); ++d_it) {
    _variables.push_back(d_it->first);
    _dynamics.push_back(d_it->second);
  }

  build_variable_index();
  validate();
}

template<typename T>
DynamicalSystem<T>::DynamicalSystem(
    const std::map<SymbolicAlgebra::Symbol<T>, SymbolicAlgebra::Expression<T>>
        &dynamics,
    const std::vector<SymbolicAlgebra::Symbol<T>> &parameters):
    _variables(),
    _parameters(parameters), _dynamics(), _var_index()
{
  _variables.reserve(dynamics.size());
  _dynamics.reserve(dynamics.size());

  for (const auto &dynamic: dynamics) {
    _variables.push_back(dynamic.first);
    _dynamics.push_back(dynamic.second);
  }

  build_variable_index();
  validate();
}

template<typename T>
DynamicalSystem<T>::DynamicalSystem(
    const std::map<SymbolicAlgebra::Symbol<T>, SymbolicAlgebra::Expression<T>>
        &&dynamics,
    std::vector<SymbolicAlgebra::Symbol<T>> &&parameters):
    _variables(),
    _parameters(std::move(parameters)), _dynamics(), _var_index()
{
  _variables.reserve(dynamics.size());
  _dynamics.reserve(dynamics.size());

  for (auto &dynamic: dynamics) {
    _variables.push_back(std::move(dynamic.first));
    _dynamics.push_back(std::move(dynamic.second));
  }

  build_variable_index();
  validate();
}

template<typename T>
DynamicalSystem<T>::DynamicalSystem(
    const std::vector<SymbolicAlgebra::Symbol<T>> &variables,
    const std::vector<SymbolicAlgebra::Expression<T>> &dynamics):
    DynamicalSystem(variables, {}, dynamics)
{
}

template<typename T>
DynamicalSystem<T>::DynamicalSystem(
    std::vector<SymbolicAlgebra::Symbol<T>> &&variables,
    std::vector<SymbolicAlgebra::Expression<T>> &&dynamics):
    DynamicalSystem(std::move(variables), {}, std::move(dynamics))
{
}

template<typename T>
DynamicalSystem<T>::DynamicalSystem(const DynamicalSystem<T> &orig):
    _variables(orig._variables), _parameters(orig._parameters),
    _dynamics(orig._dynamics), _var_index(orig._var_index)
{
}

template<typename T>
DynamicalSystem<T>::DynamicalSystem(DynamicalSystem<T> &&orig)
{
  std::swap(orig._variables, _variables);
  std::swap(orig._parameters, _parameters);
  std::swap(orig._dynamics, _dynamics);
  std::swap(orig._var_index, _var_index);
}

template<typename T>
size_t DynamicalSystem<T>::dim() const
{
  return _variables.size();
}

template<typename T>
DynamicalSystem<T> &
DynamicalSystem<T>::operator=(const DynamicalSystem<T> &orig)
{
  _variables = orig._variables;
  _parameters = orig._parameters;
  _dynamics = orig._dynamics;
  _var_index = orig._var_index;

  return *this;
}

template<typename T>
DynamicalSystem<T> &DynamicalSystem<T>::operator=(DynamicalSystem<T> &&orig)
{
  std::swap(orig._variables, _variables);
  std::swap(orig._parameters, _parameters);
  std::swap(orig._dynamics, _dynamics);
  std::swap(orig._var_index, _var_index);

  return *this;
}

template<typename T>
inline const std::vector<SymbolicAlgebra::Symbol<T>> &
DynamicalSystem<T>::variables() const
{
  return _variables;
}

template<typename T>
inline const std::vector<SymbolicAlgebra::Symbol<T>> &
DynamicalSystem<T>::parameters() const
{
  return _parameters;
}

template<typename T>
inline const std::vector<SymbolicAlgebra::Expression<T>> &
DynamicalSystem<T>::dynamics() const
{
  return _dynamics;
}

template<typename T>
inline const SymbolicAlgebra::Symbol<T> &
DynamicalSystem<T>::variable(const size_t i) const
{
  return _variables[i];
}

template<typename T>
inline const SymbolicAlgebra::Symbol<T> &
DynamicalSystem<T>::parameter(const size_t i) const
{
  return _parameters[i];
}

template<typename T>
inline const SymbolicAlgebra::Expression<T> &
DynamicalSystem<T>::dynamic(const size_t i) const
{
  return _dynamics[i];
}

template<typename T>
inline const SymbolicAlgebra::Expression<T> &DynamicalSystem<T>::operator[](
    const SymbolicAlgebra::Symbol<T> &variable) const
{
  return _dynamics[_var_index.at[variable.get_id()]];
}

template<typename T>
inline DynamicType DynamicalSystem<T>::type() const
{
  return DynamicType::UNDEFINED;
}

#endif // _DYNAMICAL_SYSTEM_