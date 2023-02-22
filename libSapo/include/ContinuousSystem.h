/**
 * @file ContinuousSystem.h
 * @author Alberto Casagrande <acasagrande@units.it>
 * @brief Represent continuous systems
 * @version 0.1
 * @date 2023-02-22
 *
 * @copyright Copyright (c) 2023
 */

#ifndef _CONTINUOUS_SYSTEM_
#define _CONTINUOUS_SYSTEM_

#include "DynamicalSystem.h"

/**
 * @brief Continuous-time parametric polynomial dynamical system
 *
 * This class represents continuous-time parametric polynomial dynamical
 * systems.
 *
 * @tparam T is the type of expression value domain
 */
template<typename T>
class ContinuousSystem : public DynamicalSystem<T>
{
  SymbolicAlgebra::Symbol<T> _time_variable; //!< the time variable

  /**
   * @brief Validate the system
   */
  void validate() const;

public:
  /**
   * @brief An empty constructor
   */
  ContinuousSystem();

  /**
   * @brief A constructor
   *
   * @param[in] variables is the vector of the variables
   * @param[in] parameters is the vector of the parameters
   * @param[in] dynamics is the vector of the dynamic laws
   * @param[in] time_variable is the variable representing time
   *        in the dynamics
   */
  ContinuousSystem(const std::vector<SymbolicAlgebra::Symbol<T>> &variables,
                   const std::vector<SymbolicAlgebra::Symbol<T>> &parameters,
                   const std::vector<SymbolicAlgebra::Expression<T>> &dynamics,
                   const SymbolicAlgebra::Symbol<T> &time_variable);

  /**
   * @brief A swap constructor
   *
   * @param[in] variables is the vector of the variables
   * @param[in] parameters is the vector of the parameters
   * @param[in] dynamics is the vector of the dynamic laws
   * @param[in] time_variable is the variable representing time
   *        in the dynamics
   */
  ContinuousSystem(std::vector<SymbolicAlgebra::Symbol<T>> &&variables,
                   std::vector<SymbolicAlgebra::Symbol<T>> &&parameters,
                   std::vector<SymbolicAlgebra::Expression<T>> &&dynamics,
                   const SymbolicAlgebra::Symbol<T> &time_variable);

  /**
   * @brief A constructor
   *
   * @param[in] dynamics is the map that relates variables and dynamics laws
   * @param[in] time_variable is the variable representing time
   *        in the dynamics
   */
  ContinuousSystem(const std::map<SymbolicAlgebra::Symbol<T>,
                                  SymbolicAlgebra::Expression<T>> &dynamics,
                   const SymbolicAlgebra::Symbol<T> &time_variable);

  /**
   * @brief A constructor
   *
   * @param[in] dynamics is the map that relates variables and dynamics laws
   * @param[in] parameters is the vector of the parameters
   * @param[in] time_variable is the variable representing time
   *        in the dynamics
   */
  ContinuousSystem(const std::map<SymbolicAlgebra::Symbol<T>,
                                  SymbolicAlgebra::Expression<T>> &dynamics,
                   const std::vector<SymbolicAlgebra::Symbol<T>> &parameters,
                   const SymbolicAlgebra::Symbol<T> &time_variable);
  /**
   * @brief A constructor
   *
   * @param[in] dynamics is the map that relates variables and dynamic laws
   * @param[in] parameters is the vector of the parameters
   * @param[in] time_variable is the variable representing time
   *        in the dynamics
   */
  ContinuousSystem(const std::map<SymbolicAlgebra::Symbol<T>,
                                  SymbolicAlgebra::Expression<T>> &&dynamics,
                   std::vector<SymbolicAlgebra::Symbol<T>> &&parameters,
                   const SymbolicAlgebra::Symbol<T> &time_variable);
  /**
   * @brief A constructor
   *
   * @param[in] variables is the vector of the variables
   * @param[in] dynamics is the vector of the dynamic laws
   * @param[in] time_variable is the variable representing time
   *        in the dynamics
   */
  ContinuousSystem(const std::vector<SymbolicAlgebra::Symbol<T>> &variables,
                   const std::vector<SymbolicAlgebra::Expression<T>> &dynamics,
                   const SymbolicAlgebra::Symbol<T> &time_variable);

  /**
   * @brief A swap constructor
   *
   * @param[in] variables is the vector of the variables
   * @param[in] dynamics is the vector of the dynamic laws
   * @param[in] time_variable is the variable representing time
   *        in the dynamics
   */
  ContinuousSystem(std::vector<SymbolicAlgebra::Symbol<T>> &&variables,
                   std::vector<SymbolicAlgebra::Expression<T>> &&dynamics,
                   const SymbolicAlgebra::Symbol<T> &time_variable);

  /**
   * @brief A copy constructor
   *
   * @param[in] orig the original dynamic laws
   */
  ContinuousSystem(const ContinuousSystem<T> &orig);

  /**
   * @brief A swap constructor
   *
   * @param[in] orig the original dynamic laws
   */
  ContinuousSystem(ContinuousSystem<T> &&orig);

  /**
   * @brief A copy operator
   *
   * @param orig the original dynamic laws
   * @return a reference to the updated object
   */
  ContinuousSystem &operator=(const ContinuousSystem<T> &orig);

  /**
   * @brief A move operator
   *
   * @param orig the original dynamic laws
   * @return a reference to the updated object
   */
  ContinuousSystem &operator=(ContinuousSystem<T> &&orig);

  /**
   * @brief Get the time variable
   *
   * @return a reference to the time variable
   */
  const SymbolicAlgebra::Symbol<T> &time_variable() const;

  /**
   * @brief Test whether the system is time independent
   *
   * @return `true` if and only if the dynamics does not
   *    contain the time variable
   */
  bool is_time_independent() const;

  /**
   * @brief Obtain the type of the dynamical system
   *
   * @return the type of the dynamical system
   */
  DynamicType type() const override;
};

template<typename T>
void ContinuousSystem<T>::validate() const
{
  using namespace SymbolicAlgebra;

  std::set<Symbol<T>> symbols;
  for (const auto &dynamic: this->_dynamics) {
    for (const auto &symbol: dynamic.get_symbols()) {
      symbols.insert(symbol);
    }
  }

  for (size_t i = 0; i < this->_variables.size(); ++i) {
    const SymbolicAlgebra::Symbol<T> &variable = this->_variables[i];
    symbols.erase(variable);

    if (variable == _time_variable) {
      SAPO_ERROR("The time variable (i.e., \"" << _time_variable
                                               << "\") cannot have a dynamic",
                 std::domain_error);
    }
  }

  for (const auto &parameter: this->_parameters) {
    symbols.erase(parameter);
  }

  symbols.erase(_time_variable);

  if (!symbols.empty()) {
    SAPO_ERROR("The symbols in " << symbols << " have not "
                                 << "been specified",
               std::domain_error);
  }
}

template<typename T>
ContinuousSystem<T>::ContinuousSystem(
    const std::vector<SymbolicAlgebra::Symbol<T>> &variables,
    const std::vector<SymbolicAlgebra::Symbol<T>> &parameters,
    const std::vector<SymbolicAlgebra::Expression<T>> &dynamics,
    const SymbolicAlgebra::Symbol<T> &time_variable):
    DynamicalSystem<T>(variables, parameters, dynamics),
    _time_variable(time_variable)
{
  validate();
}

template<typename T>
ContinuousSystem<T>::ContinuousSystem(
    std::vector<SymbolicAlgebra::Symbol<T>> &&variables,
    std::vector<SymbolicAlgebra::Symbol<T>> &&parameters,
    std::vector<SymbolicAlgebra::Expression<T>> &&dynamics,
    const SymbolicAlgebra::Symbol<T> &time_variable):
    DynamicalSystem<T>(std::move(variables), std::move(parameters),
                       std::move(dynamics)),
    _time_variable(time_variable)
{
  validate();
}

template<typename T>
ContinuousSystem<T>::ContinuousSystem(): DynamicalSystem<T>()
{
}

template<typename T>
ContinuousSystem<T>::ContinuousSystem(
    const std::map<SymbolicAlgebra::Symbol<T>, SymbolicAlgebra::Expression<T>>
        &dynamics,
    const SymbolicAlgebra::Symbol<T> &time_variable):
    DynamicalSystem<T>(dynamics),
    _time_variable(time_variable)
{
  validate();
}

template<typename T>
ContinuousSystem<T>::ContinuousSystem(
    const std::map<SymbolicAlgebra::Symbol<T>, SymbolicAlgebra::Expression<T>>
        &dynamics,
    const std::vector<SymbolicAlgebra::Symbol<T>> &parameters,
    const SymbolicAlgebra::Symbol<T> &time_variable):
    DynamicalSystem<T>(dynamics, parameters),
    _time_variable(time_variable)
{
  validate();
}

template<typename T>
ContinuousSystem<T>::ContinuousSystem(
    const std::map<SymbolicAlgebra::Symbol<T>, SymbolicAlgebra::Expression<T>>
        &&dynamics,
    std::vector<SymbolicAlgebra::Symbol<T>> &&parameters,
    const SymbolicAlgebra::Symbol<T> &time_variable):
    DynamicalSystem<T>(std::move(dynamics), std::move(parameters)),
    _time_variable(time_variable)
{
  validate();
}

template<typename T>
ContinuousSystem<T>::ContinuousSystem(
    const std::vector<SymbolicAlgebra::Symbol<T>> &variables,
    const std::vector<SymbolicAlgebra::Expression<T>> &dynamics,
    const SymbolicAlgebra::Symbol<T> &time_variable):
    DynamicalSystem<T>(variables, dynamics),
    _time_variable(time_variable)
{
  validate();
}

template<typename T>
ContinuousSystem<T>::ContinuousSystem(
    std::vector<SymbolicAlgebra::Symbol<T>> &&variables,
    std::vector<SymbolicAlgebra::Expression<T>> &&dynamics,
    const SymbolicAlgebra::Symbol<T> &time_variable):
    DynamicalSystem<T>(std::move(variables), std::move(dynamics)),
    _time_variable(time_variable)
{
  validate();
}

template<typename T>
ContinuousSystem<T>::ContinuousSystem(const ContinuousSystem<T> &orig):
    DynamicalSystem<T>(orig), _time_variable(orig._time_variable)
{
  validate();
}

template<typename T>
ContinuousSystem<T>::ContinuousSystem(ContinuousSystem<T> &&orig):
    DynamicalSystem<T>(std::move(orig)), _time_variable(orig.time_variable)
{
  validate();
}

template<typename T>
inline ContinuousSystem<T> &
ContinuousSystem<T>::operator=(const ContinuousSystem<T> &orig)
{
  DynamicalSystem<T>::operator=(orig);

  _time_variable = orig._time_variable;

  return *this;
}

template<typename T>
inline ContinuousSystem<T> &
ContinuousSystem<T>::operator=(ContinuousSystem<T> &&orig)
{
  DynamicalSystem<T>::operator=(std::move(orig));

  _time_variable = orig._time_variable;

  return *this;
}

template<typename T>
inline const SymbolicAlgebra::Symbol<T> &
ContinuousSystem<T>::time_variable() const
{
  return _time_variable;
}

template<typename T>
bool ContinuousSystem<T>::is_time_independent() const
{
  for (const auto &dynamic: this->_dynamics) {
    for (const auto &symbol: dynamic.get_symbols()) {
      if (symbol == _time_variable) {
        return false;
      }
    }
  }

  return true;
}

template<typename T>
inline DynamicType ContinuousSystem<T>::type() const
{
  return DynamicType::CONTINUOUS;
}

#endif // _CONTINUOUS_SYSTEM_