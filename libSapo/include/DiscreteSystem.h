/**
 * @file DiscreteSystem.h
 * @author Alberto Casagrande <acasagrande@units.it>
 * @brief Represent discrete systems
 * @version 0.1
 * @date 2023-02-21
 *
 * @copyright Copyright (c) 2023
 */

#ifndef _DISCRETE_SYSTEM_
#define _DISCRETE_SYSTEM_

#include "DynamicalSystem.h"

/**
 * @brief Discrete-time polynomial parametric dynamical system
 *
 * This class represents discrete-time parametric polynomial dynamical systems.
 *
 * @tparam T is the type of expression value domain
 */
template<typename T>
class DiscreteSystem : public DynamicalSystem<T>
{

  /**
   * @brief Validate the object
   */
  void validate() const;

public:
  /**
   * @brief A constructor
   *
   * @param[in] variables is the vector of the variables
   * @param[in] parameters is the vector of the parameters
   * @param[in] dynamics is the vector of the dynamic laws
   */
  DiscreteSystem(const std::vector<SymbolicAlgebra::Symbol<T>> &variables,
                 const std::vector<SymbolicAlgebra::Symbol<T>> &parameters,
                 const std::vector<SymbolicAlgebra::Expression<T>> &dynamics);

  /**
   * @brief A swap constructor
   *
   * @param[in] variables is the vector of the variables
   * @param[in] parameters is the vector of the parameters
   * @param[in] dynamics is the vector of the dynamic laws
   */
  DiscreteSystem(std::vector<SymbolicAlgebra::Symbol<T>> &&variables,
                 std::vector<SymbolicAlgebra::Symbol<T>> &&parameters,
                 std::vector<SymbolicAlgebra::Expression<T>> &&dynamics);

  /**
   * @brief An empty constructor
   */
  DiscreteSystem();

  /**
   * @brief A constructor
   *
   * @param[in] dynamics is the map that relates variables and dynamics laws
   */
  DiscreteSystem(const std::map<SymbolicAlgebra::Symbol<T>,
                                SymbolicAlgebra::Expression<T>> &dynamics);

  /**
   * @brief A constructor
   *
   * @param[in] dynamics is the map that relates variables and dynamics laws
   * @param[in] parameters is the vector of the parameters
   */
  DiscreteSystem(const std::map<SymbolicAlgebra::Symbol<T>,
                                SymbolicAlgebra::Expression<T>> &dynamics,
                 const std::vector<SymbolicAlgebra::Symbol<T>> &parameters);
  /**
   * @brief A constructor
   *
   * @param[in] dynamics is the map that relates variables and dynamic laws
   * @param[in] parameters is the vector of the parameters
   */
  DiscreteSystem(const std::map<SymbolicAlgebra::Symbol<T>,
                                SymbolicAlgebra::Expression<T>> &&dynamics,
                 std::vector<SymbolicAlgebra::Symbol<T>> &&parameters);
  /**
   * @brief A constructor
   *
   * @param[in] variables is the vector of the variables
   * @param[in] dynamics is the vector of the dynamic laws
   */
  DiscreteSystem(const std::vector<SymbolicAlgebra::Symbol<T>> &variables,
                 const std::vector<SymbolicAlgebra::Expression<T>> &dynamics);

  /**
   * @brief A swap constructor
   *
   * @param[in] variables is the vector of the variables
   * @param[in] dynamics is the vector of the dynamic laws
   */
  DiscreteSystem(std::vector<SymbolicAlgebra::Symbol<T>> &&variables,
                 std::vector<SymbolicAlgebra::Expression<T>> &&dynamics);

  /**
   * @brief A copy constructor
   *
   * @param[in] orig the original dynamic laws
   */
  DiscreteSystem(const DiscreteSystem<T> &orig);

  /**
   * @brief A swap constructor
   *
   * @param[in] orig the original dynamic laws
   */
  DiscreteSystem(DiscreteSystem<T> &&orig);

  /**
   * @brief A copy operator
   *
   * @param orig the original dynamic laws
   * @return a reference to the updated object
   */
  DiscreteSystem &operator=(const DiscreteSystem<T> &orig);

  /**
   * @brief A move operator
   *
   * @param orig the original dynamic laws
   * @return a reference to the updated object
   */
  DiscreteSystem &operator=(DiscreteSystem<T> &&orig);

  /**
   * @brief Obtain the type of the dynamical system
   *
   * @return the type of the dynamical system
   */
  DynamicType type() const;
};

template<typename T>
void DiscreteSystem<T>::validate() const
{
  using namespace SymbolicAlgebra;

  std::set<Symbol<T>> symbols;
  for (const auto &dynamic: this->_dynamics) {
    for (const auto &symbol: dynamic.get_symbols()) {
      symbols.insert(symbol);
    }
  }

  for (const auto &variable: this->_variables) {
    symbols.erase(variable);
  }

  for (const auto &parameter: this->_parameters) {
    symbols.erase(parameter);
  }

  if (!symbols.empty()) {
    SAPO_ERROR("The symbols in " << symbols << " are "
                                 << "neither variables nor parameters",
               std::domain_error);
  }
}

template<typename T>
DiscreteSystem<T>::DiscreteSystem(
    const std::vector<SymbolicAlgebra::Symbol<T>> &variables,
    const std::vector<SymbolicAlgebra::Symbol<T>> &parameters,
    const std::vector<SymbolicAlgebra::Expression<T>> &dynamics):
    DynamicalSystem<T>(variables, parameters, dynamics)
{
  validate();
}

template<typename T>
DiscreteSystem<T>::DiscreteSystem(
    std::vector<SymbolicAlgebra::Symbol<T>> &&variables,
    std::vector<SymbolicAlgebra::Symbol<T>> &&parameters,
    std::vector<SymbolicAlgebra::Expression<T>> &&dynamics):
    DynamicalSystem<T>(std::move(variables), std::move(parameters),
                       std::move(dynamics))
{
  validate();
}

template<typename T>
DiscreteSystem<T>::DiscreteSystem(): DynamicalSystem<T>()
{
}

template<typename T>
DiscreteSystem<T>::DiscreteSystem(
    const std::map<SymbolicAlgebra::Symbol<T>, SymbolicAlgebra::Expression<T>>
        &dynamics):
    DynamicalSystem<T>(dynamics)
{
  validate();
}

template<typename T>
DiscreteSystem<T>::DiscreteSystem(
    const std::map<SymbolicAlgebra::Symbol<T>, SymbolicAlgebra::Expression<T>>
        &dynamics,
    const std::vector<SymbolicAlgebra::Symbol<T>> &parameters):
    DynamicalSystem<T>(dynamics, parameters)
{
  validate();
}

template<typename T>
DiscreteSystem<T>::DiscreteSystem(
    const std::map<SymbolicAlgebra::Symbol<T>, SymbolicAlgebra::Expression<T>>
        &&dynamics,
    std::vector<SymbolicAlgebra::Symbol<T>> &&parameters):
    DynamicalSystem<T>(std::move(dynamics), std::move(parameters))
{
  validate();
}

template<typename T>
DiscreteSystem<T>::DiscreteSystem(
    const std::vector<SymbolicAlgebra::Symbol<T>> &variables,
    const std::vector<SymbolicAlgebra::Expression<T>> &dynamics):
    DynamicalSystem<T>(variables, dynamics)
{
  validate();
}

template<typename T>
DiscreteSystem<T>::DiscreteSystem(
    std::vector<SymbolicAlgebra::Symbol<T>> &&variables,
    std::vector<SymbolicAlgebra::Expression<T>> &&dynamics):
    DynamicalSystem<T>(std::move(variables), std::move(dynamics))
{
  validate();
}

template<typename T>
DiscreteSystem<T>::DiscreteSystem(const DiscreteSystem<T> &orig):
    DynamicalSystem<T>(orig)
{
}

template<typename T>
DiscreteSystem<T>::DiscreteSystem(DiscreteSystem<T> &&orig):
    DynamicalSystem<T>(std::move(orig))
{
}

template<typename T>
inline DiscreteSystem<T> &
DiscreteSystem<T>::operator=(const DiscreteSystem<T> &orig)
{
  DynamicalSystem<T>::operator=(orig);

  return *this;
}

template<typename T>
inline DiscreteSystem<T> &
DiscreteSystem<T>::operator=(DiscreteSystem<T> &&orig)
{
  DynamicalSystem<T>::operator=(std::move(orig));

  return *this;
}

template<typename T>
inline DynamicType DiscreteSystem<T>::type() const
{
  return DynamicType::DISCRETE;
}

#endif // _DISCRETE_SYSTEM_