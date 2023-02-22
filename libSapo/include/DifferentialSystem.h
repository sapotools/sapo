/**
 * @file DifferentialSystem.h
 * @author Alberto Casagrande <acasagrande@units.it>
 * @brief Represent differential systems
 * @version 0.1
 * @date 2023-02-21
 *
 * @copyright Copyright (c) 2023
 */

#ifndef _DIFFERENTIAL_SYSTEM_H_
#define _DIFFERENTIAL_SYSTEM_H_

#include "ContinuousSystem.h"
#include "SymbolicAlgebra.h"

template<typename T>
class ODE : public ContinuousSystem<T>
{
public:
  /**
   * @brief An empty constructor
   */
  ODE();

  /**
   * @brief A constructor
   *
   * @param[in] variables is the vector of the variables
   * @param[in] parameters is the vector of the parameters
   * @param[in] dynamics is the vector of the ODE
   */
  ODE(const std::vector<SymbolicAlgebra::Symbol<T>> &variables,
      const std::vector<SymbolicAlgebra::Symbol<T>> &parameters,
      const std::vector<SymbolicAlgebra::Expression<T>> &dynamics,
      const SymbolicAlgebra::Symbol<T> &time_variable);

  /**
   * @brief A swap constructor
   *
   * @param[in] variables is the vector of the variables
   * @param[in] parameters is the vector of the parameters
   * @param[in] dynamics is the vector of the dynamic laws
   */
  ODE(std::vector<SymbolicAlgebra::Symbol<T>> &&variables,
      std::vector<SymbolicAlgebra::Symbol<T>> &&parameters,
      std::vector<SymbolicAlgebra::Expression<T>> &&dynamics,
      const SymbolicAlgebra::Symbol<T> &time_variable);

  /**
   * @brief A constructor
   *
   * @param[in] dynamics is the map that relates variables and dynamics laws
   */
  ODE(const std::map<SymbolicAlgebra::Symbol<T>,
                     SymbolicAlgebra::Expression<T>> &dynamics,
      const SymbolicAlgebra::Symbol<T> &time_variable);

  /**
   * @brief A constructor
   *
   * @param[in] dynamics is the map that relates variables and dynamics laws
   * @param[in] parameters is the vector of the parameters
   */
  ODE(const std::map<SymbolicAlgebra::Symbol<T>,
                     SymbolicAlgebra::Expression<T>> &dynamics,
      const std::vector<SymbolicAlgebra::Symbol<T>> &parameters,
      const SymbolicAlgebra::Symbol<T> &time_variable);

  /**
   * @brief A constructor
   *
   * @param[in] dynamics is the map that relates variables and dynamic laws
   * @param[in] parameters is the vector of the parameters
   */
  ODE(const std::map<SymbolicAlgebra::Symbol<T>,
                     SymbolicAlgebra::Expression<T>> &&dynamics,
      std::vector<SymbolicAlgebra::Symbol<T>> &&parameters,
      const SymbolicAlgebra::Symbol<T> &time_variable);

  /**
   * @brief A constructor
   *
   * @param[in] variables is the vector of the variables
   * @param[in] dynamics is the vector of the dynamic laws
   */
  ODE(const std::vector<SymbolicAlgebra::Symbol<T>> &variables,
      const std::vector<SymbolicAlgebra::Expression<T>> &dynamics,
      const SymbolicAlgebra::Symbol<T> &time_variable);

  /**
   * @brief A swap constructor
   *
   * @param[in] variables is the vector of the variables
   * @param[in] dynamics is the vector of the dynamic laws
   */
  ODE(std::vector<SymbolicAlgebra::Symbol<T>> &&variables,
      std::vector<SymbolicAlgebra::Expression<T>> &&dynamics,
      const SymbolicAlgebra::Symbol<T> &time_variable);

  /**
   * @brief A copy constructor
   *
   * @param[in] orig the original ODE
   */
  ODE(const ODE<T> &orig);

  /**
   * @brief A swap constructor
   *
   * @param[in] orig the original ODE
   */
  ODE(ODE<T> &&orig);

  /**
   * @brief A copy operator
   *
   * @param orig the original ODE
   * @return a reference to the updated object
   */
  ODE &operator=(const ODE<T> &orig);

  /**
   * @brief A move operator
   *
   * @param orig the original ODE
   * @return a reference to the updated object
   */
  ODE &operator=(ODE<T> &&orig);
};

template<typename T>
ODE<T>::ODE(const std::vector<SymbolicAlgebra::Symbol<T>> &variables,
            const std::vector<SymbolicAlgebra::Symbol<T>> &parameters,
            const std::vector<SymbolicAlgebra::Expression<T>> &dynamics,
            const SymbolicAlgebra::Symbol<T> &time_variable):
    ContinuousSystem<T>(variables, parameters, dynamics, time_variable)
{
}

template<typename T>
ODE<T>::ODE(std::vector<SymbolicAlgebra::Symbol<T>> &&variables,
            std::vector<SymbolicAlgebra::Symbol<T>> &&parameters,
            std::vector<SymbolicAlgebra::Expression<T>> &&dynamics,
            const SymbolicAlgebra::Symbol<T> &time_variable):
    ContinuousSystem<T>(std::move(variables), std::move(parameters),
                        std::move(dynamics), time_variable)
{
}

template<typename T>
ODE<T>::ODE(): ContinuousSystem<T>()
{
}

template<typename T>
ODE<T>::ODE(const std::map<SymbolicAlgebra::Symbol<T>,
                           SymbolicAlgebra::Expression<T>> &dynamics,
            const SymbolicAlgebra::Symbol<T> &time_variable):
    ContinuousSystem<T>(dynamics, time_variable)
{
}

template<typename T>
ODE<T>::ODE(const std::map<SymbolicAlgebra::Symbol<T>,
                           SymbolicAlgebra::Expression<T>> &dynamics,
            const std::vector<SymbolicAlgebra::Symbol<T>> &parameters,
            const SymbolicAlgebra::Symbol<T> &time_variable):
    ContinuousSystem<T>(dynamics, parameters, time_variable)
{
}

template<typename T>
ODE<T>::ODE(const std::map<SymbolicAlgebra::Symbol<T>,
                           SymbolicAlgebra::Expression<T>> &&dynamics,
            std::vector<SymbolicAlgebra::Symbol<T>> &&parameters,
            const SymbolicAlgebra::Symbol<T> &time_variable):
    ContinuousSystem<T>(std::move(dynamics), std::move(parameters),
                        time_variable)
{
}

template<typename T>
ODE<T>::ODE(const std::vector<SymbolicAlgebra::Symbol<T>> &variables,
            const std::vector<SymbolicAlgebra::Expression<T>> &dynamics,
            const SymbolicAlgebra::Symbol<T> &time_variable):
    ContinuousSystem<T>(variables, dynamics, time_variable)
{
}

template<typename T>
ODE<T>::ODE(std::vector<SymbolicAlgebra::Symbol<T>> &&variables,
            std::vector<SymbolicAlgebra::Expression<T>> &&dynamics,
            const SymbolicAlgebra::Symbol<T> &time_variable):
    ContinuousSystem<T>(std::move(variables), std::move(dynamics),
                        time_variable)
{
}

template<typename T>
ODE<T>::ODE(const ODE<T> &orig): ContinuousSystem<T>(orig)
{
}

template<typename T>
ODE<T>::ODE(ODE<T> &&orig): ContinuousSystem<T>(std::move(orig))
{
}

template<typename T>
inline ODE<T> &ODE<T>::operator=(const ODE<T> &orig)
{
  ContinuousSystem<T>::operator=(orig);

  return *this;
}

template<typename T>
inline ODE<T> &ODE<T>::operator=(ODE<T> &&orig)
{
  ContinuousSystem<T>::operator=(std::move(orig));

  return *this;
}

#endif // _DIFFERENTIAL_SYSTEM_H_