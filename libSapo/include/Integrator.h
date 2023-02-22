/**
 * @file Integrator.h
 * @author Alberto Casagrande <acasagrande@units.it>
 * @brief Represent differential systems
 * @version 0.1
 * @date 2023-02-22
 *
 * @copyright Copyright (c) 2023
 */

#ifndef _INTEGRATOR_H_
#define _INTEGRATOR_H_

#include "DifferentialSystem.h"
#include "ContinuousSystem.h"
#include "DiscreteSystem.h"

class ODEIntegrator
{
public:
  ODEIntegrator() {}

  /**
   * @brief Integrate an ODE with a symbolic time step
   *
   * @tparam T is the ODE constant type
   * @param ode is the ODE to be integrated
   * @param time_step is the time step variable
   * @return the continuous system associated to the ODE
   */
  template<typename T>
  ContinuousSystem<T>
  operator()(const ODE<T> &ode,
             const SymbolicAlgebra::Symbol<T> &time_step) const;

  /**
   * @brief Integrate an ODE with a symbolic time step
   *
   * @tparam T is the  ODE constant type
   * @param ode is the ODE to be integrated
   * @param time_step is the time step variable
   * @return the continuous system associated to the ODE
   */
  template<typename T>
  DiscreteSystem<T> operator()(const ODE<T> &ode, const T &time_step) const;
};

class RungeKutta4Integrator : public ODEIntegrator
{
  /**
   * @brief Integrate an ODE by using the 4th degree Runge-Kutta method
   *
   * @tparam T is the  ODE constant type
   * @tparam DELTA_TYPE is the type of the delta step
   * @param variables is the variable vector
   * @param ode is the ODE vector
   * @param time is the time variable
   * @param time_step is the integration time step
   * @return the integration of `ode` according to the 4th degree
   *         Runge-Kutta method
   */
  template<typename T, typename DELTA_TYPE>
  static std::vector<SymbolicAlgebra::Expression<T>>
  runge_kutta4(const std::vector<SymbolicAlgebra::Symbol<T>> &variables,
               const std::vector<SymbolicAlgebra::Expression<T>> &dynamics,
               const SymbolicAlgebra::Symbol<T> &time,
               const DELTA_TYPE &time_step);

public:
  RungeKutta4Integrator(): ODEIntegrator() {}

  /**
   * @brief Integrate an ODE with a symbolic time step
   *
   * @tparam T is the  ODE constant type
   * @param ode is the ODE to be integrated
   * @param time_step is the time step variable
   * @return the continuous system associated to the ODE
   */
  template<typename T>
  ContinuousSystem<T>
  operator()(const ODE<T> &ode,
             const SymbolicAlgebra::Symbol<T> &time_step) const;

  /**
   * @brief Integrate an ODE with a symbolic time step
   *
   * @tparam T is the  ODE constant type
   * @param ode is the ODE to be integrated
   * @param time_step is the time step variable
   * @return the continuous system associated to the ODE
   */
  template<typename T>
  DiscreteSystem<T> operator()(const ODE<T> &ode, const T &time_step) const;
};

template<typename T, typename DELTA_TYPE>
std::vector<SymbolicAlgebra::Expression<T>>
RungeKutta4Integrator::runge_kutta4(
    const std::vector<SymbolicAlgebra::Symbol<T>> &variables,
    const std::vector<SymbolicAlgebra::Expression<T>> &dynamics,
    const SymbolicAlgebra::Symbol<T> &time, const DELTA_TYPE &time_step)
{
  using namespace SymbolicAlgebra;
  using namespace LinearAlgebra;

  const auto &k1 = dynamics;

  // get k2
  auto k2 = dynamics;
  {
    typename Expression<T>::replacement_type repl;
    for (size_t i = 0; i < variables.size(); ++i) {
      repl[variables[i]] = variables[i] + k1[i] * time_step / 2;
    }
    repl[time] = time + Expression<T>(time_step) / T(2);

    for (auto &k: k2) {
      k.replace(repl);
    }
  }

  // get k3
  auto k3 = dynamics;
  {
    typename Expression<T>::replacement_type repl;
    for (size_t i = 0; i < variables.size(); ++i) {
      repl[variables[i]] = variables[i] + k2[i] * time_step / 2;
    }
    repl[time] = time + Expression<T>(time_step) / T(2);

    for (auto &k: k3) {
      k.replace(repl);
    }
  }

  // get k4
  auto k4 = dynamics;
  {
    typename Expression<T>::replacement_type repl;
    for (size_t i = 0; i < variables.size(); ++i) {
      repl[variables[i]] = variables[i] + k3[i] * time_step;
    }
    repl[time] = time + time_step;

    for (auto &k: k4) {
      k.replace(repl);
    }
  }

  std::vector<Expression<T>> rk_dyn(variables.size());

  for (size_t i = 0; i < variables.size(); ++i) {
    rk_dyn[i] = variables[i]
                + (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) * time_step / 6;
  }

  return rk_dyn;
}

template<typename T>
ContinuousSystem<T> RungeKutta4Integrator::operator()(
    const ODE<T> &ode, const SymbolicAlgebra::Symbol<T> &time_step) const
{
  using namespace SymbolicAlgebra;

  if (ode.time_variable() == time_step) {
    SAPO_ERROR("ODE time variable and the time step "
                   << "must differ",
               std::domain_error);
  }

  auto dynamics = runge_kutta4(ode.variables(), ode.dynamics(),
                               ode.time_variable(), time_step);

  std::vector<Symbol<T>> variables(ode.variables());
  std::vector<Symbol<T>> parameters(ode.parameters());

  if (!ode.is_time_independent()) {
    variables.push_back(ode.time_variable());
    dynamics.push_back(ode.time_variable() + time_step);
  }

  return {std::move(variables), std::move(parameters), std::move(dynamics),
          time_step};
}

template<typename T>
DiscreteSystem<T> RungeKutta4Integrator::operator()(const ODE<T> &ode,
                                                    const T &time_step) const
{
  using namespace SymbolicAlgebra;

  auto dynamics = runge_kutta4(ode.variables(), ode.dynamics(),
                               ode.time_variable(), time_step);

  std::vector<Symbol<T>> variables(ode.variables());
  std::vector<Symbol<T>> parameters(ode.parameters());

  if (!ode.is_time_independent()) {
    variables.push_back(ode.time_variable());
    dynamics.push_back(ode.time_variable() + time_step);
  }

  return {std::move(variables), std::move(parameters), std::move(dynamics)};
}

#endif // _INTEGRATOR_H_