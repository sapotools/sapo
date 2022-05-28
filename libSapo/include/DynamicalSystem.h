/**
 * @file DynamicalSystem.h
 * Represent dynamic system laws
 *
 * @author Alberto Casagrande <acasagrande@units.it>
 * @version 0.1
 */

#ifndef DYNAMICS_H_
#define DYNAMICS_H_

#include <map>
#include <vector>
#include <sstream>

#include "Bundle.h"
#include "SymbolicAlgebra.h"
#include "Polytope.h"
#include "SetsUnion.h"

#include "LinearAlgebraIO.h"

#include "STL/Atom.h"

/**
 * @brief Approach to evaluate the image of a bundle
 *
 * The are two different approaches to evaluate the
 * image of a bundle through a polynomial function:
 * 1. One-For-One: the boundaries of any parallelotope
 *                 image in the bundle are evaluated by
 *                 exclusively considering the original
 *                 parallelotope itself
 * 2. All-For-One: the boundaries of any parallelotope
 *                 image in the bundle are evaluated by
 *                 exploiting all the bundle templates
 */
typedef enum {
  ONE_FOR_ONE, /* One-For-One */
  ALL_FOR_ONE  /* All-For-One */
} transformation_mode;

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
  const std::vector<SymbolicAlgebra::Symbol<T>> &variables() const
  {
    return _variables;
  }

  /**
   * @brief Get the system parameters
   *
   * @return a reference to the vector of the system parameters
   */
  const std::vector<SymbolicAlgebra::Symbol<T>> &parameters() const
  {
    return _parameters;
  }

  /**
   * @brief Get the system dynamic laws
   *
   * @return a reference to the system dynamic law vector
   */
  const std::vector<SymbolicAlgebra::Expression<T>> &dynamics() const
  {
    return _dynamics;
  }

  /**
   * @brief Get the i-th variable
   *
   * @return a reference to the i-th variable
   */
  const SymbolicAlgebra::Symbol<T> &variable(const unsigned int i) const
  {
    return _variables[i];
  }

  /**
   * @brief Get the i-th parameter
   *
   * @return a reference to the i-th parameter
   */
  const SymbolicAlgebra::Symbol<T> &parameter(const unsigned int i) const
  {
    return _parameters[i];
  }

  /**
   * @brief Get the i-th dynamic law
   *
   * @return a reference to the i-th dynamic law
   */
  const SymbolicAlgebra::Expression<T> &dynamic(const unsigned int i) const
  {
    return _dynamics[i];
  }

  /**
   * @brief Get the dynamic law associated to a variable
   *
   * @return a reference to the dynamic law associated to a variable
   */
  const SymbolicAlgebra::Expression<T> &
  operator[](const SymbolicAlgebra::Symbol<T> &variable) const
  {
    return _dynamics[_var_index.at[variable.get_id()]];
  }

  /**
   * @brief Transform a bundle according with the system dynamics
   *
   * @param bundle is the bundle to be transformed
   * @param t_mode is the mode used to compute the transformation
   * @return an over-approximation of the bundle transformed
   *         by the dynamic laws
   */
  inline Bundle transform(const Bundle &bundle,
                          const transformation_mode t_mode = ALL_FOR_ONE) const
  {
    if (_parameters.size() != 0) {
      throw std::domain_error("The parameter set has not been specified");
    }

    return this->transform(bundle, Polytope(), t_mode);
  }

  /**
   * @brief Transform a bundle according with the system dynamics
   *
   * @param bundle is the bundle to be transformed
   * @param parameter_set is the parameter set
   * @param t_mode is the mode used to compute the transformation
   * @return an over-approximation of the bundle transformed
   *         by the dynamic laws
   */
  Bundle transform(const Bundle &bundle, const Polytope &parameter_set,
                   const transformation_mode t_mode = ALL_FOR_ONE) const;
};

/**
 * @brief Perform the parametric synthesis for an atom
 *
 * This method computes a set of parameters such that the
 * transformation from a bundle satisfies the provided atom.
 *
 * @param[in] ds is the dynamical system
 * @param[in] bundle is the set from which the evolution occurs
 * @param[in] parameter_set is the initial parameter set
 * @param[in] atom is the specification to be satisfied
 * @return a subset of `parameter_set` such that the
 *         transformation from `bundle` through the dynamical
 *         system when the parameters are in the returned set
 *         satisfies the atomic formula
 */
SetsUnion<Polytope> synthesize(const DynamicalSystem<double> &ds,
                               const Bundle &bundle,
                               const SetsUnion<Polytope> &parameter_set,
                               const std::shared_ptr<STL::Atom> atom);

#endif // DYNAMICS_H_