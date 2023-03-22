/**
 * @file Model.h
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @author Alberto Casagrande <acasagrande@units.it>
 * @brief Discrete-time dynamical models
 * @version 0.1
 * @date 2015-10-15
 *
 * @copyright Copyright (c) 2015-2022
 */

#ifndef MODEL_H_
#define MODEL_H_

#include <sstream>

#include "SymbolicAlgebra.h"
#include "Bundle.h"
#include "Polytope.h"
#include "SetsUnion.h"
#include "STL/STL.h"

#include "DiscreteSystem.h"

/**
 * @brief Dynamical system models
 *
 * This class represents dynamical system models
 */
class Model
{

protected:
  std::shared_ptr<Bundle> _init_set; //!< The initial set
  SetsUnion<Polytope> _param_set;    //!< The parameter set

  std::shared_ptr<STL::STL> _spec; //!< a specification

  LinearSystem<double> _assumptions; //!< The assumptions
  LinearSystem<double> _invariant;   //!< The candidate invariant

  std::string _name; //!< The model name

  /**
   * @brief Validate the validate constructor parameters
   *
   * This method validate constructor parameters and throws a
   * `std::domain_error` exception if they are not consistent.
   *
   * @param variables is the vector of the variables
   * @param parameters is the vector of the parameter
   * @param dynamics is the vector of the dynamic laws
   * @param init_set is the initial set
   * @param param_set is the set of the parameter
   */
  static void validate_parameters(
      const std::vector<SymbolicAlgebra::Symbol<double>> &variables,
      const std::vector<SymbolicAlgebra::Symbol<double>> &parameters,
      const std::vector<SymbolicAlgebra::Expression<double>> &dynamics,
      const Bundle &init_set, const SetsUnion<Polytope> &parameter_set);

public:
  /**
   * @brief Constructor
   *
   * @param variables is the vector of the variables
   * @param dynamics is the vector of the dynamic laws
   * @param init_set is the initial set
   * @param name is the model name
   */
  Model(const Bundle &init_set, const std::string name = "Unknown");

  /**
   * @brief Constructor
   *
   * @param variables is the vector of the variables
   * @param parameters is the vector of the parameter
   * @param dynamics is the vector of the dynamic laws
   * @param init_set is the initial set
   * @param param_set is the set of the parameter
   * @param name is the model name
   */
  Model(const Bundle &init_set, const SetsUnion<Polytope> &param_set,
        const std::string name = "Unknown");

  /**
   * @brief Get the model name
   *
   * @return a reference to the model name
   */
  inline const std::string &name() const
  {
    return this->_name;
  }

  /**
   * @brief Get the model dynamical system
   *
   * @return a reference to the model dynamical system
   */
  virtual inline const DynamicalSystem<double> &dynamical_system() const = 0;

  /**
   * @brief Get the model variables
   *
   * @return a reference to the model variables
   */
  inline const std::vector<SymbolicAlgebra::Symbol<double>> &variables() const
  {
    return this->dynamical_system().variables();
  }

  /**
   * @brief Get the model parameters
   *
   * @return a reference to the model parameters
   */
  inline const std::vector<SymbolicAlgebra::Symbol<double>> &parameters() const
  {
    return this->dynamical_system().parameters();
  }

  /**
   * @brief Get the model dynamic laws
   *
   * @return a reference to the model dynamic laws
   */
  inline const std::vector<SymbolicAlgebra::Expression<double>> &
  dynamics() const
  {
    return this->dynamical_system().dynamics();
  }

  /**
   * @brief Get the model initial set
   *
   * @return a reference to the model initial set
   */
  inline const std::shared_ptr<Bundle> initial_set() const
  {
    return this->_init_set;
  }

  /**
   * @brief Get the model parameter set
   *
   * @return a reference to the model parameter set
   */
  inline const SetsUnion<Polytope> &parameter_set() const
  {
    return this->_param_set;
  }

  /**
   * @brief Get the model specification
   *
   * @return a reference to the model specification
   */
  inline const std::shared_ptr<STL::STL> specification() const
  {
    return this->_spec;
  }

  /**
   * @brief Set the model specification
   *
   * @return a reference to the model specification
   */
  Model &set_specification(std::shared_ptr<STL::STL> specification);

  /**
   * @brief Get the model assumptions
   *
   * @return a reference to the model assumptions
   */
  inline const LinearSystem<double> &assumptions() const
  {
    return this->_assumptions;
  }

  /**
   * @brief Get the candidate model invariant
   *
   * @return a reference to the candidate model invariant
   */
  inline const LinearSystem<double> &invariant() const
  {
    return this->_invariant;
  }

  /**
   * @brief Set the model assumptions
   *
   * @return a reference to the model assumptions
   */
  Model &set_assumptions(const LinearSystem<double> &assumptions);

  /**
   * @brief Set the candidate model invariant
   *
   * @return a reference to the candidate model invariant
   */
  Model &set_invariant(const LinearSystem<double> &invariant);

  /**
   * @brief Get the number of variables in the model
   *
   * @return the number of variables in the model
   */
  inline size_t dim() const
  {
    return variables().size();
  }

  /**
   * @brief Get the model type
   *
   * @return the model type
   */
  virtual inline DynamicType type() const
  {
    return DynamicType::UNDEFINED;
  }

  /**
   * @brief The destroyer
   */
  virtual ~Model() {}
};

/**
 * @brief Dynamical system models
 *
 * This class represents dynamical system models
 */
class DiscreteModel : public Model
{

  DiscreteSystem<double> _discrete_system; //!< A discrete system

public:
  /**
   * @brief Constructor
   *
   * @param variables is the vector of the variables
   * @param dynamics is the vector of the dynamic laws
   * @param init_set is the initial set
   * @param name is the model name
   */
  DiscreteModel(
      const std::vector<SymbolicAlgebra::Symbol<double>> &variables,
      const std::vector<SymbolicAlgebra::Expression<double>> &dynamics,
      const Bundle &init_set, const std::string name = "Unknown");

  /**
   * @brief Constructor
   *
   * @param variables is the vector of the variables
   * @param parameters is the vector of the parameter
   * @param dynamics is the vector of the dynamic laws
   * @param init_set is the initial set
   * @param param_set is the set of the parameter
   * @param name is the model name
   */
  DiscreteModel(
      const std::vector<SymbolicAlgebra::Symbol<double>> &variables,
      const std::vector<SymbolicAlgebra::Symbol<double>> &parameters,
      const std::vector<SymbolicAlgebra::Expression<double>> &dynamics,
      const Bundle &init_set, const SetsUnion<Polytope> &param_set,
      const std::string name = "Unknown");

  /**
   * @brief Get the model discrete system
   *
   * @return a reference to the model discrete system
   */
  inline const DynamicalSystem<double> &dynamical_system() const
  {
    return this->_discrete_system;
  }

  /**
   * @brief Get the model type
   *
   * @return the model type
   */
  inline DynamicType type() const
  {
    return _discrete_system.type();
  }
};

#endif /* MODEL_H_ */
