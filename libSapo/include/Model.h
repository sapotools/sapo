/**
 * @file Model.h
 * Represent a discrete-time (eventually parametric) dynamical system
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#ifndef MODEL_H_
#define MODEL_H_

#include <sstream> 

#include "SymbolicAlgebra.h"
#include "Bundle.h"
#include "PolytopesUnion.h"
#include "STL.h"

#include "LinearAlgebraIO.h"

/**
 * @brief Dynamical system models
 * 
 * This class represents dynamical system models
 */
class Model
{

protected:
  std::vector<SymbolicAlgebra::Symbol<>> _vars;     //!< Variables
  std::vector<SymbolicAlgebra::Symbol<>> _params;   //!< Parameters
  std::vector<SymbolicAlgebra::Expression<>> _dynamics; //!< Dynamic laws

  std::shared_ptr<Bundle> _init_set; //!< Initial set
  PolytopesUnion _param_set;  //!< Parameter set

  std::shared_ptr<STL::STL> _spec;  //!< Specification

  LinearSystem _assumptions;  //!< Assumptions

  std::string _name;  //!< Name of the model

public:

  /**
   * @brief Constructor
   * 
   * @param variables is the vector of the variables
   * @param dynamics is the vector of the dynamic laws
   * @param init_set is the initial set 
   * @param name is the model name
   */
  Model(const std::vector<SymbolicAlgebra::Symbol<>> &variables,
        const std::vector<SymbolicAlgebra::Expression<>> &dynamics,
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
  Model(const std::vector<SymbolicAlgebra::Symbol<>> &variables,
        const std::vector<SymbolicAlgebra::Symbol<>> &parameters,
        const std::vector<SymbolicAlgebra::Expression<>> &dynamics,
        const Bundle &init_set, const PolytopesUnion &param_set,
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
   * @brief Get the model variables
   * 
   * @return a reference to the model variables
   */
  inline const std::vector<SymbolicAlgebra::Symbol<>> &variables() const
  {
    return this->_vars;
  }

  /**
   * @brief Get the model parameters
   * 
   * @return a reference to the model parameters
   */
  inline const std::vector<SymbolicAlgebra::Symbol<>> &parameters() const
  {
    return this->_params;
  }

  /**
   * @brief Get the model dynamic laws
   * 
   * @return a reference to the model dynamic laws
   */
  inline const std::vector<SymbolicAlgebra::Expression<>> &dynamics() const
  {
    return this->_dynamics;
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
  inline const PolytopesUnion& parameter_set() const
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
  inline const LinearSystem& assumptions() const
  {
    return this->_assumptions;
  }

  /**
   * @brief Set the model assumptions 
   * 
   * @return a reference to the model assumptions
   */
  Model &set_assumptions(const LinearSystem &assumptions);
};

#endif /* MODEL_H_ */
