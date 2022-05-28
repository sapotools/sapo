/**
 * @file Sapo.h
 * Core of Sapo tool.
 * Here the reachable set and the parameter synthesis are done.
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#ifndef SAPO_H_
#define SAPO_H_

#include <string>

#include "Always.h"
#include "Atom.h"
#include "BaseConverter.h"
#include "Bundle.h"
#include "Conjunction.h"
#include "Disjunction.h"
#include "Eventually.h"
#include "Flowpipe.h"
#include "Polytope.h"
#include "SetsUnion.h"
#include "Model.h"
#include "STL.h"
#include "Until.h"

#include "DynamicalSystem.h"

#include "ProgressAccounter.h"

/**
 * @brief The computation orchestrator
 */
class Sapo
{
public:
  transformation_mode t_mode; //!< transformation mode (OFO or AFO)
  double decomp_weight;       //!< decomposition weight
  unsigned int decomp;        //!< number of decompositions (0: none, >0: yes)
  std::string plot; //!< the name of the file were to plot the reach set
  unsigned int time_horizon;      //!< the computation time horizon
  unsigned int max_param_splits;  //!< maximum number of splits in synthesis
  unsigned int num_of_pre_splits; //!< number of pre-splits in synthesis
  double max_bundle_magnitude; //!< maximum versor magnitude for single bundle

private:
  const DynamicalSystem<double>
      _dynamical_system; //!< the investigated dynamical system

  const LinearSystem assumptions;

  /**
   * Parameter synthesis w.r.t. an always formula
   *
   * @param[in] init_set is bundle with the initial set
   * @param[in] parameterSet is set of parameters
   * @param[in] sigma is STL always formula
   * @returns a refined sets of parameters
   */
  template<typename T>
  SetsUnion<Polytope> transition_and_synthesis(Bundle init_set,
                                          const SetsUnion<Polytope> &pSet,
                                          const std::shared_ptr<T> formula,
                                          const int time) const
  {
    SetsUnion<Polytope> result;

    // init_set.intersect(this->assumptions);

    for (auto p_it = pSet.begin(); p_it != pSet.end(); ++p_it) {
      // transition by using the n-th polytope of the parameter set
      Bundle reached_set = _dynamical_system.transform(init_set, *p_it);

      // guarantee the assumptions
      reached_set.intersect_with(this->assumptions);

      result.add(synthesize(reached_set, pSet, formula, time + 1));
    }

    return result;
  }

  /**
   * Parameter synthesis for atomic formulas
   *
   * @param[in] init_set is the initial set
   * @param[in] pSet the current set of parameters
   * @param[in] formula is an STL atomic formula providing the specification
   * @returns refined parameter set
   */
  SetsUnion<Polytope> synthesize(const Bundle &reachSet, const SetsUnion<Polytope> &pSet,
                            const std::shared_ptr<STL::Atom> formula) const;

  /**
   * Parameter synthesis for conjunctions
   *
   * @param[in] init_set is the initial set
   * @param[in] pSet the current parameter set
   * @param[in] conj is an STL conjunction providing the specification
   * @returns refined parameter set
   */
  SetsUnion<Polytope>
  synthesize(const Bundle &reachSet, const SetsUnion<Polytope> &pSet,
             const std::shared_ptr<STL::Conjunction> formula) const;

  /**
   * Parameter synthesis for disjunctions
   *
   * @param[in] init_set is the initial set
   * @param[in] pSet the current parameter set
   * @param[in] conj is an STL disjunction providing the specification
   * @returns refined parameter set
   */
  SetsUnion<Polytope>
  synthesize(const Bundle &reachSet, const SetsUnion<Polytope> &pSet,
             const std::shared_ptr<STL::Disjunction> formula) const;

  /**
   * Parameter synthesis for until formulas
   *
   * @param[in] init_set is the initial set
   * @param[in] pSet the current set of parameters
   * @param[in] formula is an STL until formula providing the specification
   * @param[in] time is the time of the current evaluation
   * @returns refined parameter set
   */
  SetsUnion<Polytope> synthesize(const Bundle &reachSet, const SetsUnion<Polytope> &pSet,
                            const std::shared_ptr<STL::Until> formula,
                            const int time) const;

  /**
   * Parameter synthesis for always formulas
   *
   * @param[in] init_set is the initial set
   * @param[in] pSet the current set of parameters
   * @param[in] formula is an STL always formula providing the specification
   * @param[in] time is the time of the current evaluation
   * @returns refined parameter set
   */
  SetsUnion<Polytope> synthesize(const Bundle &reachSet, const SetsUnion<Polytope> &pSet,
                            const std::shared_ptr<STL::Always> formula,
                            const int time) const;
  /**
   * Parameter synthesis for the eventually formulas
   *
   * @param[in] init_set is the initial set
   * @param[in] pSet the current parameter set
   * @param[in] formula is an STL eventually formula providing the
   * specification
   * @returns refined parameter set
   */
  SetsUnion<Polytope> synthesize(const Bundle &reachSet, const SetsUnion<Polytope> &pSet,
                            const std::shared_ptr<STL::Eventually> ev) const;

public:
  /**
   * Constructor that instantiates Sapo
   *
   * @param[in] model model to analyze
   */
  Sapo(const Model &model);

  /**
   * @brief Get the dynamical system
   *
   * @return the dynamical system
   */
  inline const DynamicalSystem<double> &dynamical_system() const
  {
    return _dynamical_system;
  }

  /**
   * Reachable set computation
   *
   * @param[in] init_set is the initial set
   * @param[in] k is the time horizon
   * @param[in,out] accounter accounts for the computation progress
   * @returns the reached flowpipe
   */
  Flowpipe reach(Bundle init_set, unsigned int k,
                 ProgressAccounter *accounter = NULL) const;

  /**
   * Reachable set computation for parametric dynamical systems
   *
   * @param[in] init_set is the initial set
   * @param[in] pSet is the set of parameters
   * @param[in] k is the time horizon
   * @param[in,out] accounter accounts for the computation progress
   * @returns the reached flowpipe
   */
  Flowpipe reach(Bundle init_set, const SetsUnion<Polytope> &pSet, unsigned int k,
                 ProgressAccounter *accounter = NULL) const;

  /**
   * Parameter synthesis method
   *
   * @param[in] init_set is the initial set
   * @param[in] pSet is the set of parameters
   * @param[in] formula is an STL specification for the model
   * @param[in,out] accounter accounts for the computation progress
   * @returns a parameter set refined according with `formula`
   */
  SetsUnion<Polytope> synthesize(Bundle init_set, const SetsUnion<Polytope> &pSet,
                            const std::shared_ptr<STL::STL> formula,
                            ProgressAccounter *accounter = NULL) const;

  /**
   * Parameter synthesis with splits
   *
   * @param[in] init_set is the initial set
   * @param[in] pSet is the current parameter sets
   * @param[in] formula is an STL formula providing the specification
   * @param[in] max_splits maximum number of splits of the original
   *                       parameter set to identify a non-null solution
   * @param[in] num_of_pre_splits is number of splits to be performed before
   *                             the computation
   * @param[in,out] accounter accounts for the computation progress
   * @returns the list of refined parameter sets
   */
  std::list<SetsUnion<Polytope>>
  synthesize(Bundle init_set, const SetsUnion<Polytope> &pSet,
             const std::shared_ptr<STL::STL> formula,
             const unsigned int max_splits,
             const unsigned int num_of_pre_splits = 0,
             ProgressAccounter *accounter = NULL) const;
};

/**
 * Test whether all the sets in a list are empty
 *
 * @param[in] sets ia a list of sets
 * @returns `true` if and only if all the sets in `sets` are empty
 */
template<class T>
bool every_set_is_empty(const std::list<T> &sets)
{
  for (auto s_it = std::begin(sets); s_it != std::end(sets); ++s_it) {
    if (!s_it->is_empty()) {
      return false;
    }
  }

  return true;
}

#endif /* SAPO_H_ */
