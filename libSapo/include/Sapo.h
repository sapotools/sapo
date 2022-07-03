/**
 * @file Sapo.h
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @author Alberto Casagrande <acasagrande@units.it>
 * @brief Sapo analysis methods
 * @version 0.1
 * @date 2015-10-14
 *
 * @copyright Copyright (c) 2015-2022
 */

#ifndef SAPO_H_
#define SAPO_H_

#include <string>
#include <climits>

#include "BaseConverter.h"
#include "Bundle.h"
#include "Flowpipe.h"
#include "Polytope.h"
#include "SetsUnion.h"
#include "Model.h"
#include "STL/Always.h"
#include "STL/Atom.h"
#include "STL/Conjunction.h"
#include "STL/Disjunction.h"
#include "STL/Eventually.h"
#include "STL/STL.h"
#include "STL/Until.h"

#include "Evolver.h"

#include "ProgressAccounter.h"

/**
 * @brief Result of an invariant validation
 */
typedef struct {
  bool validated;             //!< invariant candidate has been validated
  unsigned int k_induction;   //!< k-induction level
  Flowpipe flowpipe;          //!< system flowpipe
  SetsUnion<Bundle> r_approx; //!< over-approximation of the reached set
} InvariantValidationResult;

/**
 * @brief The computation orchestrator
 */
class Sapo
{

public:
  double decomp_weight; //!< decomposition weight
  unsigned int decomp;  //!< number of decompositions (0: none, >0: yes)
  std::string plot;     //!< the name of the file were to plot the reach set
  unsigned int time_horizon;      //!< the computation time horizon
  unsigned int max_param_splits;  //!< maximum number of splits in synthesis
  unsigned int num_of_pre_splits; //!< number of pre-splits in synthesis
  double max_bundle_magnitude; //!< maximum versor magnitude for single bundle

  /**
   * @brief Approximation used during invariant validation
   */
  enum invariantApproxType {
    NO_APPROX,  // no approximation is used (listing)
    CHAIN_JOIN, // chain-join approximation (packaging)
    FULL_JOIN   // full-join approximation (merging)
  };

  invariantApproxType inv_approximation;

private:
  Evolver<double> _evolver; //!< the dynamical system evolver

  const LinearSystem assumptions;

  /**
   * @brief Parameter synthesis
   *
   * @tparam T is the STL specification type
   * @param init_set is the initial set
   * @param pSet is the set of parameters
   * @param formula is the specification
   * @param epoch_horizon
   * @return SetsUnion<Polytope>
   */
  template<typename T>
  SetsUnion<Polytope>
  transition_and_synthesis(Bundle init_set, const SetsUnion<Polytope> &pSet,
                           const std::shared_ptr<T> formula,
                           const int epoch_horizon) const
  {
    SetsUnion<Polytope> result;

    // init_set.intersect(this->assumptions);

    for (auto p_it = pSet.begin(); p_it != pSet.end(); ++p_it) {
      // transition by using the n-th polytope of the parameter set
      Bundle reached_set = _evolver(init_set, *p_it);

      // guarantee the assumptions
      reached_set.intersect_with(this->assumptions);

      result.update(synthesize(reached_set, pSet, formula, epoch_horizon + 1));
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
  SetsUnion<Polytope>
  synthesize(const Bundle &reachSet, const SetsUnion<Polytope> &pSet,
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
  SetsUnion<Polytope> synthesize(const Bundle &reachSet,
                                 const SetsUnion<Polytope> &pSet,
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
  SetsUnion<Polytope> synthesize(const Bundle &reachSet,
                                 const SetsUnion<Polytope> &pSet,
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
  SetsUnion<Polytope>
  synthesize(const Bundle &reachSet, const SetsUnion<Polytope> &pSet,
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
    return _evolver.dynamical_system();
  }

  /**
   * @brief Get the dynamical system evolver
   *
   * @return the dynamical system evolver
   */
  inline const Evolver<double> &evolver() const
  {
    return _evolver;
  }

  /**
   * @brief Set the evolver mode
   *
   * @param mode is the to-be-set evolver mode
   */
  inline void set_evolver_mode(const Evolver<double>::evolver_mode mode)
  {
    _evolver.mode = mode;
  }

  /**
   * Reachable set computation
   *
   * @param[in] init_set is the initial set
   * @param[in] epoch_horizon is the time horizon
   * @param[in,out] accounter accounts for the computation progress
   * @returns the reached flowpipe
   */
  Flowpipe reach(Bundle init_set, unsigned int epoch_horizon,
                 ProgressAccounter *accounter = NULL) const;

  /**
   * Reachable set computation for parametric dynamical systems
   *
   * @param[in] init_set is the initial set
   * @param[in] pSet is the set of parameters
   * @param[in] epoch_horizon is the epoch horizon
   * @param[in,out] accounter accounts for the computation progress
   * @returns the reached flowpipe
   */
  Flowpipe reach(Bundle init_set, const SetsUnion<Polytope> &pSet,
                 unsigned int epoch_horizon,
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
  SetsUnion<Polytope> synthesize(Bundle init_set,
                                 const SetsUnion<Polytope> &pSet,
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

  /**
   * @brief Try to establish whether a set is an invariant
   *
   * @param init_set is the initial set
   * @param pSet is the parameter set
   * @param invariant_candidate is the candidate invariant
   * @param epoch_horizon is the maximum investigated epoch.
   *        If `epoch_horizon` is set to be 0, then the
   *        computation does not end until the correct
   *        answer has been found (potentially, forever).
   *        The default value is 0
   * @param accounter accounts for the computation progress
   * @return an object of the type `InvariantValidationResult`
   */
  InvariantValidationResult
  check_invariant(const SetsUnion<Bundle> &init_set,
                  const SetsUnion<Polytope> &pSet,
                  const LinearSystem &invariant_candidate,
                  const unsigned int epoch_horizon = 0,
                  ProgressAccounter *accounter = NULL) const;
};

/**
 * Test whether all the sets in a list are empty
 *
 * @param[in] sets is a list of sets
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
