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

#include "Bernstein.h"
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
  enum {
    DISPROVED = 0,
    PROVED = 1,
    EPOCH_LIMIT_REACHED = 2
  } result_type;              //!< result type
  Flowpipe flowpipe;          //!< system flowpipe
  Flowpipe k_induction_proof; //!< flowpipe of k-induction proof
} InvariantValidationResult;

/**
 * @brief The computation orchestrator
 */
class Sapo
{

public:
  std::string plot; //!< the name of the file were to plot the reach set
  unsigned int time_horizon;      //!< the computation time horizon
  unsigned int max_param_splits;  //!< maximum number of splits in synthesis
  unsigned int num_of_pre_splits; //!< number of pre-splits in synthesis
  double max_bundle_magnitude; //!< maximum versor magnitude for single bundle

  // invariant fields

  unsigned int max_k_induction;     //!< maximum k to be tested for k-induction
  double delta_thickness_threshold; //!< the minimum thickness difference
                                    //!< between evolution steps
  unsigned int
      missed_thickness_threshold; //!< how many times the thickness threshold
                                  //!< must be missed before an expansion

  /**
   * @brief Approximation used during invariant validation
   */
  enum joinApproxType {
    NO_APPROX,  // no approximation is used (listing)
    CHAIN_JOIN, // chain-join approximation (merging)
    FULL_JOIN   // full-join approximation (packaging)
  };

  joinApproxType join_approx; //!< join approximation type for `k`-induction
                              //!< invariant proof

private:
  Evolver<double> *_evolver; //!< the dynamical system evolver

  const LinearSystem<double> assumptions;

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
                           const int epoch_horizon)
  {
    SetsUnion<Polytope> result;

    // init_set.intersect(this->assumptions);

    for (auto p_it = pSet.begin(); p_it != pSet.end(); ++p_it) {
      // transition by using the n-th polytope of the parameter set
      Bundle reached_set = _evolver->operator()(init_set, *p_it);

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
  SetsUnion<Polytope> synthesize(const Bundle &reachSet,
                                 const SetsUnion<Polytope> &pSet,
                                 const std::shared_ptr<STL::Atom> formula);

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
             const std::shared_ptr<STL::Conjunction> formula);

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
             const std::shared_ptr<STL::Disjunction> formula);

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
                                 const int time);

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
                                 const int time);
  /**
   * Parameter synthesis for the eventually formulas
   *
   * @param[in] init_set is the initial set
   * @param[in] pSet the current parameter set
   * @param[in] formula is an STL eventually formula providing the
   * specification
   * @returns refined parameter set
   */
  SetsUnion<Polytope> synthesize(const Bundle &reachSet,
                                 const SetsUnion<Polytope> &pSet,
                                 const std::shared_ptr<STL::Eventually> ev);

public:
  /**
   * Constructor
   *
   * @param[in] model is the model to analyze
   * @param[in] cached is a flag to cache Bernstein coefficients
   */
  Sapo(const DiscreteModel &model, bool cached = true);

  /**
   * @brief Get the dynamical system
   *
   * @return the dynamical system
   */
  inline const DynamicalSystem<double> &dynamical_system() const
  {
    return _evolver->dynamical_system();
  }

  /**
   * @brief Get the dynamical system evolver
   *
   * @return the dynamical system evolver
   */
  inline Evolver<double> *evolver()
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
    _evolver->mode = mode;
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
                 ProgressAccounter *accounter = NULL);

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
                 ProgressAccounter *accounter = NULL);

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
                                 ProgressAccounter *accounter = NULL);

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
             ProgressAccounter *accounter = NULL);

  /**
   * @brief Try to establish whether a set is an invariant
   *
   * @param init_set is the initial set
   * @param pSet is the parameter set
   * @param invariant_candidate is the candidate invariant
   * @param max_k_induction is the maximum `k` for
   *        `k`-induction. When set to 0, `k` is not
   *        upper-bounded and grows as the epoch increases
   * @param accounter accounts for the computation progress
   * @return a pair<bool, unsigned>. The first value
   * is `true` if and only if the method has identified
   * a `k` such that `invariant_candidate` is a `k`-invariant
   * for the investigated model. In such a case, the second
   * returned value is such a `k`. Whenever the first value
   * is `false`, then the second value report the epoch
   * in which `invariant_candidate` has been disproved to be
   * a `k`-invariant
   */
  InvariantValidationResult
  check_invariant(const SetsUnion<Bundle> &init_set,
                  const SetsUnion<Polytope> &pSet,
                  const LinearSystem<double> &invariant_candidate,
                  ProgressAccounter *accounter = NULL);

  /**
   * @brief Destroyer
   */
  ~Sapo();
};

/**
 * Test whether all the sets in a list are empty
 *
 * @param[in] sets is a list of sets
 * @returns `true` when the function establishes that either
 *          `sets` has null size or every set in `sets` is
 *          empty. `false` if the function establishes that
 *          `sets` contains one non-empty set. `uncertain` in
 *          the remaining cases
 */
template<class T>
TriBool every_set_is_empty(const std::list<T> &sets)
{
  TriBool res{true};

  for (auto s_it = std::begin(sets); s_it != std::end(sets); ++s_it) {

    res = res && s_it->is_empty();
    if (is_false(res)) {
      return res;
    }
  }

  return res;
}

#endif /* SAPO_H_ */
