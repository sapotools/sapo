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
#include "PolytopesUnion.h"
#include "Model.h"
#include "STL.h"
#include "Until.h"

#include "ProgressAccounter.h"

class Sapo
{
public:
  Bundle::transfomation_mode tmode; //!< transformation mode (OFO or AFO)
  double decomp_weight;             //!< decomposition weight
  unsigned int decomp; //!< number of decompositions (0: none, >0: yes)
  std::string plot;    //!< the name of the file were to plot the reach set
  unsigned int time_horizon;     //!< the computation time horizon
  unsigned int max_param_splits; //!< maximum number of splits in synthesis
  unsigned int num_of_presplits; //!< number of presplits in synthesis
  double max_bundle_magnitude; //!< maximum versor magnitude for single bundle

private:
  const std::vector<SymbolicAlgebra::Expression<>>
      &dyns; //!< dynamics of the system
  const std::vector<SymbolicAlgebra::Symbol<>>
      &vars; //!< variables of the system
  const std::vector<SymbolicAlgebra::Symbol<>>
      &params; //!< parameters of the system

  // TODO: check whether the following method is really needed/usable.
  std::vector<Bundle *>
  reachWitDec(Bundle &initSet,
              int k); // reachability with template decomposition

  /**
   * Parameter synthesis w.r.t. an always formula
   *
   * @param[in] reachSet bundle with the initial set
   * @param[in] parameterSet set of parameters
   * @param[in] sigma STL always formula
   * @returns refined sets of parameters
   */
  template<typename T>
  PolytopesUnion transition_and_synthesis(const Bundle &reachSet,
                                          const PolytopesUnion &pSet,
                                          const std::shared_ptr<T> formula,
                                          const int time) const
  {
    PolytopesUnion result;

    for (auto p_it = pSet.begin(); p_it != pSet.end(); ++p_it) {
      // transition by using the n-th polytope of the parameter set
      Bundle newReachSet = reachSet.transform(this->vars, this->params,
                                              this->dyns, *p_it, this->tmode);

      result.add(synthesize(newReachSet, pSet, formula, time + 1));
    }

    return result;
  }

  /**
   * Parameter synthesis for atomic formulas
   *
   * @param[in] reachSet bundle representing the current reached set
   * @param[in] pSet the current set of parameters
   * @param[in] formula is an STL atomic formula providing the specification
   * @returns refined parameter set
   */
  PolytopesUnion synthesize(const Bundle &reachSet, const PolytopesUnion &pSet,
                            const std::shared_ptr<STL::Atom> formula) const;

  /**
   * Parmeter synthesis for conjunctions
   *
   * @param[in] reachSet bundle representing the current reached set
   * @param[in] pSet the current parameter set
   * @param[in] conj is an STL conjunction providing the specification
   * @returns refined parameter set
   */
  PolytopesUnion synthesize(const Bundle &reachSet, const PolytopesUnion &pSet,
                            const std::shared_ptr<STL::Conjunction> formula) const;

  /**
   * Parmeter synthesis for disjunctions
   *
   * @param[in] reachSet bundle representing the current reached set
   * @param[in] pSet the current parameter set
   * @param[in] conj is an STL disjunction providing the specification
   * @returns refined parameter set
   */
  PolytopesUnion synthesize(const Bundle &reachSet, const PolytopesUnion &pSet,
                            const std::shared_ptr<STL::Disjunction> formula) const;

  /**
   * Parameter synthesis for until formulas
   *
   * @param[in] reachSet bundle representing the current reached set
   * @param[in] pSet the current set of parameters
   * @param[in] formula is an STL until formula providing the specification
   * @param[in] time is the time of the current evaluation
   * @returns refined parameter set
   */
  PolytopesUnion synthesize(const Bundle &reachSet, const PolytopesUnion &pSet,
                            const std::shared_ptr<STL::Until> formula,
                            const int time) const;

  /**
   * Parameter synthesis for always formulas
   *
   * @param[in] reachSet bundle representing the current reached set
   * @param[in] pSet the current set of parameters
   * @param[in] formula is an STL always formula providing the specification
   * @param[in] time is the time of the current evaluation
   * @returns refined parameter set
   */
  PolytopesUnion synthesize(const Bundle &reachSet, const PolytopesUnion &pSet,
                            const std::shared_ptr<STL::Always> formula,
                            const int time) const;
  /**
   * Parmeter synthesis for the eventually fomulas
   *
   * @param[in] reachSet bundle representing the current reached set
   * @param[in] pSet the current parameter set
   * @param[in] formula is an STL eventually formula providing the
   * specification
   * @returns refined parameter set
   */
  PolytopesUnion synthesize(const Bundle &reachSet, const PolytopesUnion &pSet,
                            const std::shared_ptr<STL::Eventually> ev) const;

public:
  /**
   * Constructor that instantiates Sapo
   *
   * @param[in] model model to analyize
   * @param[in] sapo_opt options to tune sapo
   */
  Sapo(Model *model);

  /**
   * Reachable set computation
   *
   * @param[in] initSet bundle representing the current reached set
   * @param[in] k time horizon
   * @param[in,out] accounter acccounts for the computation progress
   * @returns the reached flowpipe
   */
  Flowpipe reach(const Bundle &initSet, unsigned int k,
                 ProgressAccounter *accounter = NULL) const;

  /**
   * Reachable set computation for parameteric dynamical systems
   *
   * @param[in] initSet bundle representing the current reached set
   * @param[in] pSet the set of parameters
   * @param[in] k time horizon
   * @param[in,out] accounter acccounts for the computation progress
   * @returns the reached flowpipe
   */
  Flowpipe reach(const Bundle &initSet, const PolytopesUnion &pSet,
                 unsigned int k, ProgressAccounter *accounter = NULL) const;

  /**
   * Parameter synthesis method
   *
   * @param[in] reachSet bundle representing the current reached set
   * @param[in] pSet the current parameter set
   * @param[in] formula is an STL specification for the model
   * @param[in,out] accounter acccounts for the computation progress
   * @returns refined parameter set
   */
  PolytopesUnion synthesize(const Bundle &reachSet, const PolytopesUnion &pSet,
                            const std::shared_ptr<STL::STL> formula,
                            ProgressAccounter *accounter = NULL) const;

  /**
   * Parameter synthesis with splits
   *
   * @param[in] reachSet bundle representing the current reached set
   * @param[in] pSet is the current parameter sets
   * @param[in] formula is an STL formula providing the specification
   * @param[in] max_splits maximum number of splits of the original
   *                       parameter set to identify a non-null solution
   * @param[in] num_of_presplits is number of splits to be performed before
   *                             the computation
   * @param[in,out] accounter acccounts for the computation progress
   * @returns the list of refined parameter sets
   */
  std::list<PolytopesUnion>
  synthesize(const Bundle &reachSet, const PolytopesUnion &pSet,
             const std::shared_ptr<STL::STL> formula, const unsigned int max_splits,
             const unsigned int num_of_presplits = 0,
             ProgressAccounter *accounter = NULL) const;
};

#endif /* SAPO_H_ */
