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

#include "Always.h"
#include "Atom.h"
#include "BaseConverter.h"
#include "Bundle.h"
#include "Common.h"
#include "Conjunction.h"
#include "ControlPointStorage.h"
#include "Disjunction.h"
#include "Eventually.h"
#include "Flowpipe.h"
#include "LinearSystem.h"
#include "LinearSystemSet.h"
#include "Model.h"
#include "STL.h"
#include "Until.h"
#include "ControlPointStorage.h"

class Sapo
{

private:
  const GiNaC::lst &dyns;   // dynamics of the system
  const GiNaC::lst &vars;   // variables of the system
  const GiNaC::lst &params; // parameters of the system
  const sapo_opt options;   // options

  // TODO: check whether the following two members are really Sapo properties
  ControlPointStorage reachControlPts; // symbolic control points
  ControlPointStorage synthControlPts; // symbolic control points

  // TODO: check whether the following method is really needed/usable.
  std::vector<Bundle *> reachWitDec(Bundle &initSet, int k); // reachability with template decomposition

  /**
   * Parameter synthesis w.r.t. an always formula
   *
   * @param[in] reachSet bundle with the initial set
   * @param[in] parameterSet set of parameters
   * @param[in] sigma STL always formula
   * @returns refined sets of parameters
   */
  template <typename T>
  LinearSystemSet transition_and_synthesis(Bundle &reachSet, 
                                const LinearSystemSet &pSet,
                                const std::shared_ptr<T> formula,
                                const int time)
  {
    LinearSystemSet result;

    for (auto p_it = pSet.begin(); p_it != pSet.end(); ++p_it) {
      // transition by using the n-th linear system of the parameter set
      Bundle newReachSet
          = reachSet.transform(this->vars, this->params, this->dyns, *p_it,
                                this->reachControlPts, this->options.trans);

      // TODO: Check whether the object tmpLSset can be removed
      result.unionWith(
          synthesize(newReachSet, pSet, formula, time + 1));
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
  LinearSystemSet synthesize(Bundle &reachSet,
                                   const LinearSystemSet &pSet,
                                   const std::shared_ptr<Atom> formula);

  /**
   * Parmeter synthesis for conjunctions
   *
   * @param[in] reachSet bundle representing the current reached set
   * @param[in] pSet the current parameter set
   * @param[in] conj is an STL conjunction providing the specification
   * @returns refined parameter set
   */
  LinearSystemSet synthesize(Bundle &reachSet,
                                  const LinearSystemSet &pSet,
                                  const std::shared_ptr<Conjunction> formula);

  /**
   * Parmeter synthesis for disjunctions
   *
   * @param[in] reachSet bundle representing the current reached set
   * @param[in] pSet the current parameter set
   * @param[in] conj is an STL disjunction providing the specification
   * @returns refined parameter set
   */
  LinearSystemSet synthesize(Bundle &reachSet,
                                  const LinearSystemSet &pSet,
                                  const std::shared_ptr<Disjunction> formula);

  /**
   * Parameter synthesis for until formulas
   *
   * @param[in] reachSet bundle representing the current reached set
   * @param[in] pSet the current set of parameters
   * @param[in] formula is an STL until formula providing the specification
   * @param[in] time is the time of the current evaluation
   * @returns refined parameter set
   */
  LinearSystemSet synthesize(Bundle &reachSet,
                                  const LinearSystemSet &pSet,
                                  const std::shared_ptr<Until> formula,
                                  const int time);

  /**
   * Parameter synthesis for always formulas
   *
   * @param[in] reachSet bundle representing the current reached set
   * @param[in] pSet the current set of parameters
   * @param[in] formula is an STL always formula providing the specification
   * @param[in] time is the time of the current evaluation
   * @returns refined parameter set
   */
  LinearSystemSet synthesize(Bundle &reachSet,
                                   const LinearSystemSet &pSet,
                                   const std::shared_ptr<Always> formula,
                                   const int time);
  /**
   * Parmeter synthesis for the eventually fomulas
   *
   * @param[in] reachSet bundle representing the current reached set
   * @param[in] pSet the current parameter set
   * @param[in] formula is an STL eventually formula providing the specification
   * @returns refined parameter set
   */
  LinearSystemSet synthesize(Bundle &reachSet,
                                  const LinearSystemSet &pSet,
                                  const std::shared_ptr<Eventually> ev);
public:
  /**
   * Constructor that instantiates Sapo
   *
   * @param[in] model model to analyize
   * @param[in] sapo_opt options to tune sapo
   */
  Sapo(Model *model, sapo_opt options);

  /**
   * Reachable set computation
   *
   * @param[in] initSet bundle representing the current reached set
   * @param[in] k time horizon
   * @returns the reached flowpipe
   */
  Flowpipe reach(const Bundle &initSet, unsigned int k);

  /**
   * Reachable set computation for parameteric dynamical systems
   *
   * @param[in] initSet bundle representing the current reached set
   * @param[in] pSet the set of parameters
   * @param[in] k time horizon
   * @returns the reached flowpipe
   */
  Flowpipe reach(const Bundle &initSet, const LinearSystemSet &pSet,
                 unsigned int k);

  /**
   * Parameter synthesis method
   *
   * @param[in] reachSet bundle representing the current reached set
   * @param[in] pSet the current parameter set
   * @param[in] formula is an STL specification for the model
   * @returns refined parameter set
   */
  LinearSystemSet synthesize(Bundle &reachSet, const LinearSystemSet &pSet,
                             const std::shared_ptr<STL> formula);
  
  /**
   * Parameter synthesis with splits
   *
   * @param[in] reachSet bundle representing the current reached set
   * @param[in] pSet is the current parameter sets
   * @param[in] formula is an STL formula providing the specification
   * @param[in] max_splits maximum number of splits of the original
   *                       parameter set to identify a non-null solution
   * @returns the list of refined parameter sets
   */
  std::list<LinearSystemSet> synthesize(Bundle &reachSet, const LinearSystemSet& pSet,
                             const std::shared_ptr<STL> formula,
                             const unsigned int max_splits); // parameter synthesis
  virtual ~Sapo();
};

#endif /* SAPO_H_ */
