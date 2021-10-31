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

  std::vector<Bundle *> reachWitDec(Bundle &initSet, int k); // reachability with template decomposition
  LinearSystemSet synthesize_unpack(Bundle &reachSet, LinearSystemSet &parameterSet,
                                    const std::shared_ptr<STL> formula);
  LinearSystemSet synthesizeSTL(Bundle &reachSet,
                                LinearSystem &parameterSet,
                                const std::shared_ptr<STL> formula);
  LinearSystemSet refineParameters(Bundle &reachSet,
                                   LinearSystem &parameterSet,
                                   const std::shared_ptr<Atom> formula);
  LinearSystemSet synthesizeUntil(Bundle &reachSet,
                                  LinearSystem &parameterSet,
                                  const std::shared_ptr<Until> formula,
                                  const int time);
  LinearSystemSet synthesizeAlways(Bundle &reachSet,
                                   LinearSystem &parameterSet,
                                   const std::shared_ptr<Always> formula,
                                   const int time);

public:
  Sapo(Model *model, sapo_opt options);

  Flowpipe reach(const Bundle &initSet, unsigned int k); // reachability
  Flowpipe reach(const Bundle &initSet, LinearSystemSet &paraSet,
                 unsigned int k); // parameteric reachability
  LinearSystemSet synthesize(Bundle &reachSet, LinearSystemSet &parameterSet,
                             const std::shared_ptr<STL> formula,
                             const unsigned int max_splits
                             = 4); // parameter synthesis

  virtual ~Sapo();
};

#endif /* SAPO_H_ */
