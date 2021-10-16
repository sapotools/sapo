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

#include "Common.h"
#include "BaseConverter.h"
#include "STL.h"
#include "Atom.h"
#include "Conjunction.h"
#include "Disjunction.h"
#include "Always.h"
#include "Eventually.h"
#include "Until.h"
#include "LinearSystem.h"
#include "LinearSystemSet.h"
#include "Bundle.h"
#include "Model.h"
#include "Flowpipe.h"

class Sapo {

private:
	const lst& dyns;			// dynamics of the system
	const lst& vars;			// variables of the system
	const lst& params;			// parameters of the system
	const sapo_opt options;	// options
	map< vector<int>,pair<lst,lst> > reachControlPts;		// symbolic control points
	map< vector<int>,pair<lst,lst> > synthControlPts;		// symbolic control points

	vector<Bundle*> reachWitDec(Bundle* initSet, int k);	// reachability with template decomposition
	LinearSystemSet* synthesizeSTL(Bundle *reachSet, LinearSystemSet *parameterSet, const std::shared_ptr<STL> formula);
	LinearSystemSet* refineParameters(Bundle *reachSet, LinearSystemSet *parameterSet, const std::shared_ptr<Atom> formula);
	LinearSystemSet* synthesizeUntil(Bundle *reachSet, LinearSystemSet *parameterSet,
									 const std::shared_ptr<Until> formula, const int time);
	LinearSystemSet* synthesizeAlways(Bundle *reachSet, LinearSystemSet *parameterSet,
									  const std::shared_ptr<Always> formula, const int time);

public:
	Sapo(Model *model, sapo_opt options);

	Flowpipe* reach(Bundle* initSet, int k);							// reachability
	Flowpipe* reach(Bundle* initSet, LinearSystem* paraSet, int k);		// parameteric reachability
	LinearSystemSet* synthesize(Bundle *reachSet, LinearSystemSet *parameterSet, 
	                            const std::shared_ptr<STL> formula, const unsigned int max_splits=4);	// parameter synthesis

	virtual ~Sapo();
};

#endif /* SAPO_H_ */
