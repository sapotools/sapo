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
#include "LinearSystem.h"
#include "LinearSystemSet.h"
#include "Bundle.h"
#include "Model.h"

class Sapo {

private:
	lst dyns;			// dynamics of the system
	lst vars;			// variables of the system
	lst params;			// parameters of the system
	sapo_opt options;	// options
	map< vector<int>,pair<lst,lst> > reachControlPts;		// symbolic control points
	map< vector<int>,pair<lst,lst> > synthControlPts;		// symbolic control points

	vector<Bundle*> reachWitDec(Bundle* initSet, int k);	// reachability with template decomposition
	LinearSystemSet* synthesizeUntil(Bundle *reachSet, LinearSystemSet *parameterSet, STL *sigma);
	LinearSystemSet* synthesizeAlways(Bundle *reachSet, LinearSystemSet *parameterSet, STL *formula);

public:
	Sapo(Model *model, sapo_opt options);

	vector<Bundle*> reach(Bundle* initSet, int k);							// reachability
	vector<Bundle*> reach(Bundle* initSet, LinearSystem* paraSet, int k);	// parameteric reachability

	// synthesis
	LinearSystemSet* synthesize(Bundle *reachSet, LinearSystemSet *parameterSet, STL *formula);
	LinearSystemSet* refineParameters(Bundle *reachSet, LinearSystemSet *parameterSet, STL *sigma);

	virtual ~Sapo();
};

#endif /* SAPO_H_ */
