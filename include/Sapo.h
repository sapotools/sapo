/*
 * Sapo.h
 *
 *  Created on: Oct 12, 2015
 *      Author: dreossi
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

class Sapo {

private:
	lst dyns;	// dynamics of the system
	lst vars;	// variables of the system
	lst params;	// parameters of the system
	sapo_opt options;	// options
	map< vector<int>,pair<lst,lst> > reachControlPts;	// symbolic control points
	map<vector<int>,lst> synthControlPts;	// symbolic control points


	vector<Bundle*> reachWitDec(Bundle* initSet, int k);

public:
	Sapo(lst vars, lst dyns, sapo_opt options);
	Sapo(lst vars, lst params, lst dyns, sapo_opt options);

	vector<Bundle*> reach(Bundle* initSet, int k);
	vector<Bundle*> reach(Bundle* initSet, LinearSystem* paraSet, int k);

	// synthesis
	LinearSystemSet* synthesize(Bundle *reachSet, LinearSystemSet *parameterSet, STL *formula);
	LinearSystemSet* refineParameters(Bundle *reachSet, LinearSystemSet *parameterSet, STL *sigma);
	LinearSystemSet* synthesizeUntil(Bundle *reachSet, LinearSystemSet *parameterSet, STL *sigma);
	LinearSystemSet* synthesizeAlways(Bundle *reachSet, LinearSystemSet *parameterSet, STL *formula);


	virtual ~Sapo();
};

#endif /* SAPO_H_ */
