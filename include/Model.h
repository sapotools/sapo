/**
 * @file Model.h
 * Represent a discrete-time (eventually parameteric) dynamical system
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#ifndef MODEL_H_
#define MODEL_H_

#include "Bundle.h"
#include "LinearSystemSet.h"
#include "STL.h"

#include "Common.h"

class Model {

protected:
	lst vars;		// variables
	lst params;		// parameters
	lst dyns;		// dynamics

	Bundle *reachSet; // Initial reach set
	LinearSystemSet *paraSet;

	STL *spec;

public:

	char name[64];

	char* getName(){ return this->name; }
	lst getVars(){ return this->vars; }
	lst getParams(){ return this->params; }
	lst getDyns(){ return this->dyns; }

	Bundle* getReachSet(){ return this->reachSet; }
	LinearSystemSet* getParaSet(){ return this->paraSet; }
	STL* getSpec(){ return this->spec; }

};

#endif /* MODEL_H_ */
