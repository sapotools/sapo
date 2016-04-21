/**
 * @file Model.h
 * Represent a discrete-time (eventually parameteric) dynamical system
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#ifndef MODEL_H_
#define MODEL_H_

#include "Common.h"

class Model {

private:
	lst vars;		// variables
	lst params;		// parameters
	lst dyns;		// dynamics

public:
	Model(lst vars, lst dyns);
	Model(lst vars, lst params, lst dyns);

	lst getVars(){ return this->vars; }
	lst getParams(){ return this->params; }
	lst getDyns(){ return this->dyns; }

	virtual ~Model();
};

#endif /* MODEL_H_ */
