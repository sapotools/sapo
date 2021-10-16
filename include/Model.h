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

	std::shared_ptr<STL> spec;

	std::string name;
public:

	inline const std::string& getName() const { return this->name; }
	inline const lst& getVars() const { return this->vars; }
	inline const lst& getParams() const { return this->params; }
	inline const lst& getDyns() const { return this->dyns; }

	inline Bundle* getReachSet(){ return this->reachSet; }
	inline LinearSystemSet* getParaSet(){ return this->paraSet; }
	inline const std::shared_ptr<STL>  getSpec() const { return this->spec; }

	~Model();
};

#endif /* MODEL_H_ */
