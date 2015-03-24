/*
 * Model.h
 *
 *  Created on: Jan 7, 2015
 *      Author: dreossi
 */

#ifndef MODEL_H_
#define MODEL_H_

#include "Common.h"
#include "Polyhedron.h"
#include "LinearSystemSet.h"
#include "STL.h"
#include "Atom.h"
#include "Conjunction.h"
#include "Disjunction.h"
#include "Until.h"
#include "Always.h"
#include "Eventually.h"

#include <string>

class Model {

protected:
	lst vars; //list of variables
	lst params; // list of parameters
	lst dyns; // list of dynamics

	Polyhedron *reach_set;			// Reachability set
	LinearSystemSet *init_para_set;	// Parameter set

	STL *spec;		// Stl specification

	string name;

public:
	Model();
	virtual ~Model();

	virtual lst getVars(){ return this->vars; }
	virtual lst getParams(){ return this->params; }
	virtual lst getDyns(){ return this->dyns; }
	virtual string getName(){ return this->name; }

	virtual Polyhedron* getReachSet(){ return this->reach_set; }
	virtual LinearSystemSet* getInitParaSet(){ return this->init_para_set; }

	virtual STL* getSpec(){ return this->spec; }
};

#endif /* MODEL_H_ */
