/*
 * Model.h
 *
 *  Created on: Jan 7, 2015
 *      Author: Tommaso Dreossi
 */

#ifndef MODEL_H_
#define MODEL_H_

#include "Common.h"
#include "Polyhedron.h"
#include "LinearSystemSet.h"
#include "STL/STL.h"
#include "STL/Atom.h"
#include "STL/Conjunction.h"
#include "STL/Disjunction.h"
#include "STL/Until.h"
#include "STL/Always.h"
#include "STL/Eventually.h"

#include <string>

class Model {

protected:
	lst vars; //list of variables
	lst params; // list of parameters
	lst dyns; // list of dynamics

	vector< Polyhedron* > reach_sets;			// Reachability sets
	vector<poly_values> initial_conditions;
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

	virtual vector<poly_values> getInitConds(){ return this->initial_conditions; }
	virtual vector<Polyhedron*> getReachSets(){ return this->reach_sets; }
	virtual LinearSystemSet* getInitParaSet(){ return this->init_para_set; }

	virtual STL* getSpec(){ return this->spec; }

	pair< vector< vector<double> >, vector< double > > normalizeVectors(vector<double> q, vector< vector<double> > vectors);

	ex genPoly( ex var, int deg );
	ex genPoly( lst vars, int deg );

};

#endif /* MODEL_H_ */
