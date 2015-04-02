/*
 * DiscreteDynamicalSystem.h
 *
 *  Created on: Nov 5, 2014
 *      Author: Tommaso Dreossi
 */

#ifndef DISCRETEDYNAMICALSYSTEM_H_
#define DISCRETEDYNAMICALSYSTEM_H_

#include "DynamicalSystem.h"
#include "BaseConverter.h"

class DiscreteDynamicalSystem {

private:
	lst vars;				// variables
	lst params;				// parameters
	lst dynamics;			// dynamics
	bool rational;			// has the system rational dynamics?

	// Multiple reach sets
	vector< Polyhedron* > reachSets;				// various templates
	vector< lst > fogs;								// dynamics composed generator functions
	vector< vector< lst > > templatesControlPts;	// templates of control points


public:
	DiscreteDynamicalSystem(lst vars, lst params, lst dynamics, bool rational);
	DiscreteDynamicalSystem(lst vars, lst params, lst dynamics, Polyhedron *reachSet, bool rational);
	DiscreteDynamicalSystem(lst vars, lst params, lst dynamics, vector< Polyhedron* > reachSet, bool rational);

	lst getVars();
	lst getParams();
	lst getDynamics();
	bool isRational(){return this->rational;};

	vector< Polyhedron* > getReachSets();
	Polyhedron* getReachSet(int i);
	vector< lst > getFogs();
	vector< vector< lst > > getTemplatesControlPts();
	virtual ~DiscreteDynamicalSystem();

};

#endif /* DISCRETEDYNAMICALSYSTEM_H_ */
