/*
 * DiscreteDynamicalSystem.h
 *
 *  Created on: Nov 5, 2014
 *      Author: dreossi
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

	Polyhedron* reachSet;				// reachability set
	lst fog;							// dynamics composed generator function
	vector< lst > templateControlPts;	// template of control points
public:
	DiscreteDynamicalSystem(lst vars, lst params, lst dynamics, bool rational);
	DiscreteDynamicalSystem(lst vars, lst params, lst dynamics, Polyhedron *reachSet, bool rational);

	lst getVars();
	lst getParams();
	lst getDynamics();
	bool isRational(){return this->rational;};

	Polyhedron* getReachSet();
	lst getFog();
	vector< lst > getTemplateControlPts();
	virtual ~DiscreteDynamicalSystem();

};

#endif /* DISCRETEDYNAMICALSYSTEM_H_ */
