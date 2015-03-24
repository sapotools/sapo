/*
 * DynamicalSystem.h
 *
 *  Created on: Oct 30, 2014
 *      Author: dreossi
 */

#ifndef DYNAMICALSYSTEM_H_
#define DYNAMICALSYSTEM_H_

#include "Common.h"

#include "Polyhedron.h"
#include "LinearSystem.h"

class DynamicalSystem {

private:

	lst vars;				// variables
	lst params;				// parameters
	lst dynamics;			// dynamics
	bool rational;			// has rational dynamics?

public:

	DynamicalSystem(lst vars, lst params, lst dynamics, bool rational);

	lst eulerDisc(double disc_step);

	lst getVars();
	lst getParams();
	lst getDynamics();
	bool isRational(){ return this->rational; };

	virtual ~DynamicalSystem();
};

#endif /* DYNAMICALSYSTEM_H_ */
