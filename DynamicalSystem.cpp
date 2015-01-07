/*
 * DynamicalSystem.cpp
 *
 *  Created on: Oct 30, 2014
 *      Author: dreossi
 */

#include "DynamicalSystem.h"

DynamicalSystem::DynamicalSystem(lst vars, lst params, lst dynamics) {

	this->vars = vars;
	this->params = params;
	this->dynamics = dynamics;

}

// Discretiazion of the dynamics with Euler's method
lst DynamicalSystem::eulerDisc(double disc_step){

	lst disc_dynamics;

	for(int i=0; i<this->dynamics.nops(); i++){
		disc_dynamics.append( this->vars[i] + disc_step*this->dynamics[i] );
	}
	return disc_dynamics;

}



lst DynamicalSystem::getVars(){ return this->vars; }
lst DynamicalSystem::getParams(){ return this->params; }
lst DynamicalSystem::getDynamics(){ return this->dynamics; }

DynamicalSystem::~DynamicalSystem() {
	// TODO Auto-generated destructor stub
}

