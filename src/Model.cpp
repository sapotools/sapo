/*
 * Model.cpp
 *
 *  Created on: Oct 15, 2015
 *      Author: dreossi
 */

#include "Model.h"

// constructor for non-parametric dynamical systems
Model::Model(lst vars, lst dyns) {

	if( vars.nops() != dyns.nops() ){
		cout<<"Model::Model : vars and dyns must have the same size";
		exit (EXIT_FAILURE);
	}

	this->vars = vars;
	this->dyns = dyns;
}

// constructor for parametric dynamical systems
Model::Model(lst vars, lst params, lst dyns) {

	if( vars.nops() != dyns.nops() ){
		cout<<"Model::Model : vars and dyns must have the same size";
		exit (EXIT_FAILURE);
	}

	this->vars = vars;
	this->params = params;
	this->dyns = dyns;
}


Model::~Model() {
	// TODO Auto-generated destructor stub
}

