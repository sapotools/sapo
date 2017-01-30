/**
 * @file Model.cpp
 * Represent a discrete-time (eventually parameteric) dynamical system
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#include "Model.h"

/**
 * Constructor that instantiates a non-parameteric dynamical system
 *
 * @param[in] vars list of variables
 * @param[in] dyns list of dynamics
 */

Model::Model() { }


/*Model::Model(lst vars, lst dyns) {

	if( vars.nops() != dyns.nops() ){
		cout<<"Model::Model : vars and dyns must have the same size";
		exit (EXIT_FAILURE);
	}

	this->vars = vars;
	this->dyns = dyns;
}*/

/**
 * Constructor that instantiates a parameteric dynamical system
 *
 * @param[in] vars list of variables
 * @param[in] params list of parameters
 * @param[in] dyns list of dynamics
 */
/*Model::Model(lst vars, lst params, lst dyns) {

	if( vars.nops() != dyns.nops() ){
		cout<<"Model::Model : vars and dyns must have the same size";
		exit (EXIT_FAILURE);
	}

	this->vars = vars;
	this->params = params;
	this->dyns = dyns;
}*/


Model::~Model() {
	// TODO Auto-generated destructor stub
}
