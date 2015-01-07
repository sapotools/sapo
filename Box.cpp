/*
 * Box.cpp
 *
 *  Created on: Nov 3, 2014
 *      Author: dreossi
 */

#include "Box.h"

Box::Box(vector<lst> vars) {

	if(vars.size() != 3){
		cout<<"Box::Box : vars must contain 3 collections of variable names (q,alpha,beta)";
		exit (EXIT_FAILURE);
	}

	this->vars.push_back(vars[0]);
	this->vars.push_back(vars[1]);
	this->vars.push_back(vars[2]);

	// get the dimension of the box
	this->dim = vars[0].nops();

	// and store its variable names
	for(int i=0; i<3; i++){
		if((signed)vars[i].nops() != this->dim){
			cout<<"Box::Box : vars["<<i<<"] must have "<<this->dim<<" variables";
			exit (EXIT_FAILURE);
		}
	}

	// extract variable names
	lst q = this->vars[0];
	lst alpha = this->vars[1];
	lst beta = this->vars[2];

	// initialize generator function
	for(int i=0; i<this->dim; i++){
		this->generator_function.append(q[i]);
	}

	// create the generation function acclumulatin the versor values
	for(int i=0; i<this->dim; i++){
		this->generator_function[i] = this->generator_function[i] + alpha[i]*beta[i];
	}

}

// convert from generator to constraints representation
// q : numeric base vertex
LinearSystem* Box::gen2const(vector<double> q, vector<double> beta){

	if((signed)q.size() != this->dim){
		cout<<"Box::gen2const : q must have dimension "<<this->dim;
		exit (EXIT_FAILURE);
	}

	vector< vector<double> > Lambda;
	vector< double > d;

	// initialize template Lambda
	for(int i=0; i<this->dim; i++){
		vector<double> Lambda_i (this->dim,0);
		Lambda_i[i] = 1;
		Lambda.push_back(Lambda_i);
		d.push_back(q[i]+beta[i]);
	}
	for(int i=0; i<this->dim; i++){
		vector<double> Lambda_i (this->dim,0);
		Lambda_i[i] = -1;
		Lambda.push_back(Lambda_i);
		d.push_back(-q[i]);
	}

	LinearSystem *LS = new LinearSystem(Lambda,d);
	return LS;
}

poly_values Box::const2gen(LinearSystem *constr){

	poly_values res;
	return res;

}


Box::~Box() {
	// TODO Auto-generated destructor stub
}

