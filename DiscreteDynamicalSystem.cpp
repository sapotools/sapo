/*
 * DiscreteDynamicalSystem.cpp
 *
 *  Created on: Nov 5, 2014
 *      Author: dreossi
 */

#include "DiscreteDynamicalSystem.h"

DiscreteDynamicalSystem::DiscreteDynamicalSystem(lst vars, lst params, lst dynamics) {

	this->vars = vars;
	this->params = params;
	this->dynamics = dynamics;

}

DiscreteDynamicalSystem::DiscreteDynamicalSystem(lst vars, lst params, lst dynamics, Polyhedron *reachSet) {

	this->vars = vars;
	this->params = params;
	this->dynamics = dynamics;
	this->reachSet = reachSet;

	// Extract the template matrix Lambda from reachSet
	int reachSetDim = this->reachSet->getDim();
	vector<double> base_vertex (reachSetDim, 0);
	vector<double> betas (reachSetDim, 1);
	LinearSystem *LS = this->reachSet->gen2const(base_vertex, betas);
	vector< vector<double> > Lambda = LS->getA();

	// Prepare the substitutions for the function composition
	lst generator_function = this->reachSet->getGeneratorFunction();
	lst sub, fog;
	for(int i=0; i<(signed)generator_function.nops(); i++){
		sub.append(vars[i] == generator_function[i]);
	}

	// Compute fog
	for(int i=0; i<(signed)dynamics.nops(); i++){
		fog.append(dynamics[i].subs(sub));
	}

	this->fog = fog; // store the composition

	cout<<"Computing template...\n";
	vector< ex > Lambda_fog;
	// Compute Lambda_i(fog)
	for( int i=0; i<(signed)Lambda.size(); i++ ){
		ex Lambda_fog_i = 0;
		for( int j=0; j<(signed)Lambda[i].size(); j++ ){
			Lambda_fog_i = Lambda_fog_i + Lambda[i][j]*fog[j];
		}
		Lambda_fog.push_back(Lambda_fog_i);
		//cout<<Lambda_fog[i]<<"\n";
	}

	cout<<"Computing control points...\n";
	// Compute the template control points
	vector< lst > templateControlPts;
	for( int i=0; i<(signed)Lambda_fog.size(); i++){
		BaseConverter *BC = new BaseConverter(this->reachSet->getAlpha(), Lambda_fog[i]);
		templateControlPts.push_back(BC->getBernCoeffs());								// compute control points
	}
	this->templateControlPts = templateControlPts;

}

lst DiscreteDynamicalSystem::getVars(){ return this->vars; }
lst DiscreteDynamicalSystem::getParams(){ return this->params; }
lst DiscreteDynamicalSystem::getDynamics(){ return this->dynamics; }

Polyhedron* DiscreteDynamicalSystem::getReachSet(){ return this->reachSet; }
lst DiscreteDynamicalSystem::getFog(){ return this->fog; }
vector< lst > DiscreteDynamicalSystem::getTemplateControlPts(){ return this->templateControlPts; }

DiscreteDynamicalSystem::~DiscreteDynamicalSystem() {
	// TODO Auto-generated destructor stub
}

