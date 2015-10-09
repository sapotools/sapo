/*
 * DiscreteDynamicalSystem.cpp
 *
 *  Created on: Nov 5, 2014
 *      Author: Tommaso Dreossi
 */

#include "DiscreteDynamicalSystem.h"

DiscreteDynamicalSystem::DiscreteDynamicalSystem(lst vars, lst params, lst dynamics, bool rational) {

	this->vars = vars;
	this->params = params;
	this->dynamics = dynamics;
	this->rational = rational;

}

DiscreteDynamicalSystem::DiscreteDynamicalSystem(lst vars, lst params, lst dynamics, Polyhedron *reachSet, bool rational) {

	vector< Polyhedron* > reachSets (1,reachSet);
	DiscreteDynamicalSystem(vars,params,dynamics, reachSets, rational);

}

DiscreteDynamicalSystem::DiscreteDynamicalSystem(lst vars, lst params, lst dynamics, vector<Polyhedron*> reachSets, bool rational) {

	this->vars = vars;
	this->params = params;
	this->dynamics = dynamics;
	this->reachSets = reachSets;
	this->rational = rational;

	// Extract the template matrix Lambda from reachSet
	int reachSetDim = this->reachSets[0]->getDim();
	vector<double> base_vertex (reachSetDim, 0);
	vector<double> betas (reachSetDim, 1);

	vector<LinearSystem *> LSs;
	vector< vector< vector<double> > > Lambdas;
	for(int i=0; i<(signed)this->reachSets.size(); i++){
		LSs.push_back(this->reachSets[i]->gen2const(base_vertex,betas));
		Lambdas.push_back(LSs[i]->getA());
	}

	// Prepare the substitutions for the function composition

	vector<lst> generator_functions;
	vector<lst> fogs;
	vector< vector< ex > > Lambda_fogs;

	for(int i=0; i<(signed)this->reachSets.size(); i++){

		generator_functions.push_back(this->reachSets[i]->getGeneratorFunction());
		//cout<<"\nDD gen fun: "<<generator_functions[i]<<"\n";

		lst sub;
		for(int j=0; j<(signed)generator_functions[i].nops(); j++){
			sub.append(vars[j] == generator_functions[i][j]);
		}

		// Compute fog
		lst fog;
		for(int j=0; j<(signed)dynamics.nops(); j++){
			fog.append(dynamics[j].subs(sub));
		}

		//cout<<"fog: "<<fog<<"\n";

		this->fogs.push_back(fog); // store the composition

		cout<<"Computing template (n. "<<i<<")...\n";
		vector< ex > Lambda_fog;
		// Compute Lambda_i(fog)
		for( int j=0; j<(signed)Lambdas[i].size(); j++ ){
			ex Lambda_fog_j = 0;
			for( int k=0; k<(signed)Lambdas[i][j].size(); k++ ){
				Lambda_fog_j = Lambda_fog_j + Lambdas[i][j][k]*fog[k];
			}
			Lambda_fog.push_back(Lambda_fog_j);
			//cout<<Lambda_fog[i]<<"\n";
		}
		Lambda_fogs.push_back(Lambda_fog);
	}


	for(int i=0; i<(signed)this->reachSets.size(); i++){

		cout<<"Computing control points (n. "<<i<<")...\n";
		// Compute the template control points
		vector< lst > templateControlPts;
		if(!this->rational){		// check if the dynamics are rational
			for( int j=0; j<(signed)Lambda_fogs[i].size(); j++){
				//cout<<Lambda_fogs[i][j]<<"\n";
				BaseConverter *BC = new BaseConverter(this->reachSets[i]->getAlpha(), Lambda_fogs[i][j]);
				templateControlPts.push_back(BC->getBernCoeffs());								// compute control points
			}
		}else{
			for( int j=0; j<(signed)Lambda_fogs[i].size(); j++){
				BaseConverter *BC = new BaseConverter(this->reachSets[i]->getAlpha(), Lambda_fogs[i][j].numer(),Lambda_fogs[i][j].denom());
				templateControlPts.push_back(BC->getRationalBernCoeffs());								// compute control points
			}

		}
		this->templatesControlPts.push_back(templateControlPts);
	}

}

// get i-th reachset/template
Polyhedron* DiscreteDynamicalSystem::getReachSet(int i){
	if( i < 0 || i >= this->reachSets.size() ){
		cout<<"DiscreteDynamicalSystem::getReachSet : i must be between 0 and"<<this->reachSets.size()-1;
		exit (EXIT_FAILURE);
	}
	return this->reachSets[i];
}


lst DiscreteDynamicalSystem::getVars(){ return this->vars; }
lst DiscreteDynamicalSystem::getParams(){ return this->params; }
lst DiscreteDynamicalSystem::getDynamics(){ return this->dynamics; }

vector< Polyhedron* > DiscreteDynamicalSystem::getReachSets(){ return this->reachSets; }
vector< lst > DiscreteDynamicalSystem::getFogs(){ return this->fogs; }
vector< vector< lst > > DiscreteDynamicalSystem::getTemplatesControlPts(){ return this->templatesControlPts; }

DiscreteDynamicalSystem::~DiscreteDynamicalSystem() {
	// TODO Auto-generated destructor stub
}

