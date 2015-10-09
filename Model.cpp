/*
 * Model.cpp
 *
 *  Created on: Jan 7, 2015
 *      Author: Tommaso Dreossi
 */

#include "Model.h"

Model::Model(){

}

// normalizes the vectors with origin q and end points in vectors
pair< vector< vector<double> >, vector< double > > Model::normalizeVectors(vector<double> q, vector< vector<double> > vectors){


	for(int i=0; i<(signed)vectors.size(); i++){
		if( q.size() != vectors[i].size() ){
			cout<<"Model::normalizeVectors : q and vectors must have same dimensions";
			exit (EXIT_FAILURE);
		}
	}

	// shift vectors the origin and compute their norms
	vector< vector< double > > shifted_vectors;
	vector< double > norms;
	for(int i=0; i<(signed)vectors.size(); i++){
		vector< double > vi;
		double accum = 0;
		for(int j=0; j<(signed)vectors[i].size(); j++){
			vi.push_back(vectors[i][j]-q[j]);
			accum = accum + vi[j]*vi[j];
		}
		shifted_vectors.push_back(vi);
		norms.push_back(sqrt(accum));
	}

	// normalize vectors
	vector< vector< double > > norm_vectors;
	for(int i=0; i<(signed)shifted_vectors.size(); i++){
		vector< double > norm_vec;
		for(int j=0; j<(signed)shifted_vectors[i].size(); j++){
			norm_vec.push_back(shifted_vectors[i][j]/norms[i]);
		}
		norm_vectors.push_back(norm_vec);
	}

	pair< vector< vector<double> >, vector< double > > result (norm_vectors,norms);
	return result;

}

// generate a polynomial of degree deg
ex Model::genPoly( ex var, int deg ){
	if(deg == 0){
		return 1;
	}else{
		return power(var,deg) + this->genPoly(var,deg-1);
	}
}


ex Model::genPoly( lst vars, int deg ){

	if(vars.nops() == 1){
		return this->genPoly(vars[0],deg);
	}else{
		lst sub_vars;
		for(int i=0; i<vars.nops()-1; i++){
			sub_vars.append(vars[i]);
		}
		ex sub_poly = this->genPoly(sub_vars,deg);
		ex poly;
		poly = 0;
		for(int i=0; i<=deg; i++){
			poly = poly + power(vars[vars.nops()-1],i)*sub_poly;
		}
		return poly;
	}

}

Model::~Model() {
	// TODO Auto-generated destructor stub
}

