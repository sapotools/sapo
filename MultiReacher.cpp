/*
 * MultiReacher.cpp
 *
 *  Created on: Jun 3, 2015
 *      Author: dreossi
 */

#include "MultiReacher.h"

MultiReacher::MultiReacher(lst vars, lst dynamics, MultiParallelotope *D, lst ls) {

	this->vars = vars;
	this->dynamics = dynamics;

	this->qs = D->getBaseVertexVars();
	this->as = D->getFreeVars();
	this->bs = D->getLengthsVars();
	this->us = D->getGenVars();
	this->ls = ls;

	lst generator_functions = D->getGenFuns();
	lst sub, fog;

	// compute the composition fog
	for(int i=0; i<(signed)this->vars.nops(); i++){
		sub.append(this->vars[i] == generator_functions[i]);
	}
	for(int i=0; i<(signed)this->dynamics.nops(); i++){
		fog.append(this->dynamics[i].subs(sub));
	}
	// compose with the template directions
	ex lambda_fog = 0;
	for(int i=0; i<(signed)this->ls.nops(); i++){
		lambda_fog = lambda_fog + this->ls[i]*fog[i];
	}

	cout<<"Computing symbolic control points...\n";
//	// compute the control points symbolically
//	BaseConverter *converter = new BaseConverter(this->as,lambda_fog);
//	this->dyn_bern_coeffs = converter->getBernCoeffs();

}

// find the maximum control point of the generator function composed with dynamics and template
// with numerical values specified by the user
double MultiReacher::maxBernCoeffStatic(int parall_idx, int dir_idx, vector<double> base_vertex, vector<double> lenghts){

	// prepare the substitutions
	exmap map;
	for(int i=0; i<(signed)base_vertex.size(); i++){
		map[this->qs[i]] = base_vertex[i];
		map[this->bs[i]] = lenghts[i];
	}

	double max_bern_coeff = -DBL_MAX;
	lst bern_coeffs = this->static_bern_coeffs[parall_idx][dir_idx];

	for(int i=0; i<(signed)bern_coeffs.nops(); i++){
		// substitute symb variables and convert them to double
		double num_bern_coeff = ex_to<numeric>(evalf(bern_coeffs[i].subs(map,subs_options::no_pattern))).to_double();
		// store the maximum coeff
		max_bern_coeff = max(max_bern_coeff,num_bern_coeff);
	}

	return max_bern_coeff;

}

// find the maximum control point of the generator function composed with dynamics and template
// with numerical values specified by the user
double MultiReacher::maxBernCoeffNumerical(int parall_idx, int dir_idx, vector<double> base_vertex, vector<double> lenghts){

	// prepare the substitutions
	exmap map;
	for(int i=0; i<(signed)base_vertex.size(); i++){
		map[this->qs[i]] = base_vertex[i];
		map[this->bs[i]] = lenghts[i];
	}

	double max_bern_coeff = -DBL_MAX;

	ex numeric_lambda = this->Lambda_fogs[parall_idx][dir_idx].subs(map,subs_options::no_pattern);
	BaseConverter *BC = new BaseConverter(this->as,numeric_lambda);
	BC->print();
	BC->implicitMaxIndex();
	lst bern_coeffs = BC->getBernCoeffs();

	for(int i=0; i<(signed)bern_coeffs.nops(); i++){
		// substitute symb variables and convert them to double
		double num_bern_coeff = ex_to<numeric>(bern_coeffs[i]).to_double();
		// store the maximum coeff
		max_bern_coeff = max(max_bern_coeff,num_bern_coeff);
	}

	return max_bern_coeff;

}

// find the maximum control point of the generator function composed with dynamics and template
// with numerical values specified by the user
double MultiReacher::maxBernCoeffDirs(vector<double> base_vertex, vector<double> lenghts, vector< vector<double> > generators, vector<double> temp_direction){

	// prepare the substitutions
	exmap map;
	for(int i=0; i<(signed)base_vertex.size(); i++){
		map[this->qs[i]] = base_vertex[i];
		map[this->bs[i]] = lenghts[i];
		map[this->ls[i]]= temp_direction[i];
		for(int j=0; j<(signed)generators[i].size(); j++){
			map[this->us[i][j]] = generators[j][i];
		}
	}

	double max_bern_coeff = -DBL_MAX;
	for(int i=0; i<(signed)this->dyn_bern_coeffs.nops(); i++){
		// substitute symb variables and convert them to double
		double num_bern_coeff = ex_to<numeric>(evalf(this->dyn_bern_coeffs[i].subs(map,subs_options::no_pattern))).to_double();
		// store the maximum coeff
		max_bern_coeff = max(max_bern_coeff,num_bern_coeff);
	}

	return max_bern_coeff;

}

// Perform a single with the fixed parallelotopes discribed by the template matrix
MultiParallelotope* MultiReacher::reachStepStaticTemplate(MultiParallelotope *D){

	vector< double > upper_bounds(D->getNumDirs(),DBL_MAX);
	vector< double > lower_bounds(D->getNumDirs(),DBL_MAX);

	//D = this->reachStepStaticTemplate(D);

	for(int i=0; i<(signed)D->getCard(); i++){	// for each parallelotope

		LinearSystem *LS = D->getParallelotope(i);	// get i-th parallelotope
		vector< vector< double > > A = LS->getA();	// extract directions of parallelotope

		vector<lst> set_vars;
		set_vars.push_back(this->qs); set_vars.push_back(this->as); set_vars.push_back(this->bs);
		Parallelotope *P = new Parallelotope(set_vars,LS);

		vector< double > base_vertex = P->getBaseVertex();
		vector< double > lenghts = P->getLenghts();

		// get the template idxs of D
		vector<int> template_idxs = D->getTemplate(i);

		for(int j=0; j<(signed)template_idxs.size(); j++){

			int dir_index = template_idxs[j];

			// determine current maximum upper coefficient
			double max_coeff = this->maxBernCoeffStatic(i,j,base_vertex, lenghts);

			upper_bounds[dir_index] = min(upper_bounds[dir_index],max_coeff);

			// determine current maximum upper coefficient
			max_coeff = this->maxBernCoeffStatic(i, template_idxs.size()+j, base_vertex, lenghts);
			lower_bounds[dir_index] = min(lower_bounds[dir_index],max_coeff);

		}
	}

//	for(int i=0; i<upper_bounds.size(); i++){
//		cout<<upper_bounds[i]<<" "<<lower_bounds[i]<<"\n";
//	}


	// Create new reached set
	MultiParallelotope *Dp = new MultiParallelotope(D->getDirections(), upper_bounds, lower_bounds, D->getTemplates());
	return Dp;
}

MultiParallelotope* MultiReacher::reachStepNumerical(MultiParallelotope *D){

	vector< double > upper_bounds(D->getNumDirs(),DBL_MAX);
	vector< double > lower_bounds(D->getNumDirs(),DBL_MAX);

	//D = this->reachStepStaticTemplate(D);

	for(int i=0; i<(signed)D->getCard(); i++){	// for each parallelotope

		LinearSystem *LS = D->getParallelotope(i);	// get i-th parallelotope
		vector< vector< double > > A = LS->getA();	// extract directions of parallelotope

		vector<lst> set_vars;
		set_vars.push_back(this->qs); set_vars.push_back(this->as); set_vars.push_back(this->bs);
		Parallelotope *P = new Parallelotope(set_vars,LS);

		vector< double > base_vertex = P->getBaseVertex();
		vector< double > lenghts = P->getLenghts();

		// get the template idxs of D
		vector<int> template_idxs = D->getTemplate(i);

		for(int j=0; j<(signed)template_idxs.size(); j++){

			int dir_index = template_idxs[j];

			// determine current maximum upper coefficient
			double max_coeff = this->maxBernCoeffNumerical(i,j,base_vertex, lenghts);

			upper_bounds[dir_index] = min(upper_bounds[dir_index],max_coeff);

			// determine current maximum upper coefficient
			max_coeff = this->maxBernCoeffNumerical(i, template_idxs.size()+j, base_vertex, lenghts);
			lower_bounds[dir_index] = min(lower_bounds[dir_index],max_coeff);

		}
	}

//	for(int i=0; i<upper_bounds.size(); i++){
//		cout<<upper_bounds[i]<<" "<<lower_bounds[i]<<"\n";
//	}


	// Create new reached set
	MultiParallelotope *Dp = new MultiParallelotope(D->getDirections(), upper_bounds, lower_bounds, D->getTemplates());
	return Dp;
}

// Perform a single reach step combining all the directions with all the parallelotopes
MultiParallelotope* MultiReacher::reachStepDynTemplate(MultiParallelotope *D){

	vector< double > upper_bounds(D->getNumDirs(),DBL_MAX);
	vector< double > lower_bounds(D->getNumDirs(),DBL_MAX);

	// compose all the directions with the dynamics
	for(int i=0; i<D->getCard(); i++){		// for each parallelotope

		LinearSystem *LS = D->getParallelotope(i);	// get i-th parallelotope

		for(int j=0; j<D->getNumDirs(); j++){					// for each direction

			vector<double> lambda = D->getDirection(j);		// get j-th direction
			vector<lst> set_vars;
			set_vars.push_back(this->qs); set_vars.push_back(this->as); set_vars.push_back(this->bs);
			Parallelotope *P = new Parallelotope(set_vars,LS);

			vector< double > base_vertex = P->getBaseVertex();
			vector< double > lenghts = P->getLenghts();
			vector< vector<double> > u = P->getVersors();

			// determine current maximum upper coefficient
			double max_coeff = this->maxBernCoeffDirs(base_vertex, lenghts, u ,lambda);
			upper_bounds[j] = min(upper_bounds[j],max_coeff);

			// determine current maximum upper coefficient
			max_coeff = this->maxBernCoeffDirs(base_vertex, lenghts, u ,this->negate(lambda));
			lower_bounds[j] = min(lower_bounds[j],max_coeff);
		}
	}

	for(int i=0; i<upper_bounds.size(); i++){
		cout<<upper_bounds[i]<<" "<<lower_bounds[i]<<"\n";
	}

	// Create new reached set
	MultiParallelotope *Dp = new MultiParallelotope(D->getDirections(), upper_bounds, lower_bounds, D->getTemplates());
	return Dp;
}

// compute the reachability set from D up to k steps
void MultiReacher::reach(MultiParallelotope *D, int k){

	D->printParallelotopesIntersection();

	for(int i=0; i<k; i++){
		D = this->reachStepStaticTemplate(D);
		D->printParallelotopesIntersection();
	}
}

void MultiReacher::staticReach(MultiParallelotope *D, int k){

	cout<<"Computing static control points...\n";

	for(int i=0; i<D->getCard(); i++){		// for each parallelotope

		cout<<"Template "<<i<<"...\n";

		lst generator_function;
		generator_function = D->getParallelotopeGenFun(i);

		// compose generator function with dynamics
		lst sub;
		for(int j=0; j<this->vars.nops(); j++){
			sub.append(vars[j] == generator_function[j]);
		}
		lst fog;
		for(int j=0; j<(signed)dynamics.nops(); j++){
			fog.append(dynamics[j].subs(sub));
		}

		// compose with directions of the template
		vector< vector< double > > Lambdas = D->getParallelotope(i)->getA();
		vector< ex > Lambda_fog;
		for(int j=0; j<(signed)Lambdas.size(); j++){
			ex Lambda_fog_j = 0;
			for( int k=0; k<(signed)Lambdas[j].size(); k++ ){
				Lambda_fog_j = Lambda_fog_j + Lambdas[j][k]*fog[k];
			}
			Lambda_fog.push_back(Lambda_fog_j);
		}

		vector< lst > bern_coeffs;
		for(int j=0; j<(signed)Lambda_fog.size(); j++){
			cout<<"\tControl points "<<j<<"...\n";
			BaseConverter *BC = new BaseConverter(this->as, Lambda_fog[j]);
			BC->print();
			lst coeffs;
			bern_coeffs.push_back(BC->getBernCoeffs());
		}

		this->static_bern_coeffs.push_back(bern_coeffs);
	}


	cout<<"Static Reachability...\n";

	for(int i=0; i<k; i++){
		cout<<"\tReach step "<<i<<"...\n";
		D = this->reachStepStaticTemplate(D);
		pair< vector< double >, vector< double > > shrinked_bounds = D->shrink(D->parallelotopesIntersection());
		D->setBounds(shrinked_bounds.first, shrinked_bounds.second);
		D->printParallelotopesIntersection();
	}

}

void MultiReacher::numericalReach(MultiParallelotope *D, int k){

	cout<<"Numerical Reachability...\n";

	vector< vector< ex > >  Lambda_fogs;

	for(int i=0; i<D->getCard(); i++){		// for each parallelotope

		cout<<"Template "<<i<<"...\n";

		lst generator_function;
		generator_function = D->getParallelotopeGenFun(i);

		// compose generator function with dynamics
		lst sub;
		for(int j=0; j<this->vars.nops(); j++){
			sub.append(vars[j] == generator_function[j]);
		}
		lst fog;
		for(int j=0; j<(signed)dynamics.nops(); j++){
			fog.append(dynamics[j].subs(sub));
		}

		// compose with directions of the template
		vector< vector< double > > Lambdas = D->getParallelotope(i)->getA();
		vector< ex > Lambda_fog;
		for(int j=0; j<(signed)Lambdas.size(); j++){
			ex Lambda_fog_j = 0;
			for( int k=0; k<(signed)Lambdas[j].size(); k++ ){
				Lambda_fog_j = Lambda_fog_j + Lambdas[j][k]*fog[k];
			}
			Lambda_fog.push_back(Lambda_fog_j);
		}
		this->Lambda_fogs.push_back(Lambda_fog);

	}


	for(int s=0; s<k; s++){

		cout<<"\tReach step "<<s<<"...\n";

		D = this->reachStepNumerical(D);
		pair< vector< double >, vector< double > > shrinked_bounds = D->shrink(D->parallelotopesIntersection());
		D->setBounds(shrinked_bounds.first, shrinked_bounds.second);
		D->printParallelotopesIntersection();
	}

}


// generate a random direction of dimension n
vector<double> MultiReacher::generateRndDir(int n){

	vector<double> dir;

	// generate random vector
	for(int i=0; i<n; i++){
		dir.push_back((rand()%50)-25);
	}

	// return its normalization
	return this->normalize(dir);
}

// generate m directions of dimension n
vector< vector<double> > MultiReacher::generateRndDirs(int m, int n){
	vector< vector<double> > directions;
	for(int i=0; i<m; i++){
		directions.push_back(this->generateRndDir(n));
	}
	return directions;
}

// compute the norm of v
double MultiReacher::norm(vector<double> v){
	double norm = 0;
	for(int i=0; i<(signed)v.size(); i++){
		norm = norm + v[i]*v[i];
	}
	return sqrt(norm);
}

// normalize v
vector<double> MultiReacher::normalize(vector<double> v){
	double norm = this->norm(v);
	for(int i=0; i<(signed)v.size(); i++){
		v[i] = v[i]/norm;
	}
	return v;
}

// change sign to all elements of a vector
vector< double > MultiReacher::negate(vector< double > v){

	vector< double > minus_v;
	for(int i=0; i<v.size(); i++){
		minus_v.push_back(-v[i]);
	}
	return minus_v;
}

MultiReacher::~MultiReacher() {
	// TODO Auto-generated destructor stub
}

