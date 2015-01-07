/*
 * LinearSystem.cpp
 *
 *  Created on: Oct 24, 2014
 *      Author: dreossi
 */

#include "LinearSystem.h"

LinearSystem::LinearSystem(vector< vector<double> > A, vector< double > b){

	this->A = A;
	this->b = b;

	this->n_vars = this->A[0].size();
}

LinearSystem::LinearSystem(lst vars, lst constraints) {

	this->vars = vars;
	this->constraints = constraints;

	this->n_vars = this->vars.nops();

	initLS();	// initialize Linear System
}

// Initialize the Linear System
// computing A and b
void LinearSystem::initLS(){

	for(int i=0; i<(signed)this->constraints.nops(); i++){

		vector<double> Ai;
		ex const_term = this->constraints[i];

		for(int j=0; j<(signed)this->vars.nops(); j++){

			// Extract the coefficient of the i-th variable (grade 1)
			double coeff = ex_to<numeric>(evalf(this->constraints[i].coeff(this->vars[j],1))).to_double();
			Ai.push_back(coeff);

			// Project to obtain the constant term
			const_term = const_term.coeff(this->vars[j],0);
		}

		double bi = ex_to<numeric>(evalf(const_term)).to_double();

		this->A.push_back(Ai);
		this->b.push_back(-bi);

	}

}

// Get A
vector< vector<double> > LinearSystem::getA(){
	return this->A;
}

// Get b
vector<double> LinearSystem::getb(){
	return this->b;
}

// Get the i-j element of A
double LinearSystem::getA(int i, int j){
	if(( 0<= i ) && (i < (signed)this->A.size())){
		if(( 0<= j ) && (j < (signed)this->A[j].size())){
			return this->A[i][j];
		}
	}
	cout<<"LinearSystem::getA : i and j must be within the LS->A size";
	exit (EXIT_FAILURE);
}

// Get the i-th element of b
double LinearSystem::getb(int i){
	if(( 0<= i ) && (i < (signed)this->b.size())){
			return this->b[i];
	}
	cout<<"LinearSystem::getb : i and j must be within the LS->b size";
	exit (EXIT_FAILURE);
}

// Determine whether this linear system is empty or not
bool LinearSystem::isEmpty(){

	vector< vector< double > > extA = this->A;
	vector< double > obj_fun (this->n_vars, 0);
	obj_fun.push_back(1);

	// Add an extra variable to the linear system
	for(int i=0; i<extA.size(); i++){
		extA[i].push_back(-1);
	}

	double z = this->solveLinearSystem(extA,this->b,obj_fun,GLP_MIN);

	return (z>=0);

}

// Solve a linear system
double LinearSystem::solveLinearSystem(vector< vector< double > > A, vector< double > b, vector< double > obj_fun, int min_max){

	int num_rows = A.size();
	int num_cols = obj_fun.size();
	int size_lp = num_rows*num_cols;

	int ia[size_lp+1], ja[size_lp+1];
	double ar[size_lp+1], z;

	glp_prob *lp;
	lp = glp_create_prob();
	glp_set_obj_dir(lp, min_max);

	// Turn off verbose mode
	glp_smcp lp_param;
	glp_init_smcp(&lp_param);
	lp_param.msg_lev = GLP_MSG_ERR;

	glp_add_rows(lp, num_rows);
	for(int i=0; i<num_rows; i++){
		glp_set_row_bnds(lp, i+1, GLP_UP, 0.0, b[i]);
	}

	glp_add_cols(lp, num_cols);
	for(int i=0; i<num_cols; i++){
		glp_set_col_bnds(lp, i+1, GLP_FR, 0.0, 0.0);
	}

	for(int i=0; i<num_cols; i++){
		glp_set_obj_coef(lp, i+1, obj_fun[i]);
	}

	int k=1;
	for(int i=0; i<num_rows; i++){
		for(int j=0; j<num_cols; j++){
			ia[k] = i+1, ja[k] = j+1, ar[k] = A[i][j]; /* a[i+1,j+1] = A[i][j] */
			k++;
		}
	}

	glp_load_matrix(lp, size_lp, ia, ja, ar);
	glp_simplex(lp, &lp_param);
	return glp_get_obj_val(lp);

}

// minimize the function obj_fun
double LinearSystem::minLinearSystem(lst vars, ex obj_fun){

	vector< double > obj_fun_coeffs;
	ex const_term = obj_fun;

	// Extract the coefficient of the i-th variable (grade 1)
	for(int i=0; i<(signed)vars.nops(); i++){

		double coeff = ex_to<numeric>(evalf(obj_fun.coeff(vars[i],1))).to_double();

		obj_fun_coeffs.push_back(coeff);
		const_term = const_term.coeff(vars[i],0);

	}

	double c = ex_to<numeric>(evalf(const_term)).to_double();
	double min = this->solveLinearSystem(this->A,this->b,obj_fun_coeffs,GLP_MIN);

	return (min+c);

}

// maximize the function obj_fun
double LinearSystem::maxLinearSystem(lst vars, ex obj_fun){

	vector< double > obj_fun_coeffs;
	ex const_term = obj_fun;

	// Extract the coefficient of the i-th variable (grade 1)
	for(int i=0; i<(signed)vars.nops(); i++){

		double coeff = ex_to<numeric>(evalf(obj_fun.coeff(vars[i],1))).to_double();

		obj_fun_coeffs.push_back(coeff);
		const_term = const_term.coeff(vars[i],0);

	}

	double c = ex_to<numeric>(evalf(const_term)).to_double();
	double max = this->solveLinearSystem(this->A,this->b,obj_fun_coeffs,GLP_MAX);

	return (max+c);
}


// Create a new liner system by merging this LS and the specified one
LinearSystem* LinearSystem::appendLinearSystem(LinearSystem *LS){

	vector< vector<double> > newA = LS->getA();
	vector<double> newb = LS->getb();

	for(int i=0; i<this->A.size(); i++){
		newA.push_back( this->A[i] );
		newb.push_back( this->b[i] );
	}

	return new LinearSystem(newA,newb);

}

void LinearSystem::print(){
	for(int i=0; i<this->A.size(); i++){
		for(int j=0; j<this->A[i].size(); j++){
			cout<<this->A[i][j]<<" ";
		}
		cout<<" <= "<<this->b[i]<<"\n";
	}
}

LinearSystem::~LinearSystem() {
	// TODO Auto-generated destructor stub
}

