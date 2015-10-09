/*
 * SIR.cpp
 *
 *  Created on: Jan 7, 2015
 *      Author: Tommaso Dreossi
 */

#include "SIR.h"

SIR::SIR() {

	// List of state variables and parameters
	symbol s("s"), i("i"), r("r");
	symbol beta("beta"), gamma("gamma");
	lst vars, params, dyns;
	vars = s, i,r;
	params = beta,gamma;

	// System's dynamics
	ex ds = -beta*s*i;
	ex di = beta*s*i - gamma*i;
	ex dr = gamma*i;
	dyns = ds,di,dr;

	this->vars = vars;
	this->params = params;
	this->dyns = dyns;

	// Initialize the template for the rechability sets
	symbol q1("q1"), q2("q2"), q3("q3");	// base vertex variables
	symbol a1("a1"), a2("a2"), a3("a3");	// parallelotope variables
	symbol b1("b1"), b2("b2"), b3("b3");	// amplitude variables
	lst qs, as, bs;
	qs = q1,q2,q3;
	as = a1,a2,a3;
	bs = b1,b2,b3;

	vector<lst> set_vars;
	set_vars.push_back(qs);
	set_vars.push_back(as);
	set_vars.push_back(bs);

	//ciao


	// Initialize Template 1 versors
	vector<double> ui1 (3,0);
	vector< vector<double> > u1 (3,ui1);
	//u1[0][0] = 0.7071; u1[0][1] = -0.7071;
	u1[0][0] = 1;
	u1[1][1] = 1;
	u1[2][2] = 1;
	// Create Template1 parallelotope
	Polyhedron *reach_set_1 = new Parallelotope(set_vars,u1);
	this->reach_sets.push_back(reach_set_1);


//	vector<double> origin (3,0);
//	origin[0] = 0.8;
//	origin[1] = 0.15;
//	origin[2] = 0;
//
//	vector<double> vi2 (3,0);
//	vector< vector<double> > v2 (3,vi2);
//	v2[0][0] = 0.801; v2[0][1] = 0.201; v2[0][2] = 0.001;
//	v2[1][0] = 0.8-0.001; v2[1][1] = 0.201; v2[1][2] = 0.001;
//	v2[2][0] = 0.8; v2[2][1] = 0.201; v2[2][2] = 0;
//
//	pair< vector< vector<double> >, vector<double> > norm2 = this->normalizeVectors(origin,v2);
//
//	// Initialize Template 2 versors
//	vector< vector<double> > u2 = norm2.first;
//
//	// Create Template1 parallelotope
//	Polyhedron *reach_set_2 = new Parallelotope(set_vars,u2);
//	this->reach_sets.push_back(reach_set_2);

//	// Initialize Template 3 versors
//	vector<double> ui3 (3,0);
//	vector< vector<double> > u3 (3,ui3);
//	u3[0][0] = 1;	u3[0][1] = 0; u3[0][2] = 0;
//	u3[1][0] = 0.7071;	u3[1][1] = 0.7071; u3[1][2] = 0;
//	u3[2][0] = 0;	u3[2][1] = 0; u3[2][2] = 1;
//	// Create Template1 parallelotope
//	Polyhedron *reach_set_3 = new Parallelotope(set_vars,u3);
//	this->reach_sets.push_back(reach_set_3);


	// Declare initial conditions Template 1
	// base vertex
	vector< double > base_vertex_1 (3,0);
	base_vertex_1[0] = 0.7;
	base_vertex_1[1] = 0.3;
	base_vertex_1[2] = 0.0;
	// and versors lenghts
	vector< double > lenghts_1 (3,0);
	lenghts_1[0] = 0.001;
	lenghts_1[1] = 0.001;
	lenghts_1[2] = 0.001;

	poly_values ini_cond_1;
	ini_cond_1.base_vertex = base_vertex_1;
	ini_cond_1.lenghts = lenghts_1;
	this->initial_conditions.push_back(ini_cond_1);

//	// Declare initial conditions Template 2
//	// base vertex
//	vector< double > base_vertex_2 = origin;
//	// and versors lenghts
//	vector< double > lenghts_2 = norm2.second;
//
//	poly_values ini_cond_2;
//	ini_cond_2.base_vertex = base_vertex_2;
//	ini_cond_2.lenghts = lenghts_2;
//	this->initial_conditions.push_back(ini_cond_2);

//	// Declare initial conditions Template 3
//	// base vertex
//	vector< double > base_vertex_3 (3,0);
//	base_vertex_3[0] = 8-0.001;
//	base_vertex_3[1] = 2;
//	base_vertex_3[0] = 0;
//	// and versors lenghts
//	vector< double > lenghts_3 (3,0);
//	lenghts_3[0] = 0.002;
//	lenghts_3[1] = 0.0014;
//	lenghts_3[2] = 0.001;
//
//	poly_values ini_cond_3;
//	ini_cond_3.base_vertex = base_vertex_3;
//	ini_cond_3.lenghts = lenghts_3;
//	this->initial_conditions.push_back(ini_cond_3);


	// Declare the initial parameter set as a linear system
	vector<double> pAi (2,0);
	vector< vector<double> > pA (4,pAi);
	vector<double> pb (4,0);
	pA[0][0] = 1; pA[0][1] = 0; pb[0] = 0.36;
	pA[1][0] = -1; pA[1][1] = 0; pb[1] = -0.35;
	pA[2][0] = 0; pA[2][1] = 1; pb[2] = 0.06;
	pA[3][0] = 0; pA[3][1] = -1; pb[3] = -0.05;

	LinearSystem *parameters = new LinearSystem(pA,pb);
	LinearSystemSet *parameter_set = new LinearSystemSet(parameters);

	this->init_para_set = parameter_set;

	// Set the specification
	ex constraint1 = i-0.673; 		//atoms
	STL *phi1 = new Atom(constraint1);
	STL* phi = new Always(0,50,phi1);

	this->spec = phi;

	this->name = "SI";
}

SIR::~SIR() {
	// TODO Auto-generated destructor stub
}

