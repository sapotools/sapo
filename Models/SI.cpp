/*
 * SI.cpp
 *
 *  Created on: Jan 7, 2015
 *      Author: Tommaso Dreossi
 */

#include "SI.h"

SI::SI() {

	// List of state variables and parameters
	symbol s("s"), i("i");
	symbol beta("beta");
	lst vars, params, dyns;
	vars = s, i;
	params = beta;

	// System's dynamics
	ex ds = -beta*s*i;
	ex di = beta*s*i;
	dyns = ds,di;

	this->vars = vars;
	this->params = params;
	this->dyns = dyns;

	// Initialize the template for the rechability sets
	symbol q1("q1"), q2("q2");	// base vertex variables
	symbol a1("a1"), a2("a2");	// parallelotope variables
	symbol b1("b1"), b2("b2");	// amplitude variables
	lst qs, as, bs;
	qs = q1,q2;
	as = a1,a2;
	bs = b1,b2;

	vector<lst> set_vars;
	set_vars.push_back(qs);
	set_vars.push_back(as);
	set_vars.push_back(bs);

	// Initialize Template 1 versors
	vector<double> ui1 (2,0);
	vector< vector<double> > u1 (2,ui1);
	u1[0][0] = 1;
	u1[1][1] = 1;
	// Create Template1 parallelotope
	Polyhedron *reach_set_1 = new Parallelotope(set_vars,u1);
	this->reach_sets.push_back(reach_set_1);

	// Initialize Template 2 versors
	vector<double> ui2 (2,0);
	vector< vector<double> > u2 (2,ui2);
	u2[0][0] = 1;	u2[0][1] = 0;
	u2[1][0] = 0.4472;	u2[1][1] = 0.8999;
	// Create Template1 parallelotope
	Polyhedron *reach_set_2 = new Parallelotope(set_vars,u2);
	this->reach_sets.push_back(reach_set_2);

	// Initialize Template 3 versors
	vector<double> ui3 (2,0);
	vector< vector<double> > u3 (2,ui3);
	u3[0][0] = 1;	u3[0][1] = -0.5;
	u3[1][0] = 0;	u3[1][1] = 1;
	// Create Template1 parallelotope
	Polyhedron *reach_set_3 = new Parallelotope(set_vars,u3);
	this->reach_sets.push_back(reach_set_3);


	// Declare initial conditions Template 1
	// base vertex
	vector< double > base_vertex_1 (2,0);
	base_vertex_1[0] = 8;
	base_vertex_1[1] = 1;
	// and versors lenghts
	vector< double > lenghts_1 (2,0);
	lenghts_1[0] = 1;
	lenghts_1[1] = 1;

	poly_values ini_cond_1;
	ini_cond_1.base_vertex = base_vertex_1;
	ini_cond_1.lenghts = lenghts_1;
	this->initial_conditions.push_back(ini_cond_1);

	// Declare initial conditions Template 2
	// base vertex
	vector< double > base_vertex_2 (2,0);
	base_vertex_2[0] = 7.5;
	base_vertex_2[1] = 1;
	// and versors lenghts
	vector< double > lenghts_2 (2,0);
	lenghts_2[0] = 1.5;
	lenghts_2[1] = 1.1180;

	poly_values ini_cond_2;
	ini_cond_2.base_vertex = base_vertex_2;
	ini_cond_2.lenghts = lenghts_2;
	this->initial_conditions.push_back(ini_cond_2);

	// Declare initial conditions Template 3
	// base vertex
	vector< double > base_vertex_3 (2,0);
	base_vertex_3[0] = 8;
	base_vertex_3[1] = 1;
	// and versors lenghts
	vector< double > lenghts_3 (2,0);
	lenghts_3[0] = 1.1180;
	lenghts_3[1] = 1.5;

	poly_values ini_cond_3;
	ini_cond_3.base_vertex = base_vertex_3;
	ini_cond_3.lenghts = lenghts_3;
	this->initial_conditions.push_back(ini_cond_3);


	// Declare the initial parameter set as a linear system
	vector<double> pAi (1,0);
	vector< vector<double> > pA (2,pAi);
	vector<double> pb (4,0);
	pA[0][0] = 1; pb[0] = 0.2;
	pA[1][0] = -1; pb[1] = -0.19;

	LinearSystem *parameters = new LinearSystem(pA,pb);
	LinearSystemSet *parameter_set = new LinearSystemSet(parameters);

	this->init_para_set = parameter_set;

	// Set the specification
	ex constraint1 = s-56; 		//atoms
	ex constraint2 = s-1000;
	ex constraint3 = s-2500;
	STL *phi = new Atom(constraint1);
	STL *phi2 = new Atom(constraint2);
	STL *phi3 = new Atom(constraint3);
	//STL *phi = new Disjunction(phi2,phi3);	// phi2 or phi3
	//STL *phi = new Until(phi1,10,15,phi2);		// phi1 U (phi2 or phi3)
	//STL *phi = new Eventually(10,15,phi1);		// phi1 U (phi2 or phi3)

	this->spec = phi;

	this->name = "SI";
}

SI::~SI() {
	// TODO Auto-generated destructor stub
}

