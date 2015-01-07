/*
 * Ebola.cpp
 *
 *  Created on: Jan 7, 2015
 *      Author: dreossi
 */

#include "Ebola.h"

Ebola::Ebola() {

	// List of state variables and parameters
	symbol s("s"), e("e"), q("q"), i("i"), r("r");
	symbol beta("beta"), kappa1("kappa1"), kappa2("kappa2"), gamma1("gamma1"), gamma2("gamma2"), sigma("sigma");
	lst vars, params, dyns;
	vars = s, e, q, i, r;
	params = kappa1, gamma1;

	// System's dynamics
	ex ds = gamma1*q - (s*beta*i)/1000;
	ex de = (s*beta*i)/1000 - (kappa1+kappa2)*e;
	ex dq = kappa1*e - (gamma1+gamma2)*q;
	ex di = gamma2*q + kappa2*e - sigma*i;
	ex dr = sigma*i;
	dyns = ds,de,dq,di,dr;

	// Initialize uncontrollable parameters
	lst uncontr_params;
	uncontr_params.append(beta==0.9);
	uncontr_params.append(kappa2 == 0.5);
	uncontr_params.append(gamma2 == 0.5);
	uncontr_params.append(sigma == 0.28);
	dyns[0] = dyns[0].subs( uncontr_params );
	dyns[1] = dyns[1].subs( uncontr_params );
	dyns[2] = dyns[2].subs( uncontr_params );
	dyns[3] = dyns[3].subs( uncontr_params );
	dyns[4] = dyns[4].subs( uncontr_params );

	this->vars = vars;
	this->params = params;
	this->dyns = dyns;

	// Initialize the template for the rechability set
	symbol q1("q1"), q2("q2"), q3("q3"), q4("q4"), q5("q5");	// base vertex variables
	symbol a1("a1"), a2("a2"), a3("a3"), a4("a4"), a5("a5");	// parallelotope variables
	symbol b1("b1"), b2("b2"), b3("b3"), b4("b4"), b5("b5");	// amplitude variables
	lst qs, as, bs;
	qs = q1,q2,q3,q4,q5;
	as = a1,a2,a3,a4,a5;
	bs = b1,b2,b3,b4,b5;

	vector<lst> set_vars;
	set_vars.push_back(qs);
	set_vars.push_back(as);
	set_vars.push_back(bs);

	// Initialize parallelotope versors
	vector<double> ui (5,0);
	vector< vector<double> > u (5,ui);
	u[0][0] = 1;
	u[1][1] = 1;
	u[2][2] = 1;
	u[3][3] = 1;
	u[4][4] = 1;

	// Create the parallelotope
	Polyhedron *reach_set = new Parallelotope(set_vars,u);

	this->reach_set = reach_set;

	// Declare the initial parameter set as a linear system
	vector<double> pAi (2,0);
	vector< vector<double> > pA (4,pAi);
	vector<double> pb (4,0);
	pA[0][0] = 1; pA[0][1] = 0;		pb[0] = 0.3;
	pA[1][0] = -1; pA[1][1] = 0;	pb[1] = -0.2;
	pA[2][0] = 0; pA[2][1] = 1;		pb[2] = 0.5;
	pA[3][0] = 0; pA[3][1] = -1;	pb[3] = -0.2;

	LinearSystem *parameters = new LinearSystem(pA,pb);
	LinearSystemSet *parameter_set = new LinearSystemSet(parameters);

	this->init_para_set = parameter_set;

	// Set the specification
	ex constraint1 = q-50; 		//atoms
	ex constraint2 = -e+100;
	ex constraint3 = -q+25;
	STL *phi1 = new Atom(constraint1);
	STL *phi2 = new Atom(constraint2);
	STL *phi3 = new Atom(constraint3);
	STL *phi23 = new Disjunction(phi2,phi3);	// phi2 or phi3
	STL *phi = new Until(phi1,5,15,phi23);		// phi1 U (phi2 or phi3)

	this->spec = phi;

	this->name = "Ebola";
}

Ebola::~Ebola() {
	// TODO Auto-generated destructor stub
}

