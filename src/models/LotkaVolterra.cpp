/**
 * @file LotkaVolterra.cpp
 * Lotka-Volterra model
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#include "LotkaVolterra.h"

/**
 * Constructor that instantiates a parallelotope
 *
 * @param[in] vars vector with the list of variables used for the generator functions
 * @param[in] u collection of generator versors
 */

 LotkaVolterra::LotkaVolterra(){


   ///// The dynamical system /////

 	// System dimension (number of variables)
 	int dim_sys = 5;
 	// List of state variables
 	symbol x1("x1"), x2("x2"), x3("x3"), x4("x4"), x5("x5");
 	lst vars;
 	vars = {x1, x2, x3, x4, x5};

 	ex alpha = 0.85; ex beta = 0.5;
 	ex delta = 0.01;

 	// System's dynamics
 	ex dx1 = x1 + ( x1*(1 - (x1 + alpha*x2 + beta*x5)) )*delta;
 	ex dx2 = x2 + (x2*(1 - (x2 + alpha*x3 + beta*x1)) )*delta;
 	ex dx3 = x3 + (x3*(1 - (x3 + alpha*x4 + beta*x2)) )*delta;
 	ex dx4 = x4 + (x4*(1 - (x4 + alpha*x5 + beta*x3)) )*delta;
 	ex dx5 = x5 + (x5*(1 - (x5 + alpha*x1 + beta*x4)) )*delta;
 	lst dyns;
 	dyns = {dx1,dx2,dx3,dx4,dx5};

 	this->vars = vars;
  this->dyns = dyns;


 	///// Parallelotope bundle for reachable set representation /////

 	int num_dirs = 7;		// number of bundle directions
 	int num_temps = 3;		// number of bundle templates

 	// Directions matrix
 	vector< double > Li (dim_sys,0);
 	vector< vector< double > > L (num_dirs,Li);
 	for( int i=0; i<dim_sys; i++ ){
 		L[i][i] = 1;
 	}
 	L[5][0] = 1; L[5][1] = 1; L[5][2] = 1;
 	L[6][3] = -1; L[6][4] = 1; L[6][0] = -1;

 	// Template matrix
 	vector< int > Ti (dim_sys,0);
 	vector< vector< int > > T (num_temps,Ti);
 	T[0][0] = 0; T[0][1] = 1; T[0][2] = 2; T[0][3] = 3; T[0][4] = 4;
 	T[1][0] = 5; T[1][1] = 6; T[1][2] = 1; T[1][3] = 2; T[1][4] = 3;
 	T[2][0] = 5; T[2][1] = 6; T[2][2] = 2; T[2][3] = 3; T[2][4] = 4;

 	// Offsets for the set of initial conditions
 	vector< double > offp (num_dirs,0);
 	vector< double > offm (num_dirs,0);
 	offp[0] = 1.0; offm[0] = -0.95;
 	offp[1] = 1.0; offm[1] = -0.95;
 	offp[2] = 1.0; offm[2] = -0.95;
 	offp[3] = 1.0; offm[3] = -0.95;
 	offp[4] = 1.0; offm[4] = -0.95;

 	offp[5] = 10.0; offm[5] = 1;
 	offp[6] = 10.0; offm[6] = 1;


 	Bundle *B = new Bundle(L,offp,offm,T);
  this->reachSet = B;

 }
