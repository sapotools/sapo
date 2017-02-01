/**
 * @file SIR.cpp
 * Van Der Pol model
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#include "VanDerPol.h"

 VanDerPol::VanDerPol(){


   ///// The dynamical system /////

 	// System dimension (number of variables)
 	int dim_sys = 2;
 	// List of state variables
 	symbol x("x"), y("y");
  this->vars = {x, y};

 	// System's dynamics
 	ex dx = x + (y)*0.02;
 	ex dy = y + (0.5*(1-x*x)*y - x)*0.02;
 	lst dyns;
 	this->dyns = {dx,dy};


 	///// Parallelotope bundle for reachable set representation /////

 	int num_dirs = 4;		// number of bundle directions
 	int num_temps = 6;		// number of bundle templates

 	// Directions matrix
 	vector< double > Li (dim_sys,0);
 	vector< vector< double > > L (num_dirs,Li);
 	L[0][0] = 1;
 	L[1][1] = 1;
 	L[2][0] = -1; L[2][1] = 1;
 	L[3][0] = 1; L[3][1] = 1;

 	// Template matrix
 	vector< int > Ti (dim_sys,0);
 	vector< vector< int > > T (num_temps,Ti);
 	T[0][0] = 0; T[0][1] = 1;
 	T[1][0] = 0; T[1][1] = 2;
 	T[2][0] = 0; T[2][1] = 3;
 	T[3][0] = 1; T[3][1] = 2;
 	T[4][0] = 1; T[4][1] = 3;
 	T[5][0] = 2; T[5][1] = 3;

 	// Offsets for the set of initial conditions
 	vector< double > offp (num_dirs,0);
 	vector< double > offm (num_dirs,0);
 	offp[0] = 0.01; offm[0] = 0;
 	offp[1] = 2; offm[1] = -1.99;
 	offp[2] = 10; offm[2] = 10;
 	offp[3] = 10; offm[3] = 10;

 	Bundle *B = new Bundle(L,offp,offm,T);
  this->reachSet = B;

 }
