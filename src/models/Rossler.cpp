/**
 * @file Rossler.cpp
 * Rossler model
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#include "Rossler.h"

/**
 * Constructor that instantiates a parallelotope
 *
 * @param[in] vars vector with the list of variables used for the generator functions
 * @param[in] u collection of generator versors
 */

 Rossler::Rossler(){


   ///// The dynamical system /////

 	// System dimension (number of variables)
 	int dim_sys = 3;
 	// List of state variables
 	symbol x("x"), y("y"), z("z");
 	lst vars;
 	vars = {x, y, z};

 	// System's dynamics
 	ex dx = x + (-y-z)*0.025;
 	ex dy = y + (x + 0.1*y)*0.025;
 	ex dz = z + (0.1 + z*(x-14))*0.025;
 	lst dyns;
 	dyns = {dx,dy,dz};

  this->vars = vars;
  this->dyns = dyns;

 	///// Parallelotope bundle for reachable set representation /////

 	int num_dirs = 5;		// number of bundle directions
 	int num_temps = 3;		// number of bundle templates

 	// Directions matrix
 	vector< double > Li (dim_sys,0);
 	vector< vector< double > > L (num_dirs,Li);
 	L[0][0] = 1;
 	L[1][1] = 1;
 	L[2][2] = 1;
 	L[3][0] = 1; L[3][1] = 0.5;
 	L[4][0] = 0.5; L[4][2] = 0.5;

 	// Template matrix
 	vector< int > Ti (dim_sys,0);
 	vector< vector< int > > T (num_temps,Ti);
 	T[0][0] = 0; T[0][1] = 1; T[0][2] = 2;
 	T[1][0] = 1; T[1][1] = 2; T[1][2] = 3;
 	T[2][0] = 2; T[2][1] = 3; T[2][2] = 4;

 	// Offsets for the set of initial conditions
 	vector< double > offp (num_dirs,0);
 	vector< double > offm (num_dirs,0);
 	offp[0] = 0.1; offm[0] = -0.09;
 	offp[1] = 5; offm[1] = -4.99;
 	offp[2] = 0.1; offm[2] = -0.09;
 	offp[3] = 10; offm[3] = 0;
 	offp[4] = 10; offm[4] = 0;

 	Bundle *B = new Bundle(L,offp,offm,T);
  this->reachSet = B;

 }
