/**
 * @file Phosphorelay.cpp
 * Phosphorelay model
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#include "Phosphorelay.h"

 Phosphorelay::Phosphorelay(){

   ///// The dynamical system /////

 	// System dimension (number of variables)
  strcpy(this->name,"Phospohorelay");
 	int dim_sys = 7;
 	// List of state variables
 	symbol x1("x1"), x2("x2"), x3("x3"), x4("x4"), x5("x5"), x6("x6"), x7("x7");
 	lst vars;
 	vars = {x1, x2, x3, x4, x5, x6, x7};

 	ex delta = 0.01;

 	// System's dynamics
 	ex dx1 = x1 + ( -0.4*x1 + 5*x3*x4 )*delta;
 	ex dx2 = x2 + ( 0.4*x1 - x2 )*delta;
 	ex dx3 = x3 + ( x2-5*x3*x4 )*delta;
 	ex dx4 = x4 + ( 5*x5*x6 - 5*x3*x4 )*delta;
 	ex dx5 = x5 + ( -5*x5*x6 + 5*x3*x4 )*delta;
 	ex dx6 = x6 + ( 0.5*x7 - 5*x5*x6 )*delta;
 	ex dx7 = x7 + ( -0.5*x7 + 5*x5*x6 )*delta;
 	lst dyns;
 	dyns = {dx1,dx2,dx3,dx4,dx5,dx6,dx7};

 	this->vars = vars;
  this->dyns = dyns;


 	///// Parallelotope bundle for reachable set representation /////

 	int num_dirs = 10;		// number of bundle directions
 	int num_temps = 4;		// number of bundle templates

 	// Directions matrix
 	vector< double > Li (dim_sys,0);
 	vector< vector< double > > L (num_dirs,Li);
 	for( int i=0; i<dim_sys; i++ ){
 		L[i][i] = 1;
 	}
 	L[7][2] = 1; L[7][3] = 1;
 	L[8][4] = 1; L[8][5] = 1;
 	L[9][2] = 1; L[9][3] = 1; L[9][4] = 1; L[9][5] = 1;

 	// Template matrix
 	vector< int > Ti (dim_sys,0);
 	vector< vector< int > > T (num_temps,Ti);
 	T[0][0] = 0; T[0][1] = 1; T[0][2] = 2; T[0][3] = 3; T[0][4] = 4; T[0][5] = 5; T[0][6] = 6;
 	T[1][0] = 0; T[1][1] = 1; T[1][2] = 2; T[1][3] = 7; T[1][4] = 4; T[1][5] = 5; T[1][6] = 6;
 	T[2][0] = 0; T[2][1] = 1; T[2][2] = 2; T[2][3] = 7; T[2][4] = 8; T[2][5] = 5; T[2][6] = 6;
 	T[3][0] = 0; T[3][1] = 1; T[3][2] = 2; T[3][3] = 7; T[3][4] = 9; T[3][5] = 5; T[3][6] = 6;

 	// Offsets for the set of initial conditions
 	vector< double > offp (num_dirs,0);
 	vector< double > offm (num_dirs,0);
 	offp[0] = 1.01; offm[0] = -1.00;
 	offp[1] = 1.01; offm[1] = -1.00;
 	offp[2] = 1.01; offm[2] = -1.00;
 	offp[3] = 1.01; offm[3] = -1.00;
 	offp[4] = 1.01; offm[4] = -1.00;
 	offp[5] = 1.01; offm[5] = -1.00;
 	offp[6] = 1.01; offm[6] = -1.00;

 	offp[7] = 100; offm[7] = 100;
 	offp[8] = 100; offm[8] = 100;
 	offp[9] = 100; offm[9] = 100;
 	//offp[10] = 100; offm[10] = 100;


 	Bundle *B = new Bundle(L,offp,offm,T);
  this->reachSet = B;

 }
