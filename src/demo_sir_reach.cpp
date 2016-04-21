/**
 * @file demo_sir_reach.cpp
 * Reachability analysis on SIR epidemic model
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#include <stdio.h>
#include <iostream>

#include "Common.h"
#include "Bundle.h"
#include "Sapo.h"

using namespace std;

int main(int argc,char** argv){

	//The model
	int dim_sys = 3;

	// List of state variables and parameters
	symbol s("s"), i("i"), r("r"), beta("beta"), gamma("gamma");
	lst vars, dyns, params;
	vars = s, i, r;
	params = beta, gamma;

	// System's dynamics
	ex ds = s - (0.34*s*i)*0.1;				// susceptible
	ex di = i + (0.34*s*i - 0.05*i)*0.1;	// infected
	ex dr = r + 0.05*i*0.1;					// removed
	dyns = ds,di,dr;
	Model *sir = new Model(vars,params,dyns);


	// The initial set
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
//	L[5][0] = 1; L[5][1] = 0.5; L[5][2] = 0.5;
//	L[6][0] = 0; L[6][1] = 0.75; L[6][2] = 0.75;
//	L[7][0] = 1; L[7][1] = 0.2; L[7][2] = 0;

	// Offsets
	vector< double > offp (num_dirs,0);
	vector< double > offm (num_dirs,0);
	offp[0] = 0.8; offm[0] = -0.79;
	offp[1] = 0.2; offm[1] = -0.19;
	offp[2] = 0.0001; offm[2] = -0.000099;
	offp[3] = 1; offm[3] = 0;
	offp[4] = 1; offm[4] = 0;
//	offp[5] = 1; offm[5] = 0;
//	offp[6] = 1; offm[6] = 0;
//	offp[7] = 1; offm[7] = 0;

	// Template matrix
	vector< int > Ti (dim_sys,0);
	vector< vector< int > > T (num_temps,Ti);

	T[0][0] = 0; T[0][1] = 1; T[0][2] = 2;
	T[1][0] = 1; T[1][1] = 2; T[1][2] = 3;
	T[2][0] = 2; T[2][1] = 3; T[2][2] = 4;
//	T[3][0] = 5; T[3][1] = 2; T[3][2] = 1;

//	// Declare the initial parameter set as a linear system
//	vector<double> pAi (2,0);
//	vector< vector<double> > pA (4,pAi);
//	vector<double> pb (4,0);
//	pA[0][0] = 1; pA[0][1] = 0; pb[0] = 0.36;
//	pA[1][0] = -1; pA[1][1] = 0; pb[1] = -0.35;
//	pA[2][0] = 0; pA[2][1] = 1; pb[2] = 0.06;
//	pA[3][0] = 0; pA[3][1] = -1; pb[3] = -0.05;

	Bundle *B = new Bundle(L,offp,offm,T);


	// Options for Sapo
	sapo_opt options;
	options.trans = 1;
	options.alpha = 0.5;
	options.decomp = 0;
	options.verbose = false;

	Sapo *sapo = new Sapo(sir,options);
	Flowpipe* flowpipe = sapo->reach(B,300);

	// Store the constructed flowpipe in file sir.m
	char file_name[] = "sir.m";
	flowpipe->plotRegionToFile(file_name,'b');

}


