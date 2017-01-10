/**
 * @file demo_sir_reach.cpp
 * Demo: Reachability analysis of SIR epidemic model
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

/**
 * Main function
 *
 */
int main(int argc,char** argv){


	///// The dynamical system /////

	// System dimension (number of variables)
	int dim_sys = 3;
	// List of state variables
	symbol s("s"), i("i"), r("r");
	lst vars;
	vars = s, i, r;

	// System's dynamics
	ex ds = s - (0.34*s*i)*0.1;				// susceptible
	ex di = i + (0.34*s*i - 0.05*i)*0.1;	// infected
	ex dr = r + 0.05*i*0.1;					// removed
	lst dyns;
	dyns = ds,di,dr;

	Model *sir = new Model(vars,dyns);


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
	offp[0] = 0.8; offm[0] = -0.79;
	offp[1] = 0.2; offm[1] = -0.19;
	offp[2] = 0.0001; offm[2] = -0.000099;
	offp[3] = 1; offm[3] = 0;
	offp[4] = 1; offm[4] = 0;

	Bundle *B = new Bundle(L,offp,offm,T);


	///// SAPO core /////

	// Sapo's options
	sapo_opt options;
	options.trans = 1;			// Set transformation (0=OFO, 1=AFO)
	options.decomp = 0;			// Template decomposition (0=no, 1=yes)
	//options.alpha = 0.5;		// Weight for bundle size/orthgonal proximity
	options.verbose = false;

	Sapo *sapo = new Sapo(sir,options);
	int reach_steps = 300;
	Flowpipe* flowpipe = sapo->reach(B,reach_steps);	// reachability analysis

	// Store the constructed flowpipe in file sir.m (in Matlab format)
	char file_name[] = "sir_flowpipe.m";
	flowpipe->plotRegionToFile(file_name,'w');

}


