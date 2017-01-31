/**
 * @file Influenza.cpp
 * Parametric influenza epidemic model
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#include "Influenza.h"

/**
 * Constructor that instantiates a parallelotope
 *
 * @param[in] vars vector with the list of variables used for the generator functions
 * @param[in] u collection of generator versors
 */

 Influenza::Influenza(){

   int dim_sys = 4;

 	// List of state variables and parameters
 	symbol s("s"), i("i"), t("t"), r("r"), tau("tau"), dist("dist");
 	lst vars, dyns, params;
 	vars = s, i, t, r;
 	params = tau, dist;

 	// Systems's fixed parameters
 	ex sigma1 = 0.1428;
 	ex sigma2 = 0.2;
 	ex epsilon = 0.7;
 	ex delta = 0.00008;
 	ex rho = 0.5;
 	ex g = rho*(1-dist)*(i+epsilon*t);

 	// System's dynamics
 	ex ds = s*(1-g);
 	ex di = (1-tau)*(1-sigma1)*(1-delta)*i +s*g;
 	ex dt = (1-sigma2)*t + tau*(1-sigma1)*(1-delta)*i;
 	ex dr = r + sigma1*(1-delta)*i +sigma2*t;
 	dyns = ds,di,dt,dr;

 	this->vars = vars;
  this->params = params;
  this->dyns = dyns;

 	// The initial set
 	int num_dirs = 4;

 	// Directions matrix
 	vector< double > Li (dim_sys,0);
 	vector< vector< double > > L (num_dirs,Li);
 	L[0][0] = 0.7053; L[0][1] = 0.7053; L[0][2] = 0.7053; L[0][3] = 0;
 	L[1][0] = 0; L[1][1] = 0.9806; L[1][2] = 0.1961; L[1][3] = 0;
 	L[2][0] = 0; L[2][1] = 0; L[2][2] = 1; L[2][3] = 0;
 	L[3][0] = 0; L[3][1] = 0.7071; L[3][2] = 0; L[3][3] = 0.7071;

 	// Offsets
 	vector< double > offp (num_dirs,0);
 	vector< double > offm (num_dirs,0);
 	offp[0] = 0.7053; offm[0] = -0.6912;
 	offp[1] = 0.0981; offm[1] = -0.0883;
 	offp[2] = 0.0000; offm[2] = 0.00000;
 	offp[3] = 0.0707; offm[3] = -0.0636;

 	// Template matrix
 	vector< int > Ti (dim_sys,0);
 	vector< vector< int > > T (1,Ti);

 	T[0][0] = 0; T[0][1] = 1; T[0][2] = 2; T[0][3] = 3;

 	// The initial parameter set
 	vector<double> pAi (2,0);
 	vector< vector<double> > pA (4,pAi);
 	vector<double> pb (4,0);
 	pA[0][0] = 1; pA[0][1] = 0; pb[0] = 0.002;
 	pA[1][0] = -1; pA[1][1] = 0; pb[1] = -0.001;
 	pA[2][0] = 0; pA[2][1] = 1; pb[2] = 0.01;
 	pA[3][0] = 0; pA[3][1] = -1; pb[3] = -0.005;

 	LinearSystem *parameters = new LinearSystem(pA,pb);
 	LinearSystemSet *parameter_set = new LinearSystemSet(parameters);

  Bundle *B = new Bundle(L,offp,offm,T);

  this->reachSet = B;
  this->paraSet = parameter_set;


 	// The specification
 	ex constraint = i-0.4235;
 	Atom *sigma = new Atom(constraint,0);
 	STL *phi = new Always(0,50,sigma);

  this->spec = phi;

 }
