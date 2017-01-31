/**
 * @file pSIR.cpp
 * Parametric SIR epidemic model
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#include "SIRp.h"

/**
 * Constructor that instantiates a parallelotope
 *
 * @param[in] vars vector with the list of variables used for the generator functions
 * @param[in] u collection of generator versors
 */

 SIRp::SIRp(){


   ///// The dynamical system /////

 	// System dimension (number of variables)
 	int dim_sys = 3;

 	// List of state variables and parameters
 	symbol s("s"), i("i"), r("r");
 	symbol beta("beta"), gamma("gamma");

 	lst vars, params;
 	vars = s, i, r;
 	params = beta, gamma;

 	// System's dynamics
 	ex ds = s - (beta*s*i)*0.1;				// susceptible
 	ex di = i + (beta*s*i - gamma*i)*0.1;	// infected
 	ex dr = r + gamma*i*0.1;				// removed
 	lst dyns;
 	dyns = ds,di,dr;

 	this->vars = vars;
  this->params = params;
  this->dyns = dyns;


 	///// Reachable set representation (single box or parallelotope) /////

 	int num_dirs = 3;		// number of set directions

 	// Directions matrix
 	vector< double > Li (dim_sys,0);
 	vector< vector< double > > L (num_dirs,Li);
 	L[0][0] = 0.7071; L[0][1] = 0.7071;
 	L[1][0] = -0.7071; L[1][1] = 0.7071;
 	L[2][2] = 1;

 	// Template matrix
 	vector< int > Ti (dim_sys,0);
 	vector< vector< int > > T (1,Ti);
 	T[0][0] = 0; T[0][1] = 1; T[0][2] = 2;

 	// Offsets for the set of initial conditions
 	vector< double > offp (num_dirs,0);
 	vector< double > offm (num_dirs,0);
 	offp[0] = 0.7071; offm[0] = -0.6930;
 	offp[1] = -0.4172; offm[1] = 0.4313;
 	offp[2] = 0.0000; offm[2] = 0.00000;

 	Bundle *B = new Bundle(L,offp,offm,T);


 	///// Initial parameter set (polytope) /////

 	vector<double> pAi (2,0);
 	vector< vector<double> > pA (4,pAi);
 	vector<double> pb (4,0);
 	// Template matrix (pA) and offset (pb)
 	pA[0][0] = 1; pA[0][1] = 0; pb[0] = 0.2;
 	pA[1][0] = -1; pA[1][1] = 0; pb[1] = -0.18;
 	pA[2][0] = 0; pA[2][1] = 1; pb[2] = 0.06;
 	pA[3][0] = 0; pA[3][1] = -1; pb[3] = -0.05;
 	LinearSystem *parameters = new LinearSystem(pA,pb);
 	LinearSystemSet *parameter_set = new LinearSystemSet(parameters);

  this->reachSet = B;
  this->paraSet = parameter_set;

  ex constraint = i-0.4405;
	Atom *sigma = new Atom(constraint,0);
	STL *phi = new Always(50,100,sigma);

  this->spec = phi;


 }
