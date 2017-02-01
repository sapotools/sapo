/**
 * @file Ebola.cpp
 * Parametric Ebola epidemic model
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#include "Ebola.h"

 Ebola::Ebola(){

   /////////// THE Model /////////////
   strcpy(this->name,"Ebola");
   int dim_sys = 5;

   // List of state variables and parameters
   symbol s("s"), e("e"), q("q"), i("i"), r("r"), kappa1("kappa1"), gamma1("gamma1");
   lst vars, dyns, params;
   vars = {s, e,q, i, r};
   params = {kappa1, gamma1};

   // Systems's fixed parameters
   ex beta = 0.35;
   ex kappa2 = 0.3;
   ex gamma2 = 0.6;
   ex sigma = 0.28;
   ex Delta = 0.5;

   // System's dynamics
   ex ds = s - ((s*beta*i) + gamma1*q)*Delta;
   ex de = e + ((s*beta*i) - (kappa1+kappa2)*e)*Delta;
   ex dq = q + (kappa1*e - (gamma1+gamma2)*q)*Delta;
   ex di = i + (gamma2*q + kappa2*e - sigma*i)*Delta;
   ex dr = r + (sigma*i)*Delta;

   dyns = {ds,de,dq,di,dr};

   this->vars = vars;
   this->params = params;
   this->dyns = dyns;

   // The initial set
   int num_dirs = 5;

   // Directions matrix
   vector< double > Li (dim_sys,0);
   vector< vector< double > > L (num_dirs,Li);
   L[0][0] = 1;
   L[1][1] = 1;
   L[2][2] = 1;
   L[3][3] = 1;
   L[4][4] = 1;

   // Offsets
   vector< double > offp (num_dirs,0);
   vector< double > offm (num_dirs,0);
   offp[0] = 0.8; offm[0] = -0.79;
   offp[1] = 0; offm[1] = 0;
   offp[2] = 0.0000; offm[2] = 0.00000;
   offp[3] = 0.2; offm[3] = -0.19;
   offp[4] = 0; offm[4] = 0;

   // Template matrix
   vector< int > Ti (dim_sys,0);
   vector< vector< int > > T (1,Ti);

   T[0][0] = 0; T[0][1] = 1; T[0][2] = 2; T[0][3] = 3; T[0][4] = 4;

   // The initial parameter set
   vector<double> pAi (2,0);
   vector< vector<double> > pA (4,pAi);
   vector<double> pb (4,0);
   pA[0][0] = 1; pA[0][1] = 0; pb[0] = 0.3;
   pA[1][0] = -1; pA[1][1] = 0; pb[1] = -0.2;
   pA[2][0] = 0; pA[2][1] = 1; pb[2] = 0.5;
   pA[3][0] = 0; pA[3][1] = -1; pb[3] = -0.2;

   LinearSystem *parameters = new LinearSystem(pA,pb);
   LinearSystemSet *parameter_set = new LinearSystemSet(parameters);

   Bundle *B = new Bundle(L,offp,offm,T);

   this->reachSet = B;
   this->paraSet = parameter_set;


   // The specification
 	ex constraint1 = -q+0.04;
 	ex constraint2 = -i+0.27;
 	Atom *phi1 = new Atom(constraint1,0);
 	Atom *phi2 = new Atom(constraint2,1);

 	STL *phi = new Until(phi1, 10, 15, phi2);

   this->spec = phi;

}
