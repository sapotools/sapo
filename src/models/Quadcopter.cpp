/**
 * @file Quadcopter.cpp
 * Quadcopter model
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#include "Quadcopter.h"

 Quadcopter::Quadcopter(){

   ///// The dynamical system /////

   // System dimension (number of variables)
   strcpy(this->name,"Quadcopter");
   int dim_sys = 17;
   // List of state variables and parameters
   symbol pn("pn"), pe("pe"), h("h"), u("u"), v("v"), w("w"), q0v("q0v"), q1v("q1v"), q2v("q2v"), q3v("q3v"), p("p"), q("q"), r("r"), hI("hI"), uI("uI"), vI("vI"), psiI("psiI");
   lst vars, dyns;
   vars = {pn, pe, h, u, v, w, q0v, q1v, q2v, q3v, p, q, r, hI, uI, vI, psiI};

   ex M = 0.0015;
   ex mr = 0.001;
   ex R = 0.020;
   ex l = 0.045;
   ex g = 9.81;
   ex m = M + 4*mr;
   ex Jx = (2*M*pow(R,2))/5 + 2*pow(l,2)*mr;
   ex Jy = (2*M*pow(R,2))/5 + 2*pow(l,2)*mr;
   ex Jz = (2*M*pow(R,2))/5 + 4*pow(l,2)*mr;

   // Reference values
   ex ur = 0;
   ex vr = 0;
   ex psir = 0;
   ex hr = 1;

   // Outputs
   ex phi = 2*q1v;
   ex theta = 2*q2v;
   ex psi = 2*q3v;

   ex delta = 0.01;

   ex F = 0.0361*hI + 0.0694*h + 0.0603*w;
   ex tauphi = -0.0003*vI - 0.0005*v - 0.0018*phi - 0.0004*p;
   ex tautheta = 0.0003*uI + 0.0005*u - 0.0018*theta - 0.0004*q;
   ex taupsi = -0.0003*psiI - 0.0006*psi - 0.0003*r;

   // System's dynamics
   ex dpn = pn + (u*(2*pow(q0v,2) + 2*pow(q1v,2) - 1) - v*(2*q0v*q3v - 2*q1v*q2v ) + w*(2*q0v*q2v + 2*q1v*q3v ))*delta;
   ex dpe = pe + (v*(2*pow(q0v,2) + 2*pow(q2v,2) - 1) + u*(2*q0v*q3v + 2*q1v*q2v ) - w*(2*q0v*q1v - 2*q2v*q3v ))*delta;
   ex dh = h + (w*(2*pow(q0v,2) + 2*pow(q3v,2) - 1) - u*(2*q0v*q2v - 2*q1v*q3v ) + v*(2*q0v*q1v + 2*q2v*q3v ))*delta;

   ex du = u + (r*v - q*w - g*(2*q0v*q2v - 2*q1v*q3v ))*delta;
   ex dv = v + (p*w - r*u + g*(2*q0v*q1v + 2*q2v*q3v ))*delta;
   ex dw = w + (q*u - p*v -F/m + g*(2*pow(q0v,2) + 2*pow(q3v,2) - 1 ))*delta;

   ex dq0v = q0v +(-(q1v/2)*p - (q2v/2)*q - (q3v/2)*r)*delta;
   ex dq1v = q1v + ((q0v/2)*p - (q3v/2)*q + (q2v/2)*r)*delta;
   ex dq2v = q2v + ((q3v/2)*p + (q0v/2)*q - (q1v/2)*r)*delta;
   ex dq3v = q3v + ((q1v/2)*q - (q2v/2)*p + (q0v/2)*r)*delta;

   ex dp = p + ((1/Jx)*tauphi + ((Jy - Jz)/Jx)*q*r)*delta;
   ex dq = q + ((1/Jy)*tautheta - ((Jx - Jz)/Jy)*p*r)*delta;
   ex dr = r + ((1/Jz)*taupsi + ((Jx - Jy)/Jz)*p*q)*delta;

   // Controller
   ex dhI = hI + (h - hr)*delta;
   ex duI = uI +(u - ur)*delta;
   ex dvI = vI + (v - vr)*delta;
   ex dpsiI = psiI + (psi - psir)*delta;

   dyns = {dpn,dpe,dh,du,dv,dw,dq0v,dq1v,dq2v,dq3v,dp,dq,dr,dhI,duI,dvI,dpsiI};

   this->vars = vars;
   this->dyns = dyns;

   ///// Parallelotope bundle for reachable set representation /////
   int num_dirs = 18;
   int num_temp = 2;

   // Directions matrix
   vector< double > Li (dim_sys,0);
   vector< vector< double > > L (num_dirs,Li);
   //initialize box
   for(int i=0; i<dim_sys; i++){
     L[i][i] = 1;
   }
   L[17][2] = 0.5; L[17][5] = 0.5; L[17][6] = 0.5; L[17][15] = 0.25;

   // Template matrix
   vector< int > Ti (dim_sys,0);
   vector< vector< int > > T (num_temp,Ti);
   for(int i=0; i<dim_sys; i++){
     T[0][i] = i;
     T[1][i] = i;
   }
   T[1][5] = 17;

   vector< double > offp (num_dirs,0);
   vector< double > offm (num_dirs,0);

   // Offsets for the set of initial conditions
   offp[2] = 0.21; offm[2] = -0.20;		// h
   offp[6] = 1; offm[6] = -1;				// init quaternion
   offp[17] = 100; offm[17] = 100;		// init quaternion



   Bundle *B = new Bundle(L,offp,offm,T);
   this->reachSet = B;

 }
