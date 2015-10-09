///*
// * main.cpp
// *
// *  Created on: May 20, 2014
// *      Author: Tommaso Dreossi
// */
//
//#include <stdio.h>
//#include <iostream>
//
//#include "Common.h"
//#include "BaseConverter.h"
//#include "LinearSystem.h"
//#include "Polyhedron.h"
//#include "Box.h"
//
//#include "MultiParallelotope.h"
//#include "Bundle.h"
//
//#include "DynamicalSystem.h"
//#include "DiscreteDynamicalSystem.h"
//
//#include "ParameterSynthesizer.h"
//#include "MultiReacher.h"
//
//#include "VarsGenerator.h"
//
//#include <Eigen/Dense>
//using namespace Eigen;
//
//using namespace std;
//
//int main(int argc,char** argv){
//
//
//
//	/////////// THE Model /////////////
//
//	int dim_sys = 17;
//	lst vars, dyns;
//
//	// List of state variables and parameters
//	symbol pn("pn"), pe("pe"), h("h"), u("u"), v("v"), w("w"), q0v("q0v"), q1v("q1v"), q2v("q2v"), q3v("q3v"), p("p"), q("q"), r("r"), hI("hI"), uI("uI"), vI("vI"), psiI("psiI");
//
//	ex M = 0.0015;
//	ex mr = 0.001;
//	ex R = 0.020;
//	ex l = 0.045;
//	ex g = 9.81;
//	ex m = M + 4*mr;
//	ex Jx = (2*M*pow(R,2))/5 + 2*pow(l,2)*mr;
//	ex Jy = (2*M*pow(R,2))/5 + 2*pow(l,2)*mr;
//	ex Jz = (2*M*pow(R,2))/5 + 4*pow(l,2)*mr;
//
//	// Reference values
//	ex ur = 0;
//	ex vr = 0;
//	ex psir = 0;
//	ex hr = 1;
//
//	// Outputs
//	ex phi = 2*q1v;
//	ex theta = 2*q2v;
//	ex psi = 2*q3v;
//
//	ex delta = 0.01;
//
//	vars = pn, pe, h, u, v, w, q0v, q1v, q2v, q3v, p, q, r, hI, uI, vI, psiI;
//
//	ex F = 0.0361*hI + 0.0694*h + 0.0603*w;
//	ex tauphi = -0.0003*vI - 0.0005*v - 0.0018*phi - 0.0004*p;
//	ex tautheta = 0.0003*uI + 0.0005*u - 0.0018*theta - 0.0004*q;
//	ex taupsi = -0.0003*psiI - 0.0006*psi - 0.0003*r;
//
//	// System's dynamics
//	ex dpn = pn + (u*(2*pow(q0v,2) + 2*pow(q1v,2) - 1) - v*(2*q0v*q3v - 2*q1v*q2v ) + w*(2*q0v*q2v + 2*q1v*q3v ))*delta;
//	ex dpe = pe + (v*(2*pow(q0v,2) + 2*pow(q2v,2) - 1) + u*(2*q0v*q3v + 2*q1v*q2v ) - w*(2*q0v*q1v - 2*q2v*q3v ))*delta;
//	ex dh = h + (w*(2*pow(q0v,2) + 2*pow(q3v,2) - 1) - u*(2*q0v*q2v - 2*q1v*q3v ) + v*(2*q0v*q1v + 2*q2v*q3v ))*delta;
//
//	ex du = u + (r*v - q*w - g*(2*q0v*q2v - 2*q1v*q3v ))*delta;
//	ex dv = v + (p*w - r*u + g*(2*q0v*q1v + 2*q2v*q3v ))*delta;
//	ex dw = w + (q*u - p*v -F/m + g*(2*pow(q0v,2) + 2*pow(q3v,2) - 1 ))*delta;
//
//	ex dq0v = q0v +(-(q1v/2)*p - (q2v/2)*q - (q3v/2)*r)*delta;
//	ex dq1v = q1v + ((q0v/2)*p - (q3v/2)*q + (q2v/2)*r)*delta;
//	ex dq2v = q2v + ((q3v/2)*p + (q0v/2)*q - (q1v/2)*r)*delta;
//	ex dq3v = q3v + ((q1v/2)*q - (q2v/2)*p + (q0v/2)*r)*delta;
//
//	ex dp = p + ((1/Jx)*tauphi + ((Jy - Jz)/Jx)*q*r)*delta;
//	ex dq = q + ((1/Jy)*tautheta - ((Jx - Jz)/Jy)*p*r)*delta;
//	ex dr = r + ((1/Jz)*taupsi + ((Jx - Jy)/Jz)*p*q)*delta;
//
//	// Controller
//	ex dhI = hI + (h - hr)*delta;
//	ex duI = uI +(u - ur)*delta;
//	ex dvI = vI + (v - vr)*delta;
//	ex dpsiI = psiI + (psi - psir)*delta;
//
//	dyns = dpn,dpe,dh,du,dv,dw,dq0v,dq1v,dq2v,dq3v,dp,dq,dr,dhI,duI,dvI,dpsiI;
//
//	VarsGenerator *varsGen = new VarsGenerator(dim_sys);
//	lst qs, as, bs, ls;
//	vector<lst> us;
//
//	qs = varsGen->getBaseVertex();
//	as = varsGen->getFreeVars();
//	bs = varsGen->getLenghts();
//	ls = varsGen->getDirections();
//	us = varsGen->getVersors();
//
//	int num_dirs = 18;
//
//	vector< double > Li (dim_sys,0);
//	vector< vector< double > > L (num_dirs,Li);
//	//initialize box
//	for(int i=0; i<dim_sys; i++){
//		L[i][i] = 1;
//	}
//
//	L[17][2] = 0.5; L[17][5] = 0.5; L[17][6] = 0.5; L[17][15] = 0.25;
//
//	vector< double > offp (num_dirs,0);
//	vector< double > offm (num_dirs,0);
//
//	offp[2] = 0.21; offm[2] = -0.20;		// h
//	offp[6] = 1; offm[6] = -1;		// init quaternion
//	//offp[7] = 22.5; offm[7] = -22.5;		// init quaternion
//	//offp[8] = 50; offm[8] = -50;		// init quaternion
//	//offp[9] = 30; offm[9] = -30;		// init quaternion
//
//	offp[17] = 100; offm[17] = 100;		// init quaternion
////	offp[18] = 1; offm[18] = 0;		// init quaternion
//
//
//
//	vector< int > Ti (dim_sys,0);
//	vector< vector< int > > T (2,Ti);
//	for(int i=0; i<dim_sys; i++){
//		T[0][i] = i;
//		T[1][i] = i;
//	}
//
//	T[1][5] = 17;
////	T[2][2] = 18;
//
//
//	vector<lst> paraVars;
//	paraVars.push_back(qs);
//	paraVars.push_back(as);
//	paraVars.push_back(bs);
//
//	Bundle *B = new Bundle(paraVars,L,offp,offm,T);
//
//	vector<double> plothp;
//	vector<double> plothm;
//	vector<double> plotwp;
//	vector<double> plotwm;
//	vector<double> plott;
//	vector<double> plothip;
//	vector<double> plothim;
//
//	plothp.push_back(B->getOffp(2));
//	plothm.push_back(-B->getOffm(2));
//	plotwp.push_back(B->getOffp(5));
//	plotwm.push_back(-B->getOffm(5));
//	plothip.push_back(B->getOffp(13));
//	plothim.push_back(-B->getOffm(13));
//	plott.push_back(0);
//
//	clock_t tStart = clock();
//	for(int i=0; i<300; i++){
//
//		//B->getBundle()->print();
//
//		pair< vector<double>, vector<double> > newOffsets = B->transform(vars,dyns,true);
//		B->setOffsetP(newOffsets.first);
//		B->setOffsetM(newOffsets.second);
//
//		plothp.push_back(B->getOffp(2));
//		plothm.push_back(-B->getOffm(2));
//		plotwp.push_back(B->getOffp(5));
//		plotwm.push_back(-B->getOffm(5));
//		plothip.push_back(B->getOffp(13));
//		plothim.push_back(-B->getOffm(13));
//
//		plott.push_back((i+1)*0.01);
//
//	}
//
//	printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
//
//	cout<<"t = [";
//	for(int i=0; i<plott.size(); i++){
//		cout<<plott[i]<<" ";
//	}
//	cout<<"]\n";
//	cout<<"hp = [";
//	for(int i=0; i<plothp.size(); i++){
//		cout<<plothp[i]<<" ";
//	}
//	cout<<"]\n";
//	cout<<"hm = [";
//	for(int i=0; i<plothm.size(); i++){
//		cout<<plothm[i]<<" ";
//	}
//	cout<<"]\n";
//	cout<<"wp = [";
//	for(int i=0; i<plotwp.size(); i++){
//		cout<<plotwp[i]<<" ";
//	}
//	cout<<"]\n";
//	cout<<"wm = [";
//	for(int i=0; i<plotwm.size(); i++){
//		cout<<plotwm[i]<<" ";
//	}
//	cout<<"]\n";
//	cout<<"hip = [";
//	for(int i=0; i<plotwp.size(); i++){
//		cout<<plothip[i]<<" ";
//	}
//	cout<<"]\n";
//	cout<<"him = [";
//	for(int i=0; i<plotwm.size(); i++){
//		cout<<plothim[i]<<" ";
//	}
//	cout<<"]\n";
//
//
//}
//
//
