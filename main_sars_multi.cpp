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
//#include "MultiParallelotope.h"
//
//#include "DynamicalSystem.h"
//#include "DiscreteDynamicalSystem.h"
//
//#include "ParameterSynthesizer.h"
//#include "Bundle.h"
//#include "MultiReacher.h"
//
//#include "STL/STL.h"
//#include "STL/Atom.h"
//#include "STL/Conjunction.h"
//#include "STL/Disjunction.h"
//#include "STL/Until.h"
//
//#include "Models/SI.h"
//
//using namespace std;
//
//int main(int argc,char** argv){
//
//
//
//	/////////// THE Model /////////////
//	int dim = 6;
//
//	// List of state variables and parameters
//	symbol s("s"), e("e"), q("q"), i("i"), j("j"), r("r");
//	symbol p("p"), pi("pi"), beta("beta"), epsE("epsE"), epsQ("epsQ"), epsJ("epsJ"), N("N"), mu("mu"), gamma1("gamma1"), gamma2("gamma2"), kappa1("kappa1"), kappa2("kappa2"), d1("d1"), d2("d2"), sigma1("sigma1"), sigma2("sigma2");
//	symbol h("h");
//	lst vars, dyns;
//	vars = s, e, q, i, j, r;
//
//	// System's dynamics
//	ex ds = s + (pi - (s*(beta*i+epsE*beta*e+epsQ*beta*q+epsJ*beta*j))/N - mu*s)*h;
//	ex de = e + (p + (s*(beta*i+epsE*beta*e+epsQ*beta*q+epsJ*beta*j))/N - (gamma1+kappa1+mu)*e)*h;
//	ex dq = q + (gamma1*e-(kappa2+mu)*q)*h;
//	ex di = i + (kappa1*e-(gamma2+d1+sigma1+mu)*i)*h;
//	ex dj = j + (gamma2*i+kappa2*q-(sigma2+d2+mu)*j)*h;
//	ex dr = r + (sigma1*i+sigma2*j-mu*r)*h;
//
//	lst sub;
//	sub.append(pi == 221);
//	sub.append(p == 0);
//	sub.append(N == 131.5);
//	sub.append(mu == 0.05);
//	sub.append(sigma1 == 0.0337);
//	sub.append(sigma2 == 0.0386);
//	sub.append(d1 == 0.015);
//	sub.append(d2 == 0.0068);
//	sub.append(epsE == 0);
//	sub.append(epsQ == 0);
//	sub.append(epsJ == 0.5);
//	sub.append(beta == 0.15);
//	sub.append(gamma1 == 0.1);
//	sub.append(gamma2 == 2.50);
//	sub.append(kappa1 == 0.1);
//	sub.append(kappa2 == 0.1);
//	sub.append(h == 0.25);
//
//	dyns.append(ds.subs(sub));
//	dyns.append(de.subs(sub));
//	dyns.append(dq.subs(sub));
//	dyns.append(di.subs(sub));
//	dyns.append(dj.subs(sub));
//	dyns.append(dr.subs(sub));
//
//
//	//	// Initialize the template for the rechability sets
//	symbol q1("q1"), q2("q2"), q3("q3"), q4("q4"), q5("q5"), q6("q6");	// base vertex variables
//	symbol a1("a1"), a2("a2"), a3("a3"), a4("a4"), a5("a5"), a6("a6");	// parallelotope variables
//	symbol b1("b1"), b2("b2"), b3("b3"), b4("b4"), b5("b5"), b6("b6");	// amplitude variables
//
//	symbol u11("u11"), u12("u12"), u13("u13"), u14("u14"), u15("u15"), u16("u16");
//	symbol u21("u21"), u22("u22"), u23("u23"), u24("u24"), u25("u25"), u26("u26");
//	symbol u31("u31"), u32("u32"), u33("u33"), u34("u34"), u35("u35"), u36("u36");
//	symbol u41("u41"), u42("u42"), u43("u43"), u44("u44"), u45("u45"), u46("u46");
//	symbol u51("u51"), u52("u52"), u53("u53"), u54("u54"), u55("u55"), u56("u56");
//	symbol u61("u61"), u62("u62"), u63("u63"), u64("u64"), u65("u65"), u66("u66");
//
//	symbol l1("l1"), l2("l2"), l3("l3"), l4("l4"), l5("l5"), l6("l6");	// lambda variables
//	lst qs, as, bs, u1s, u2s, u3s, u4s, u5s, u6s, ls;
//	qs = q1,q2,q3,q4,q5,q6;
//	as = a1,a2,a3,a4,a5,a6;
//	bs = b1,b2,b3,b4,b5,b6;
//	u1s = u11, u12, u13, u14, u15, u16;
//	u2s = u21, u22, u23, u24, u25, u26;
//	u3s = u31, u32, u33, u34, u35, u36;
//	u4s = u41, u42, u43, u44, u45, u46;
//	u5s = u51, u52, u53, u54, u55, u56;
//	u6s = u61, u62, u63, u64, u65, u66;
//
//	ls = l1,l2,l3,l4,l5,l6;
//	vector< lst > us;
//	us.push_back(u1s); us.push_back(u2s); us.push_back(u3s);
//	us.push_back(u4s); us.push_back(u5s); us.push_back(u6s);
//
//	vector<lst> set_vars;
//	set_vars.push_back(qs);
//	set_vars.push_back(as);
//	set_vars.push_back(bs);
//
//	int num_dirs = 6;
//
//	vector< double > Li (dim,0);
//	vector< vector< double > > L (num_dirs,Li);
//	L[0][0] = 1;
//	L[1][1] = 1;
//	L[2][2] = 1;
//	L[3][3] = 1;
//	L[4][4] = 1;
//	L[5][5] = 1;
//
//	vector< double > offp (num_dirs,0);
//	vector< double > offm (num_dirs,0);
//	offp[0] = 6.5; 			offm[0] = -6.49;		//S
//	offp[1] = 124; 			offm[1] = -123.99;		//E
//	offp[2] = 0.0001; 		offm[2] = -0.000099;	//Q
//	offp[3] = 1; 			offm[3] = -0.99;			//I
//	offp[4] = 0.0001;; 			offm[4] = -0.000099;			//J
//	offp[5] = 0.0001;; 			offm[5] = -0.000099;			//R
//
//	vector< int > Ti (dim,0);
//	vector< vector< int > > T (1,Ti);
//
//	T[0][0] = 0; T[0][1] = 1; T[0][2] = 2; T[0][3] = 3; T[0][4] = 4; T[0][5] = 5;
//	//T[1][0] = 1; T[1][1] = 2; T[1][2] = 3;
//	//T[2][0] = 2; T[2][1] = 3; T[2][2] = 4;
//	//T[3][0] = 5; T[3][1] = 2; T[3][2] = 1;
//
//	vector<lst> paraVars;
//	paraVars.push_back(qs);
//	paraVars.push_back(as);
//	paraVars.push_back(bs);
//
//	Bundle *B = new Bundle(paraVars,L,offp,offm,T);
//
//	vector<int> projx (6,0); vector<int> projy (3,0);
//	projx[0] = 0; projy[0] = 0;
//	projx[1] = 1; projy[1] = 1;
//	projx[2] = 2; projy[2] = 2;
//	projx[3] = 6;
//	projx[4] = 7;
//	projx[5] = 8;
//
//	clock_t tStart = clock();
//	for(int i=0; i<200; i++){
//
//		B->getBundle()->plotRegion(projx,projy);
//
//		//vector< vector<int> > T = B->decomposeRand();
//		//B->setTemplate(T);
//
//		pair< vector<double>, vector<double> > newOffsets = B->transform(vars,dyns);
//		B->setOffsetP(newOffsets.first);
//		B->setOffsetM(newOffsets.second);
//	}
//	//B->getBundle()->plotRegion();
//
//	printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
//
//	cout<<"done";
//
//}
