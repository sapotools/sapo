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
//	int dim_sys = 2;
//
//	// List of state variables and parameters
//	symbol x("x"), y("y");//, beta("beta"), gamma("gamma"), qr("qr");
//	lst vars, dyns, params;
//	vars = x, y;//, beta, gamma, qr;
//
//	// System's dynamics
//	//ex dx = x + (x-0.33*pow(x,3)-y+0.5)*0.1;
//	//ex dy = y + (x+0.7-0.8*y)*0.1;
//
//	ex dx = x + (0.5*pow(x,2) - 0.5*pow(y,2))*0.01;
//	ex dy = y + (2*x*y)*0.01;
//
//	dyns = dx,dy;//,dbeta, dgamma, dqr;
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
//	//D->decompose(Ab,rnd_directions);
//	//reacher->numericalReach(D,1);
//
//	int num_dirs = 4;
//	int dim  = 2;
//
//	vector< double > Li (dim,0);
//	vector< vector< double > > L (num_dirs,Li);
//	L[0][0] = 1;
//	L[1][1] = 1;
//	L[2][0] = -1; L[2][1] = 1;
//	L[3][0] = 1; L[3][1] = 1;
////	L[4][0] = -0.33; L[4][1] = 0.7;
////	L[5][0] = 0.6; L[5][1] = 0.8;
//
//	vector< double > offp (num_dirs,0);
//	vector< double > offm (num_dirs,0);
////	offp[0] = 1.1; offm[0] = -0.999;
////	offp[1] = 0; offm[1] = 0.001;
////	offp[2] = 2; offm[2] = 20;
////	offp[3] = 2; offm[3] = 11;
////	offp[4] = 2; offm[4] = 11;
////	offp[5] = 2; offm[5] = 11;
//
//	offp[0] = 0.1; offm[0] = -0.05;
//	offp[1] = 1; offm[1] = -0.99;
//	offp[2] = 2; offm[2] = 20;
//	offp[3] = 2; offm[3] = 11;
////	offp[4] = 2; offm[4] = 11;
////	offp[5] = 2; offm[5] = 11;
//
//	vector< int > Ti (dim,0);
//	vector< vector< int > > T (2,Ti);
//
//	T[0][0] = 0; T[0][1] = 2;
//	T[1][0] = 1; T[1][1] = 3;
//	//T[2][0] = 0; T[2][1] = 2;
//	//T[3][0] = 1; T[3][1] = 3;
//	//T[4][0] = 0; T[4][1] = 3;
//
//	vector<lst> paraVars;
//	paraVars.push_back(qs);
//	paraVars.push_back(as);
//	paraVars.push_back(bs);
//
//	Bundle *B = new Bundle(paraVars,L,offp,offm,T);
//
//	B->getBundle()->plotRegionT(0);
//
//	clock_t tStart = clock();
//	int steps = 25;
//	for(int i=0; i<steps; i++){
//
//		pair< vector<double>, vector<double> > newOffsets = B->transform(vars,dyns,true);
//		B->setOffsetP(newOffsets.first);
//		B->setOffsetM(newOffsets.second);
//
//		B->getBundle()->plotRegionT((i+1)*0.01);
//
//		B->setTemplate(B->decomposeTotalRand());
//
//	}
////	B->getBundle()->plotRegionT(steps*0.01);
//
//
//	printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
//
//	cout<<"done";
//
//}
//
//
