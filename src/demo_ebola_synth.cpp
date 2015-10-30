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
//#include "LinearSystem.h"
//#include "Bundle.h"
//#include "Sapo.h"
//#include "STL.h"
//#include "Atom.h"
//#include "Always.h"
//
//using namespace std;
//
//int main(int argc,char** argv){
//
//	/////////// THE Model /////////////
//	int dim_sys = 4;
//
//	// List of state variables and parameters
//	symbol s("s"), i("i"), t("t"), r("r"), tau("tau"), dist("dist");
//	lst vars, dyns, params;
//	vars = s, i, t, r;
//	params = tau, dist;
//
//	// Systems's fixed parameters
//	ex sigma1 = 0.1428;
//	ex sigma2 = 0.2;
//	ex epsilon = 0.7;
//	ex delta = 0.00008;
//	ex rho = 0.5;
//	ex g = rho*(1-dist)*(i+epsilon*t);
//
//	// System's dynamics
//	ex ds = s*(1-g);
//	ex di = (1-tau)*(1-sigma1)*(1-delta)*i +s*g;
//	ex dt = (1-sigma2)*t + tau*(1-sigma1)*(1-delta)*i;
//	ex dr = r + sigma1*(1-delta)*i +sigma2*t;
//	dyns = ds,di,dt,dr;
//
//	Model *sitr = new Model(vars,params,dyns);
//
//	// The initial set
//	int num_dirs = 4;
//
//	// Directions matrix
//	vector< double > Li (dim_sys,0);
//	vector< vector< double > > L (num_dirs,Li);
//	L[0][0] = 0.7053; L[0][1] = 0.7053; L[0][2] = 0.7053; L[0][3] = 0;
//	L[1][0] = 0; L[1][1] = 0.9806; L[1][2] = 0.1961; L[1][3] = 0;
//	L[2][0] = 0; L[2][1] = 0; L[2][2] = 1; L[2][3] = 0;
//	L[3][0] = 0; L[3][1] = 0.7071; L[3][2] = 0; L[3][3] = 0.7071;
//
//	// Offsets
//	vector< double > offp (num_dirs,0);
//	vector< double > offm (num_dirs,0);
//	offp[0] = 0.7053; offm[0] = -0.6912;
//	offp[1] = 0.0981; offm[1] = -0.0883;
//	offp[2] = 0.0000; offm[2] = 0.00000;
//	offp[3] = 0.0707; offm[3] = -0.0636;
////	offp[4] = 1; offm[4] = 0;
////	offp[5] = 1; offm[5] = 0;
////	offp[6] = 1; offm[6] = 0;
////	offp[7] = 1; offm[7] = 0;
//
//	// Template matrix
//	vector< int > Ti (dim_sys,0);
//	vector< vector< int > > T (1,Ti);
//
//	T[0][0] = 0; T[0][1] = 1; T[0][2] = 2; T[0][3] = 3;
////	T[1][0] = 1; T[1][1] = 2; T[1][2] = 3;
////	T[2][0] = 2; T[2][1] = 3; T[2][2] = 4;
////	T[3][0] = 5; T[3][1] = 2; T[3][2] = 1;
//
//	// The initial parameter set
//	vector<double> pAi (2,0);
//	vector< vector<double> > pA (4,pAi);
//	vector<double> pb (4,0);
//	pA[0][0] = 1; pA[0][1] = 0; pb[0] = 0.002;
//	pA[1][0] = -1; pA[1][1] = 0; pb[1] = -0.001;
//	pA[2][0] = 0; pA[2][1] = 1; pb[2] = 0.01;
//	pA[3][0] = 0; pA[3][1] = -1; pb[3] = -0.005;
//
//	LinearSystem *parameters = new LinearSystem(pA,pb);
//	LinearSystemSet *parameter_set = new LinearSystemSet(parameters);
//
//	Bundle *B = new Bundle(L,offp,offm,T);
//
//
//	// The specification
//	ex constraint = i-0.4235;
//	Atom *sigma = new Atom(constraint,0);
//	STL *phi = new Always(0,50,sigma);
//
//	// Options for Sapo
//	sapo_opt options;
//	options.trans = 1;
//	options.alpha = 0.5;
//	options.decomp = 0;
//	options.verbose = false;
//
//	Sapo *sapo = new Sapo(sitr,options);
//
//	parameter_set->at(0)->plotRegion();
//
//	//clock_t tStart = clock();
//	//LinearSystemSet  *synthpara = sapo->synthesize(B,parameter_set,phi);
//	//printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
//
////	if(!synthpara->isEmpty()){
//		clock_t tStart = clock();
//		vector<Bundle*> flowpipe = sapo->reach(B,parameter_set->at(0),50);
//		printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
////
////		vector<double> plotsp;
////		vector<double> plotsm;
////		vector<double> plotip;
////		vector<double> plotim;
////		vector<double> plottp;
////		vector<double> plottm;
////		vector<double> plott;
////
////		for(int i=0; i<flowpipe.size(); i++){
////			plotsp.push_back(flowpipe[i]->getOffp(0));
////			plotsm.push_back(-flowpipe[i]->getOffm(0));
////			plotip.push_back(flowpipe[i]->getOffp(1));
////			plotim.push_back(-flowpipe[i]->getOffm(1));
////			plottp.push_back(flowpipe[i]->getOffp(2));
////			plottm.push_back(-flowpipe[i]->getOffm(2));
////			plott.push_back(i);
////		}
////
////		cout<<"t = [";
////		for(int i=0; i<plott.size(); i++){
////			cout<<plott[i]<<" ";
////		}
////		cout<<"]\n";
////		cout<<"sp = [";
////		for(int i=0; i<plotsp.size(); i++){
////			cout<<plotsp[i]<<" ";
////		}
////		cout<<"]\n";
////		cout<<"sm = [";
////		for(int i=0; i<plotsm.size(); i++){
////			cout<<plotsm[i]<<" ";
////		}
////		cout<<"]\n";
////		cout<<"ip = [";
////		for(int i=0; i<plotip.size(); i++){
////			cout<<plotip[i]<<" ";
////		}
////		cout<<"]\n";
////		cout<<"im = [";
////		for(int i=0; i<plotim.size(); i++){
////			cout<<plotim[i]<<" ";
////		}
////		cout<<"]\n";
////		cout<<"tp = [";
////		for(int i=0; i<plottp.size(); i++){
////			cout<<plottp[i]<<" ";
////		}
////		cout<<"]\n";
////		cout<<"tm = [";
////		for(int i=0; i<plotim.size(); i++){
////			cout<<plottm[i]<<" ";
////		}
////		cout<<"]\n";
////	}else{
////		cout<<cout<<"EMPTY";
////	}
//
//	//synthpara->at(0)->plotRegion();
//
//
//}
//
//
