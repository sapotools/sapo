/*
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
//#include "Bundle.h"
//#include "Sapo.h"
//
//using namespace std;
//
//int main(int argc,char** argv){
//
//	/////////// THE Model /////////////
//	int dim_sys = 5;
//
//	// List of state variables and parameters
//	symbol x("x"), y1("y1"), y2("y2"), z1("z1"), z2("z2");
//	lst vars, dyns, params;
//	vars = x, y1, y2, z1, z2;
//
//	ex beta1 = 0.001;
//	ex beta2 = 0.001;
//	ex gamma = 0.3;
//	ex delta = 0.5;
//	ex alpha = 0.7;
//	ex h = 0.01;
//
//	// System's dynamics
//	ex dx = x + h*(-beta1*x*y1-beta2*x*y2);
//	ex dy1 = y1 + h*(beta1*x*y1-gamma*y1+delta*beta1*y1*z1+alpha*beta1*y1*z2);
//	ex dy2 = y2 + h*(beta2*x*y2-gamma*y2+delta*beta2*y2*z2+alpha*beta2*y2*z1);
//	ex dz1 = z1 + h*(gamma*y1-delta*beta1*y1*z1-alpha*beta2*y2*z1);
//	ex dz2 = z2 + h*(gamma*y2-delta*beta2*y2*z2-alpha*beta1*y1*z2);
//	dyns = dx,dy1,dy2,dz1,dz2;
//
//	Model *sir = new Model(vars,dyns);
//
//
//	// The initial set
//	int num_dirs = 7;
//	int num_temps = 3;
//
//	// Directions matrix
//	vector< double > Li (dim_sys,0);
//	vector< vector< double > > L (num_dirs,Li);
//	L[0][0] = 1;
//	L[1][1] = 1;
//	L[2][2] = 1;
//	L[3][3] = 1;
//	L[4][4] = 1;
//	L[5][0] = 1; L[5][1] = 0.5;
//	L[6][1] = 0.5; L[6][4] = 1;
//	//L[7][2] = 0.5; L[7][3] = 0.5;
//
//	// Offsets
//	vector< double > offp (num_dirs,0);
//	vector< double > offm (num_dirs,0);
//	offp[0] = 500; offm[0] = -500;
//	offp[1] = 400; offm[1] = -390;
//	offp[2] = 100; offm[2] = -90;
//	offp[3] = 0; offm[3] = 0;
//	offp[4] = 0; offm[4] = 0;
//	offp[5] = 1000; offm[5] = 1000;
//	offp[6] = 1000; offm[6] = 1000;
//	//offp[7] = 1000; offm[7] = 1000;
//
//	// Template matrix
//	vector< int > Ti (dim_sys,0);
//	vector< vector< int > > T (num_temps,Ti);
//
//	T[0][0] = 0; T[0][1] = 1; T[0][2] = 2; T[0][3] = 3; T[0][4] = 4;
//	T[1][0] = 1; T[1][1] = 2; T[1][2] = 3; T[1][3] = 4; T[1][4] = 5;
//	T[2][0] = 0; T[2][1] = 2; T[2][2] = 3; T[2][3] = 5; T[2][4] = 6;
//
//
//	Bundle *B = new Bundle(L,offp,offm,T);
//
//	// Options for Sapo
//	sapo_opt options;
//	options.trans = 1;
//	options.alpha = 0.5;
//	options.decomp = 0;
//	options.verbose = false;
//
//	Sapo *sapo = new Sapo(sir,options);
//
//	clock_t tStart = clock();
//	vector<Bundle*> flowpipe = sapo->reach(B,1500);
//	printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
//	cout<<"\ndone\n";
//
//		vector<double> ploty1p;
//		vector<double> ploty1m;
//		vector<double> ploty2p;
//		vector<double> ploty2m;
//		vector<double> plott;
//
//		for(int i=0; i<flowpipe.size(); i++){
//			ploty1p.push_back(flowpipe[i]->getOffp(1));
//			ploty1m.push_back(-flowpipe[i]->getOffm(1));
//			ploty2p.push_back(flowpipe[i]->getOffp(2));
//			ploty2m.push_back(-flowpipe[i]->getOffm(2));
//			plott.push_back(i*0.01);
//		}
//
//		cout<<"t = [";
//		for(int i=0; i<plott.size(); i++){
//			cout<<plott[i]<<" ";
//		}
//		cout<<"]\n";
//		cout<<"y1p = [";
//		for(int i=0; i<ploty1p.size(); i++){
//			cout<<ploty1p[i]<<" ";
//		}
//		cout<<"]\n";
//		cout<<"y1m = [";
//		for(int i=0; i<ploty1m.size(); i++){
//			cout<<ploty1m[i]<<" ";
//		}
//		cout<<"]\n";
//		cout<<"y2p = [";
//		for(int i=0; i<ploty2p.size(); i++){
//			cout<<ploty2p[i]<<" ";
//		}
//		cout<<"]\n";
//		cout<<"y2m = [";
//		for(int i=0; i<ploty2m.size(); i++){
//			cout<<ploty2m[i]<<" ";
//		}
//		cout<<"]\n";
//
//
//}
//
//
