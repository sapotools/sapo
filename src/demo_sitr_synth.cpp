/*
 * main.cpp
 *
 *  Created on: May 20, 2014
 *      Author: Tommaso Dreossi
 */

#include <stdio.h>
#include <iostream>

#include "Common.h"
#include "LinearSystem.h"
#include "Bundle.h"
#include "Sapo.h"
#include "STL.h"
#include "Atom.h"
#include "Until.h"
#include "Always.h"
#include "Conjunction.h"
#include "Disjunction.h"

using namespace std;

int main(int argc,char** argv){

	/////////// THE Model /////////////
	int dim_sys = 5;

	// List of state variables and parameters
	symbol s("s"), e("e"), q("q"), i("i"), r("r"), kappa1("kappa1"), gamma1("gamma1");
	lst vars, dyns, params;
	vars = s, e,q, i, r;
	params = kappa1, gamma1;

	// Systems's fixed parameters
	ex beta = 0.9;
	ex kappa2 = 0.5;
	ex gamma2 = 0.5;
	ex sigma = 0.28;
	ex Delta = 0.5;

	// System's dynamics
	ex ds = s - ((s*beta*i) + gamma1*q)*Delta;
	ex de = e + ((s*beta*i) - (kappa1+kappa2)*e)*Delta;
	ex dq = q + (kappa1*e - (gamma1+gamma2)*q)*Delta;
	ex di = i + (gamma2*q + kappa2*e - sigma*i)*Delta;
	ex dr = r + (sigma*i)*Delta;

	dyns = ds,de,dq,di,dr;

	Model *ebola = new Model(vars,params,dyns);

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
//	offp[5] = 1; offm[5] = 0;
//	offp[6] = 1; offm[6] = 0;
//	offp[7] = 1; offm[7] = 0;

	// Template matrix
	vector< int > Ti (dim_sys,0);
	vector< vector< int > > T (1,Ti);

	T[0][0] = 0; T[0][1] = 1; T[0][2] = 2; T[0][3] = 3; T[0][4] = 4;
//	T[1][0] = 1; T[1][1] = 2; T[1][2] = 3;
//	T[2][0] = 2; T[2][1] = 3; T[2][2] = 4;
//	T[3][0] = 5; T[3][1] = 2; T[3][2] = 1;

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


	// The specification
	ex constraint1 = -q+0.01;
	ex constraint2 = q-0.05;
	ex constraint3 = -i+0.13;
	ex constraint4 = i-0.29;
	Atom *phi1 = new Atom(constraint1,0);
	Atom *phi2 = new Atom(constraint2,1);
	Atom *phi3 = new Atom(constraint3,2);
	Atom *phi4 = new Atom(constraint4,3);
	Conjunction *phi12 = new Conjunction(phi1,phi2);
	Conjunction *phi34 = new Conjunction(phi3,phi4);

	STL *phi = new Until(phi12, 10, 15, phi34);

	// Options for Sapo
	sapo_opt options;
	options.trans = 1;
	options.alpha = 0.5;
	options.decomp = 0;
	options.verbose = false;

	Sapo *sapo = new Sapo(ebola,options);


	clock_t tStart = clock();
	LinearSystemSet  *synthpara = sapo->synthesize(B,parameter_set,phi);
	printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

	cout<<"Size: "<<synthpara->size()<<"\n";


	for(int i=0; i<synthpara->size(); i++){
		synthpara->at(i)->plotRegion();
	}

	if(synthpara->size() != 0){
		//clock_t tStart = clock();
		vector<Bundle*> flowpipe = sapo->reach(B,synthpara->at(0),15);
		//printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

		vector<double> plotsp;
		vector<double> plotsm;
		vector<double> plotip;
		vector<double> plotim;
		vector<double> plotqp;
		vector<double> plotqm;
		vector<double> plotrp;
		vector<double> plotrm;
		vector<double> plott;

		for(int i=0; i<flowpipe.size(); i++){
			plotsp.push_back(flowpipe[i]->getOffp(0));
			plotsm.push_back(-flowpipe[i]->getOffm(0));
			plotqp.push_back(flowpipe[i]->getOffp(2));
			plotqm.push_back(-flowpipe[i]->getOffm(2));
			plotip.push_back(flowpipe[i]->getOffp(3));
			plotim.push_back(-flowpipe[i]->getOffm(3));
			plotrp.push_back(flowpipe[i]->getOffp(4));
			plotrm.push_back(-flowpipe[i]->getOffm(4));
			plott.push_back(i);
		}

		cout<<"t = [";
		for(int i=0; i<plott.size(); i++){
			cout<<plott[i]<<" ";
		}
		cout<<"]\n";
		cout<<"sp = [";
		for(int i=0; i<plotsp.size(); i++){
			cout<<plotsp[i]<<" ";
		}
		cout<<"]\n";
		cout<<"sm = [";
		for(int i=0; i<plotsm.size(); i++){
			cout<<plotsm[i]<<" ";
		}
		cout<<"]\n";
		cout<<"qp = [";
		for(int i=0; i<plotsp.size(); i++){
			cout<<plotqp[i]<<" ";
		}
		cout<<"]\n";
		cout<<"qm = [";
		for(int i=0; i<plotsm.size(); i++){
			cout<<plotqm[i]<<" ";
		}
		cout<<"]\n";
		cout<<"ip = [";
		for(int i=0; i<plotip.size(); i++){
			cout<<plotip[i]<<" ";
		}
		cout<<"]\n";
		cout<<"im = [";
		for(int i=0; i<plotim.size(); i++){
			cout<<plotim[i]<<" ";
		}
		cout<<"]\n";
//		cout<<"tp = [";
//		for(int i=0; i<plottp.size(); i++){
//			cout<<plottp[i]<<" ";
//		}
//		cout<<"]\n";
//		cout<<"tm = [";
//		for(int i=0; i<plotim.size(); i++){
//			cout<<plottm[i]<<" ";
//		}
//		cout<<"]\n";
		cout<<"rp = [";
		for(int i=0; i<plotip.size(); i++){
			cout<<plotrp[i]<<" ";
		}
		cout<<"]\n";
		cout<<"rm = [";
		for(int i=0; i<plotim.size(); i++){
			cout<<plotrm[i]<<" ";
		}
		cout<<"]\n";



	}else{
		cout<<cout<<"EMPTY";
	}
//
//	//synthpara->at(0)->plotRegion();


}


