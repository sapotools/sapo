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
//#include "Parallelotope.h"
//
//#include "DynamicalSystem.h"
//#include "DiscreteDynamicalSystem.h"
//
//#include "ParameterSynthesizer.h"
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
//	Model *M = new Model();
//
//
//	symbol x1("x1"),x2("x2"),x3("x3"),x4("x4"),x5("x5"),x6("x6"), x7("x7"), x8("x8"), x9("x9"), x10("x10"), x11("x11"), x12("x12"), x13("x13"), x14("x14"), x15("x15");
//	lst vars;
//
//	vars = x1,x2,x3,x4,x5,x6,x7;
//	int max_deg = 10;
//
//	ex poly = M->genPoly(vars,max_deg);
//
//	// System's dynamics
//
//	BaseConverter *BC = new BaseConverter(vars,poly);
//	clock_t tStart_1 = clock();
//	//cout<<BC->getBernCoeffs().nops();
//	cout<<BC->getBernCoeffsMatrix().nops();
//	double tEnd_1 = (double)(clock() - tStart_1)/CLOCKS_PER_SEC;
//	printf("\nTime taken: %.2fs\n\n", tEnd_1);
//
//	cout<<"\ndone";
//
////	vector< double > Ai (3,0);
////	vector< vector< double > > A (6,Ai);
////	vector< double > b (6,0);
////
////	A[0][0] = -1;					b[0] = -0.8;
////	A[1][0] = 1;					b[1] = 0.85;
////	A[2][0] = -1; A[2][1] = -1;		b[2] = -0.95;
////	A[3][0] = 1; A[3][1] = 1;		b[3] = 1.00;
////	A[4][2] = -1;					b[4] = 0.0;
////	A[5][2] = 1;					b[5] = 0.3;
////
////	LinearSystem* LS = new LinearSystem(A,b);
////	LS->print();
////
////	symbol q1("q1"), q2("q2"), q3("q3"), q4("q4"), q5("q5");	// base vertex variables
////		symbol a1("a1"), a2("a2"), a3("a3"), a4("a4"), a5("a5");	// parallelotope variables
////		symbol b1("b1"), b2("b2"), b3("b3"), b4("b4"), b5("b5");	// amplitude variables
////		lst qs, as, bs;
////		qs = q1,q2,q3;
////		as = a1,a2,a3;
////		bs = b1,b2,b3;
////
////		vector<lst> set_vars;
////		set_vars.push_back(qs);
////		set_vars.push_back(as);
////		set_vars.push_back(bs);
////
////		vector<double> ui (3,0);
////		vector< vector<double> > u (3,ui);
////
////		u[0][0] = 0.7071; u[0][1] = -0.7071;
////		u[1][1] = 1;
////		u[2][2] = 1;
////
////
////		Parallelotope *P = new Parallelotope(set_vars,u);
////		cout<<P->getGeneratorFunction();
////		vector< vector< double > > temp = P->getTemplate();
////		for(int i=0; i<temp.size(); i++){
////			for(int j=0; j<temp[i].size(); j++){
////				cout<<temp[i][j]<<" ";
////			}
////			cout<<"\n";
////		}
////
////	cout<<"done";
//
//}
