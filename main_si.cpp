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
//	// SI model
//	Model *model = new SI();
//
//	lst vars, params, dyns;
//	vars = model->getVars();
//	params = model->getParams();
//	dyns = model->getDyns();
//
//	cout<<"Model name: "<<model->getName()<<"\n";
//	cout<<"Dynamics:\n";
//	for(int i=0; i<(signed)dyns.nops(); i++){
//		cout<<"\td"<<vars[i]<<" = "<<dyns[i]<<"\n";
//	}
//	cout<<"Specification:";
//	model->getSpec()->print();
//	cout<<"\n\n";
//
//
//	// Create a new dynamic system with model's dynamics
//	bool rational_dyns = false;
//	DynamicalSystem *cont_SI = new DynamicalSystem(vars,params,dyns,rational_dyns);
//
//	vector<Polyhedron*> multi_template = model->getReachSets();
//	DiscreteDynamicalSystem *disc_SI = new DiscreteDynamicalSystem(vars,params,cont_SI->eulerDisc(0.1), multi_template, false);
//
//
//	// Create the parameter synthesizer
//	synthesizer_opt options;
//	options.largest_para_set = true;
//
//	ParameterSynthesizer *PS = new ParameterSynthesizer(disc_SI,model->getSpec(),options);
//
//	//PS->reach(model->getInitConds(),model->getInitParaSet()->at(0),1);
//
//	// Synthesize
//	cout<<"\nStart synthesis...\n\n";
//	clock_t tStart = clock();
//
//	LinearSystemSet *result = PS->synthesizeSTL(model->getInitConds(),model->getInitParaSet());
//
//	printf("\nTime taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
//
//	if(result->isEmpty()){
//		cout<<"Result: Empty parameter set.";
//	}else{
//		cout<<"Result: Parameter set composed by "<<result->size()<<" convex polytopes. ";
//		cout<<"Print result? (y/n): ";
//		char print_res;
//		cin >> print_res;
//		if(print_res == 'y'){
//			result->print();
//		}
//	}
//
//	cout<<"done";
//
//}
//
