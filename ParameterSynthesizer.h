/*
 * ParameterSynthesizer.h
 *
 *  Created on: Nov 5, 2014
 *      Author: Tommaso Dreossi
 */

#ifndef PARAMETERSYNTHESIZER_H_
#define PARAMETERSYNTHESIZER_H_

#include "Common.h"
#include "DiscreteDynamicalSystem.h"
#include "BaseConverter.h"
#include "Polyhedron.h"
#include "Parallelotope.h"
#include "STL/STL.h"
#include "LinearSystemSet.h"

#include <iostream>
#include <fstream>

class ParameterSynthesizer {

private:

	DiscreteDynamicalSystem *dynamicalSystem;
	ex constraint; 	// constraint to impose over the system
	STL *stl_constraint;
	vector< vector< lst > > dynamicalSystemControlPts;
	lst constraintControlPts;
	synthesizer_opt options;


	STL* initConstraintControlPts(STL *formula);

	LinearSystemSet* refineParameters(vector< poly_values > reach_sets, LinearSystemSet *parameterSet, vector< lst > constraintControlPts);
	LinearSystemSet* synthesizeUntil(vector<poly_values> reach_sets, LinearSystemSet *parameterSet, STL *formula);
	LinearSystemSet* synthesizeAlways(vector<poly_values> reach_sets, LinearSystemSet *parameterSet, STL *formula);
	LinearSystemSet* synthesizeEventually(vector<poly_values> reach_sets, LinearSystemSet *parameterSet, STL *formula);

	vector<poly_values> reachStep(vector<poly_values> reach_sets, LinearSystem *parameterSet, char color);

	double maxVec( vector<double> vec);

public:
	ParameterSynthesizer(DiscreteDynamicalSystem *dynSys, STL *stl_constraint, synthesizer_opt options);
	LinearSystemSet* synthesizeSTL(poly_values reach_set, LinearSystemSet *parameterSet);
	LinearSystemSet* synthesizeSTL(vector<poly_values> reach_sets, LinearSystemSet *parameterSet);
	LinearSystemSet* synthesize(vector<poly_values> reach_sets, LinearSystemSet *parameterSet, STL *formula);
	vector<poly_values> reach(vector<poly_values> reach_sets, LinearSystem *parameterSet, int k );

	virtual ~ParameterSynthesizer();
};

#endif /* PARAMETERSYNTHESIZER_H_ */
