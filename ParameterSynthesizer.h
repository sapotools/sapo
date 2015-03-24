/*
 * ParameterSynthesizer.h
 *
 *  Created on: Nov 5, 2014
 *      Author: dreossi
 */

#ifndef PARAMETERSYNTHESIZER_H_
#define PARAMETERSYNTHESIZER_H_

#include "Common.h"
#include "DiscreteDynamicalSystem.h"
#include "BaseConverter.h"
#include "Polyhedron.h"
#include "STL.h"
#include "LinearSystemSet.h"

class ParameterSynthesizer {

private:

	DiscreteDynamicalSystem *dynamicalSystem;
	ex constraint; 	// constraint to impose over the system
	STL *stl_constraint;
	vector< lst > dynamicalSystemControlPts;
	lst constraintControlPts;

	LinearSystem* refineParameters(vector< double > base_v, vector< double > lenghts, LinearSystem *parameterSet);
	LinearSystemSet* refineParameters(vector< double > base_v, vector< double > lenghts, LinearSystemSet *parameterSet, lst constraintControlPts);
	poly_values reachStep(vector< double > base_v, vector< double > lenghts, LinearSystem *parameterSet);
	double maxVec( vector<double> vec);

	STL* initConstraintControlPts(STL *formula);
	LinearSystemSet* synthesizeUntil(vector< double > base_v, vector< double > lenghts, LinearSystemSet *parameterSet, STL *formula);
	LinearSystemSet* synthesizeAlways(vector< double > base_v, vector< double > lenghts, LinearSystemSet *parameterSet, STL *formula);
	LinearSystemSet* synthesizeEventually(vector< double > base_v, vector< double > lenghts, LinearSystemSet *parameterSet, STL *formula);

public:
	ParameterSynthesizer(DiscreteDynamicalSystem *dynSys, STL *stl_constraint);
	LinearSystemSet* synthesizeSTL(vector< double > base_v, vector< double > lenghts, LinearSystemSet *parameterSet);
	LinearSystemSet* synthesize(vector< double > base_v, vector< double > lenghts, LinearSystemSet *parameterSet, STL *formula);

	virtual ~ParameterSynthesizer();
};

#endif /* PARAMETERSYNTHESIZER_H_ */
