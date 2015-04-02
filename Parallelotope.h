/*
 * Parallelotope.h
 *
 *  Created on: Oct 30, 2014
 *      Author: Tommaso Dreossi
 */

#ifndef PARALLELOTOPE_H_
#define PARALLELOTOPE_H_

#include "Common.h"
#include "Polyhedron.h"
#include "LinearSystem.h"

class Parallelotope : public Polyhedron{

private:

	vector<double> hyperplaneThroughPts(vector< vector<double> > pts);	// find hyper plane passing through pts
	vector<double> lst2vec(ex list);									// convert a lst to a vector of doubles
	double euclidNorm(vector<double> v);								// compute euclidean norm

public:

	Parallelotope(vector<lst> vars, vector< vector<double> > u);

	// Representation conversion
	LinearSystem* gen2const(vector<double> q, vector<double> beta);		// from generator to constraints
	poly_values const2gen(LinearSystem *LS);							// from constraints to generators

	virtual ~Parallelotope();
};

#endif /* PARALLELOTOPE_H_ */
