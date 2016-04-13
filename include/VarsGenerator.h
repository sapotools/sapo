/*
 * VarsGenerator.h
 *
 *  Created on: Jan 7, 2015
 *      Author: Tommaso Dreossi
 */

#ifndef VARSGENERATOR_H_
#define VARSGENERATOR_H_

#include "Common.h"
#include "LinearSystem.h"

class VarsGenerator {

private:

	int dim;

	lst qs;	// base vertex
	lst as; // alphas (free variables)
	lst bs; // betas (lenghts)
	lst ls; // directions (lambdas)
	vector<lst> us; // genertor versors

public:
	VarsGenerator(int dim);

	lst getBaseVertex(){ return this->qs; };
	lst getFreeVars(){ return this->as; };
	lst getLenghts(){ return this->bs; };
	lst getDirections(){ return this->ls; };
	vector<lst> getVersors(){ return this->us; };

	LinearSystem* genBox(vector<double> b);

	virtual ~VarsGenerator();
};

#endif /* VARSGENERATOR_H_*/
