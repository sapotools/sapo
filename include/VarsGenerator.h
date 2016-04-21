/**
 * @file VarsGenerator.h
 * Automatically generate variables for paralleltope generator functions.
 * For high dimensions declaring manually the variables can be tedious...
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
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
