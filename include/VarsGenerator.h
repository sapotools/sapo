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
	unsigned int dim;

	lst qs;	// base vertex
	lst as; // alphas (free variables)
	lst bs; // betas (lenghts)
	lst ls; // directions (lambdas)
	vector<lst> us; // genertor versors

public:
	VarsGenerator(const unsigned int dim);

	inline const lst& getBaseVertex() const { return this->qs; };
	inline const lst& getFreeVars() const { return this->as; };
	inline const lst& getLenghts() const { return this->bs; };
	inline const lst& getDirections() const { return this->ls; };
	inline const vector<lst>& getVersors() const { return this->us; };

	LinearSystem genBox(const vector<double>& b) const;

	virtual ~VarsGenerator();
};

#endif /* VARSGENERATOR_H_*/
