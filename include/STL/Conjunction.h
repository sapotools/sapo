/**
 * @file Conjunction.h
 * Conjunction STL formula
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#ifndef CONJUNCTION_H_
#define CONJUNCTION_H_

#include "STL.h"

class Conjunction : public STL {

private:
	const STL * f1, * f2;		// subformulas

public:

	// constructor
	Conjunction(const STL * f1, const STL * f2);

	inline const STL * getLeftSubFormula() const {return f1;}
	inline const STL * getRightSubFormula() const {return f2;}

	void print() const;

	~Conjunction();
};

#endif /* CONJUNCTION_H_ */
