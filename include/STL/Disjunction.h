/**
 * @file Disjunction.h
 * Disjunction STL formula
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#ifndef DISJUNCTION_H_
#define DISJUNCTION_H_

#include "STL.h"

class Disjunction : public STL {

private:
	const STL * f1, * f2;	// subformulas

public:

	Disjunction(const STL * f1, const STL * f2);

	inline const STL * getLeftSubFormula() const {return f1;}
	inline const STL * getRightSubFormula() const {return f2;}

	void print() const;

	~Disjunction();
};

#endif /* CONJUNCTION_H_ */
