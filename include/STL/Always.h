/**
 * @file Always.h
 * Always STL formula
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#ifndef ALWAYS_H_
#define ALWAYS_H_

#include "STL.h"

class Always : public STL {

private:

	const STL * f;			// subformula
	int a, b;			// temporal interval bounds

public:

	Always(int a, int b, const STL * f);

	inline const STL * getSubFormula() const { return this->f; }

	inline const int getA() const { return this->a;}
	inline const int getB() const { return this->b; }

	inline void setA(const int a){ this->a = a; };
	inline void setB(const int b){ this->b = b; };

	void print() const;

	virtual ~Always();
};

#endif /* ALWAYS_H */
