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

	STL * f;			// subformula
	int a, b;			// temporal interval bounds

public:

	Always(int a, int b, STL * f);

	STL * getSubFormula(){return f;};

	int getA(){return a;};
	int getB(){return b;};

	void setA(int a){this->a = a;};
	void setB(int b){this->b = b;};

	void print();

	virtual ~Always();
};

#endif /* ALWAYS_H */
