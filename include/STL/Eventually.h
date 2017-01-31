/**
 * @file Eventually.h
 * Eventually STL formula
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#ifndef EVENTUALLY_H_
#define EVENTUALLY_H_

#include "STL.h"

class Eventually : public STL {

private:

	STL * f;			// subformula
	int a, b;			// interval bounds

public:

	Eventually(int a, int b, STL * f);

	STL * getSubFormula(){return f;};

	int getA(){return a;};
	int getB(){return b;};

	void setA(int a){this->a = a;};
	void setB(int b){this->b = b;};

	void print();

	virtual ~Eventually();
};

#endif /* Eventually_H */
