/**
 * @file Until.h
 * Until STL formula
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#ifndef UNTIL_H_
#define UNTIL_H_

#include "STL.h"

class Until : public STL {

private:

	STL * f1, * f2;		// subformulas
	int a, b;			// interval bounds

public:

	Until(STL * f1, int a, int b, STL * f2);

	STL * getLeftSubFormula(){return f1;};
	STL * getRightSubFormula(){return f2;};

	int getA(){return a;};
	int getB(){return b;};

	void setA(int a){this->a = a;};
	void setB(int b){this->b = b;};

	void print();

	virtual ~Until();
};

#endif /* UNTIL_H */
