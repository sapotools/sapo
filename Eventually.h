/*
 * Eventually.h
 *
 *  Created on: May 20, 2014
 *      Author: Tommaso Dreossi
 */

#ifndef EVENTUALLY_H_
#define EVENTUALLY_H_

#include "STL.h"

class Eventually : public STL {

private:

	STL * f;
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
