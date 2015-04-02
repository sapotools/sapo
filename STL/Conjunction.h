/*
 * Conjunction.h
 *
 *  Created on: May 20, 2014
 *      Author: Tommaso Dreossi
 */

#ifndef CONJUNCTION_H_
#define CONJUNCTION_H_

#include "STL.h"

class Conjunction : public STL {

private:
	STL * f1, * f2;

public:

	Conjunction(STL * f1, STL * f2);

	STL * getLeftSubFormula(){return f1;};
	STL * getRightSubFormula(){return f2;};

	void print();

	virtual ~Conjunction();
};

#endif /* CONJUNCTION_H_ */
