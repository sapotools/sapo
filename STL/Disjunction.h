/*
 * Disjunction.h
 *
 *  Created on: May 20, 2014
 *      Author: Tommaso Dreossi
 */

#ifndef DISJUNCTION_H_
#define DISJUNCTION_H_

#include "STL.h"

class Disjunction : public STL {

private:
	STL * f1, * f2;

public:

	Disjunction(STL * f1, STL * f2);

	STL * getLeftSubFormula(){return f1;};
	STL * getRightSubFormula(){return f2;};

	void print();

	virtual ~Disjunction();
};

#endif /* CONJUNCTION_H_ */
