/*
 * LinearSystemSet.h
 *
 *  Created on: Nov 17, 2014
 *      Author: dreossi
 */

#ifndef LINEARSYSTEMSET_H_
#define LINEARSYSTEMSET_H_

#include "LinearSystem.h"

class LinearSystemSet {

private:

	vector<LinearSystem*> set;

public:

	LinearSystemSet();
	LinearSystemSet(LinearSystem *LS);
	LinearSystemSet(vector<LinearSystem*> set);

	vector<LinearSystem*> getSet();
	void add(LinearSystem *LS);

	LinearSystemSet* intersectWith(LinearSystemSet *LSset);
	LinearSystemSet* unionWith(LinearSystemSet *LSset);
	LinearSystemSet* boundedUnionWith(LinearSystemSet *LSset, int bound);

	double boundingVol();
	int size();
	LinearSystem* at(int i);
	bool isEmpty();
	void print();

	virtual ~LinearSystemSet();
};

#endif /* LINEARSYSTEMSET_H_ */
