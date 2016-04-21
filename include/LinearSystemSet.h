/**
 * @file LinearSystemSet.h
 * Represent and manipulate a set of linear systems
 * It can be used to represent a symbolic union of polytopes
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#ifndef LINEARSYSTEMSET_H_
#define LINEARSYSTEMSET_H_

#include "LinearSystem.h"

class LinearSystemSet {

private:

	vector<LinearSystem*> set;	// set of linear systems

public:

	LinearSystemSet();
	LinearSystemSet(LinearSystem *LS);
	LinearSystemSet(vector<LinearSystem*> set);

	vector<LinearSystem*> getSet();
	void add(LinearSystem *LS);

	// operations on set
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
