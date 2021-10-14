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

#define MINIMIZE_LS_SET_REPRESENTATION true


class LinearSystemSet {

private:

	vector<LinearSystem*> set;	// set of linear systems

public:

	LinearSystemSet();
	LinearSystemSet(LinearSystem *LS);
	LinearSystemSet(const vector<LinearSystem*>& set);

	const vector<LinearSystem*>& getSet() const;
	void add(LinearSystem *LS);

	// operations on set
	LinearSystemSet* intersectWith(const LinearSystemSet *LSset) const;
	LinearSystemSet* unionWith(LinearSystemSet *LSset);
	LinearSystemSet* boundedUnionWith(LinearSystemSet *LSset, int bound);

	double boundingVol() const;

	/**
	 * Get the size of this set, i.e,
	 * the number of linear systems
	 *
	 * @returns number of linear systems in the set
	 */
	inline int size() const{ return this->set.size(); }

	LinearSystem* at(int i);
	const LinearSystem* at(int i) const;

	bool isEmpty() const;
	void print() const;

	~LinearSystemSet();
};

#endif /* LINEARSYSTEMSET_H_ */
