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

	std::vector<LinearSystem*> set;	// set of linear systems

public:
	LinearSystemSet();
	LinearSystemSet(const LinearSystem& ls);
	LinearSystemSet(LinearSystem* ls);
	LinearSystemSet(const std::vector<LinearSystem*>& set);

	void add(LinearSystem *LS);

	LinearSystemSet* get_a_finer_covering() const;

	// operations on set
	LinearSystemSet* getIntersectionWith(const LinearSystemSet *LSset) const;

	// in-place set operations
	LinearSystemSet& unionWith(LinearSystemSet *LSset);
	LinearSystemSet& boundedUnionWith(LinearSystemSet *LSset, const int bound);

	double boundingVol() const;

	/**
	 * Get the size of this set, i.e,
	 * the number of linear systems
	 *
	 * @returns number of linear systems in the set
	 */
	inline int size() const{ return this->set.size(); }

	/**
	 * Get the set of linear systems
	 *
	 * @returns the current collection of linear systems
	 */
	inline const std::vector<LinearSystem*>& get_set() const { return set; }

	LinearSystem* at(int i);
	const LinearSystem* at(int i) const;

	bool isEmpty() const;
	void print() const;

	~LinearSystemSet();
};

#endif /* LINEARSYSTEMSET_H_ */
