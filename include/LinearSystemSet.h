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

	void add(LinearSystem *ls);
	void add(const LinearSystem& ls);
	void add(LinearSystem&& ls);

	LinearSystemSet& simplify();

	LinearSystemSet* get_a_finer_covering() const;

	// operations on set
	LinearSystemSet* getIntersectionWith(const LinearSystemSet *LSset) const;

	// in-place set operations
	LinearSystemSet& unionWith(LinearSystemSet *LSset);
	LinearSystemSet& boundedUnionWith(LinearSystemSet *LSset, const unsigned int bound);

	double boundingVol() const;

	/**
	 * Get the size of this set, i.e,
	 * the number of linear systems
	 *
	 * @returns number of linear systems in the set
	 */
	inline unsigned int size() const{ return this->set.size(); }

	/**
	 * Get the number of variables
	 * 
	 * @returns number of columns of linear systems in the set
	 */
	unsigned int dim() const;

	/**
	 * Get the set of linear systems
	 *
	 * @returns the current collection of linear systems
	 */
	inline const std::vector<LinearSystem*>& get_set() const { return set; }

	LinearSystem* at(int i);
	const LinearSystem* at(int i) const;

	bool isEmpty() const;

	/**
	 * Print the set of linear systems
	 */
	void print() const;

	/**
	 * Print the linear system set in Matlab format (for plotregion script)
	 * 
	 * @param[in] os is the output stream
	 * @param[in] color color of the polytope to plot
	 */
	void plotRegion(std::ostream& os=std::cout, const char color=' ') const;

	~LinearSystemSet();
};

std::ostream& operator<<(std::ostream& out, const LinearSystemSet& ls);

#endif /* LINEARSYSTEMSET_H_ */
