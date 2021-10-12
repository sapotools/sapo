/**
 * @file Atom.h
 * Atomic STL formula
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#ifndef ATOM_H_
#define ATOM_H_

#include "STL.h"

class Atom : public STL {
private:
	
	static int num_of_atoms;
	
	ex predicate;						// predicate
	vector<lst> predicateControlPts;	// control points associated to this atom
	int id;								// identifier

public:

	Atom(const ex& predicate);

	inline const ex& getPredicate() const { return predicate; };

	/**
	 * Returns the control points associated with this atom
	 *
	 * @ returns vector of control points
	 */
	inline const vector<lst>& getPredicateControlPts() const { return this->predicateControlPts; }
	
	/**
	 * Associate a vector of control points to this atom
	 *
	 * @param[in] predicateControlPts vector of control points
	 */
	inline void setPredicateControlPts(const vector<lst>&  predicateControlPts) { this->predicateControlPts = predicateControlPts; }

	inline void print() const { cout<<this->predicate<<" <= 0"; }
	inline int getID() const { return this->id; }

	~Atom();
};

#endif /* ATOM_H_ */
