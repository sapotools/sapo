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

	ex predicate;						// predicate
	vector<lst> predicateControlPts;	// control points associated to this atom
	int id;								// identifier

public:

	Atom(ex predicate, int id);

	ex getPredicate(){ return predicate; };
	vector<lst> getPredicateControlPts();
	void setPredicateControlPts(vector<lst>  predicateControlPts);
	void print(){ cout<<this->predicate<<" <= 0"; }
	int getID(){ return this->id; }

	virtual ~Atom();


};

#endif /* ATOM_H_ */
