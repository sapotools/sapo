/*
 * Atom.h
 *
 *  Created on: May 20, 2014
 *      Author: Tommaso Dreossi
 */

#ifndef ATOM_H_
#define ATOM_H_

#include "STL.h"

class Atom : public STL {
private:

	ex predicate;
	vector<lst> predicateControlPts;
	int id;

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
