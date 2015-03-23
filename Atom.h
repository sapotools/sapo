/*
 * Atom.h
 *
 *  Created on: May 20, 2014
 *      Author: nicola
 */

#ifndef ATOM_H_
#define ATOM_H_

#include "STL.h"

class Atom : public STL {
private:

	ex predicate;
	lst predicateControlPts;

public:

	Atom(ex predicate);

	ex getPredicate(){ return predicate; };
	lst getPredicateControlPts();
	void setPredicateControlPts(lst predicateControlPts);
	void print(){ cout<<this->predicate<<" <= 0"; }

	virtual ~Atom();


};

#endif /* ATOM_H_ */
