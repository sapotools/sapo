/**
 * @file Atom.cpp
 * Atomic STL formula
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#include "Atom.h"

 /**
  * Constructor that instantiates a atomic predicate
  *
  * @param[in] predicate a symbolic expression
  * @param[in] id an identifier for the atomic formula
  */
Atom::Atom(ex predicate, int id) {
	type=ATOM;
	this->predicate = predicate;
	this->id = id;
}

/**
 * Returns the control points associated with this atom
 *
 * @ returns vector of control points
 */
vector<lst> Atom::getPredicateControlPts(){
	return this->predicateControlPts;
}

/**
 * Associate a vector of control points to this atom
 *
 * @param[in] predicateControlPts vector of control points
 */
void Atom::setPredicateControlPts(vector<lst> predicateControlPts){
	this->predicateControlPts = predicateControlPts;
}


Atom::~Atom() {
	// TODO Auto-generated destructor stub
}
