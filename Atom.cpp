/*
 * Atom.cpp
 *
 *  Created on: May 20, 2014
 *      Author: nicola
 */

#include "Atom.h"


Atom::Atom(ex predicate) {

	type=ATOM;
	this->predicate = predicate;

}

lst Atom::getPredicateControlPts(){
	return this->predicateControlPts;
}

void Atom::setPredicateControlPts(lst predicateControlPts){

	this->predicateControlPts = predicateControlPts;

}


Atom::~Atom() {
	// TODO Auto-generated destructor stub
}
