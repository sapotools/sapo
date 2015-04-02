/*
 * Atom.cpp
 *
 *  Created on: May 20, 2014
 *      Author: Tommaso Dreossi
 */

#include "Atom.h"


Atom::Atom(ex predicate) {

	type=ATOM;
	this->predicate = predicate;

}

vector<lst> Atom::getPredicateControlPts(){
	return this->predicateControlPts;
}

void Atom::setPredicateControlPts(vector<lst> predicateControlPts){
	this->predicateControlPts = predicateControlPts;
}


Atom::~Atom() {
	// TODO Auto-generated destructor stub
}
