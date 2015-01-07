/*
 * Conjunction.cpp
 *
 *  Created on: May 20, 2014
 *      Author: nicola
 */

#include "Conjunction.h"

Conjunction::Conjunction(STL * f1, STL * f2){

	this->f1=f1;
	this->f2=f2;
	type=CONJUNCTION;

};

void Conjunction::print(){
	cout<<"(";
	this->f1->print();
	cout<<") and (";
	this->f2->print();
	cout<<")";
}

Conjunction::~Conjunction() {}
