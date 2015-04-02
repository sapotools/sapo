/*
 * Disjunction.cpp
 *
 *  Created on: May 20, 2014
 *      Author: Tommaso Dreossi
 */

#include "Disjunction.h"

Disjunction::Disjunction(STL * f1, STL * f2){

	this->f1=f1;
	this->f2=f2;
	type=DISJUNCTION;

};

void Disjunction::print(){
	cout<<"(";
	this->f1->print();
	cout<<") or (";
	this->f2->print();
	cout<<")";
}

Disjunction::~Disjunction() {}
