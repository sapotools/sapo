/*
 * Until.cpp
 *
 *  Created on: May 20, 2014
 *      Author: nicola
 */

#include "Until.h"

Until::Until(STL * f1, int a, int b, STL * f2){

	this->f1=f1;
	this->f2=f2;
	this->a = a;
	this->b = b;
	type=UNTIL;

};

void Until::print(){
	cout<<"(";
	this->f1->print();
	cout<<") until_["<<this->a<<","<<this->b<<"] (";
	this->f2->print();
	cout<<")";
}

Until::~Until() {}
