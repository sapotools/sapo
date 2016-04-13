/*
 * Always.cpp
 *
 *  Created on: May 20, 2014
 *      Author: Tommaso Dreossi
 */

#include "Always.h"

Always::Always(int a, int b, STL * f){

	this->f=f;
	this->a = a;
	this->b = b;
	type=ALWAYS;

};

void Always::print(){
	cout<<"always_["<<this->a<<","<<this->b<<"] (";
	this->f->print();
	cout<<")";
}

Always::~Always() {}
