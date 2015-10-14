/*
 * Eventually.cpp
 *
 *  Created on: May 20, 2014
 *      Author: Tommaso Dreossi
 */

#include "Eventually.h"

Eventually::Eventually(int a, int b, STL * f){

	this->f=f;
	this->a = a;
	this->b = b;
	type=EVENTUALLY;

};

void Eventually::print(){
	cout<<"Eventually_["<<this->a<<","<<this->b<<"] (";
	this->f->print();
	cout<<")";
}

Eventually::~Eventually() {}
