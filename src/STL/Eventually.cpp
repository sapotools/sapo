/**
 * @file Eventually.cpp
 * Eventually STL formula
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#include "Eventually.h"

/**
 * Constructor that instantiates an Eventually formula (F_[a,b]f)
 *
 * @param[in] a beginning of temporal interval
 * @param[in] b end of temporal interval
 * @param[in] f subformula
 */
Eventually::Eventually(int a, int b, STL * f){
	this->f=f;
	this->a = a;
	this->b = b;
	type=EVENTUALLY;
};

/**
 * Print the formula
 */
void Eventually::print(){
	cout<<"Eventually_["<<this->a<<","<<this->b<<"] (";
	this->f->print();
	cout<<")";
}

Eventually::~Eventually() {}
