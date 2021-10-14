/**
 * @file Until.cpp
 * Until STL formula
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#include "Until.h"

/**
 * Constructor that instantiates an Until formula (f1 U_[a,b] f2)
 *
 * @param[in] f1 left subformula
 * @param[in] a beginning of temporal interval
 * @param[in] b end of temporal interval
 * @param[in] f2 right subformula
 */
Until::Until(const std::shared_ptr<STL> f1, int a, int b, const std::shared_ptr<STL> f2): 
	STL(UNTIL), f1(f1), f2(f2), t_itvl(a,b) {}

/**
 * Print the formula
 */
void Until::print() const {
	cout<<"(";
	this->f1->print();
	cout<<") until_"<< this->t_itvl <<" (";
	this->f2->print();
	cout<<")";
}

Until::~Until() 
{}
