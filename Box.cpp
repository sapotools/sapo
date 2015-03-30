/*
 * Box.cpp
 *
 *  Created on: Nov 3, 2014
 *      Author: dreossi
 */

#include "Box.h"

Box::Box(vector< vector<double> > A, vector< double > b) {
	this->A = A;
	this->b = b;
}

// computes the volume of the box
// IMPORTANT: it's assumed that the constraint of parallel facets are next to each others
double Box::volume(){
	double vol = 1;
	for( int i=0; i<(signed)this->b.size(); i++ ){
		vol = vol * (this->b[2*i] + this->b[2*i + 1]);
	}
	return vol;
}


Box::~Box() {
	// TODO Auto-generated destructor stub
}

