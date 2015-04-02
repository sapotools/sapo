/*
 * Polyhedron.cpp
 *
 *  Created on: Oct 30, 2014
 *      Author: Tommaso Dreossi
 */

#include "Polyhedron.h"

Polyhedron::Polyhedron() {
	// TODO Auto-generated constructor stub
}

// Get functions
lst Polyhedron::getQ(){ return this->vars[0]; }
lst Polyhedron::getAlpha(){ return this->vars[1]; }
lst Polyhedron::getBeta(){  return this->vars[2]; }
lst Polyhedron::getGeneratorFunction(){	return this->generator_function; }
int Polyhedron::getDim(){ return this->dim; }
vector< vector< double > > Polyhedron::getTemplate(){ return this->template_matrix; }

Polyhedron::~Polyhedron() {
	// TODO Auto-generated destructor stub
}

