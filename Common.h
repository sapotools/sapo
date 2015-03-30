/*
 * Common.hpp
 *
 *  Created on: Jul 25, 2014
 *      Author: Tommaso Dreossi
 */

#ifndef COMMON_HPP_
#define COMMON_HPP_

#include <stdio.h>
#include <iostream>
#include <vector>
#include <math.h>
#include <ginac/ginac.h>

using namespace std;
using namespace GiNaC;

enum formula_type {ATOM,CONJUNCTION,DISJUNCTION,UNTIL,ALWAYS,EVENTUALLY};

struct synthesizer_opt{
	bool largest_para_set;	// largest parameter set
};



#endif /* COMMON_HPP_ */
