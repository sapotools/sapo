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
#include <algorithm>

using namespace std;
using namespace GiNaC;

enum formula_type {ATOM,CONJUNCTION,DISJUNCTION,UNTIL,ALWAYS,EVENTUALLY};

struct synthesizer_opt{
	bool largest_para_set;	// largest parameter set
};

struct classcomp {
  bool operator() (vector<int> v1, vector<int> v2) const {
	  return lexicographical_compare(v1.begin(),v1.end(),v2.begin(),v2.end());
  }
};



#endif /* COMMON_HPP_ */
