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

struct sapo_opt{
	int trans;				// transformation (0: static, 1: dynamic)
	double alpha;			// decomposition weight
	int decomp;				// number of decompositions (0: none, >0: yes)
};

struct classcomp {
  bool operator() (vector<int> v1, vector<int> v2) const {
	  return lexicographical_compare(v1.begin(),v1.end(),v2.begin(),v2.end());
  }
};

struct poly_values{
	vector<double> base_vertex;
	vector<double> lenghts;
};

#endif /* COMMON_HPP_ */
