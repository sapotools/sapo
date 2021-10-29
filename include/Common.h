/**
 * @file Common.h
 * Libraries, structures, and data-type shared by different classes.
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#ifndef COMMON_HPP_
#define COMMON_HPP_

#include <algorithm>
#include <ginac/ginac.h>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <vector>

using namespace std;
using namespace GiNaC;

enum formula_type {
  ATOM,
  CONJUNCTION,
  DISJUNCTION,
  UNTIL,
  ALWAYS,
  EVENTUALLY
};

struct synthesizer_opt {
  bool largest_para_set; // largest parameter set
};

struct sapo_opt {
  unsigned char trans; // transformation (0: static, 1: dynamic)
  double alpha;        // decomposition weight
  unsigned int decomp; // number of decompositions (0: none, >0: yes)
  string plot;         // the name of the file were to plot the reach set
  bool verbose;        // display info
};

struct poly_values { // numerical values for polytopes
  vector<double> base_vertex;
  vector<double> lenghts;
};

#endif /* COMMON_HPP_ */
