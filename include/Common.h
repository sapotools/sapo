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

// TODO: Move this enum to STL.h
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

struct poly_values { // numerical values for polytopes
  std::vector<double> base_vertex;
  std::vector<double> lenghts;
};

#endif /* COMMON_HPP_ */
