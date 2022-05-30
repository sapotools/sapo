/**
 * @file Flowpipe.cpp
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @author Alberto Casagrande <acasagrande@units.it>
 * @brief Represent and manipulate reachability flowpipes
 * @version 0.2
 * @date 2022-05-04
 *
 * @copyright Copyright (c) 2016-2022
 */

#include <fstream>
#include <string>

#include "Flowpipe.h"
#include "Polytope.h"
#include "SetsUnion.h"

/**
 * Constructor that instantiates flowpipes
 */
Flowpipe::Flowpipe(): std::vector<SetsUnion<Polytope>>() {}

/**
 * Return the i-th polytopes union
 *
 * @param[in] i index
 * @return i-th polytopes union
 */
const SetsUnion<Polytope> &Flowpipe::get(const unsigned int i) const
{
  if (i >= this->size()) {
    std::domain_error("Flowpipe::get: i must be included in [0,"
                      "size()-1]");
  }

  return this->operator[](i);
}

unsigned int Flowpipe::dim() const
{
  if (this->empty()) {
    return 0;
  }

  return (this->front()).dim();
}
