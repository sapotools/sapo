/**
 * @file Flowpipe.cpp
 * Represent and manipulate flowpipes of bundle
 * Used to represent the reachable set of a dynamical system
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#include <fstream>
#include <string>

#include "Flowpipe.h"

/**
 * Constructor that instantiates Flowpipe
 */
Flowpipe::Flowpipe(): std::vector<PolytopesUnion>() {}

/**
 * Return the i-th polytopes union
 *
 * @param[in] i index
 * @return i-th polytopes union
 */
const PolytopesUnion &Flowpipe::get(const unsigned int i) const
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
