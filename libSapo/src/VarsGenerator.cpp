/**
 * @file VarsGenerator.cpp
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @author Alberto Casagrande <acasagrande@units.it>
 * @brief Automatically generate vectors of variables
 * @version 0.1
 * @date 2015-10-14
 * 
 * @copyright Copyright (c) 2015-2022
 */

#include "VarsGenerator.h"

#include <sstream>

std::vector<SymbolicAlgebra::Symbol<>>
get_symbol_vector(const std::string &basename, const size_t size)
{
  using namespace SymbolicAlgebra;

  std::vector<Symbol<>> symbol_vector;
  symbol_vector.reserve(size);
  for (unsigned int i = 0; i < size; ++i) {
    std::ostringstream oss;
    oss << basename << i;

    symbol_vector.push_back(Symbol<>(oss.str()));
  }

  return symbol_vector;
}
