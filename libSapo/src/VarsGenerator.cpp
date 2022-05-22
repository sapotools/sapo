/**
 * @file VarsGenerator.cpp
 * Automatically generate variables.
 *
 * @author Alberto Casagrande <acasagrande@units.it>
 * @version 0.1
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
