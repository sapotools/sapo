/**
 * @file VarsGenerator.h
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @author Alberto Casagrande <acasagrande@units.it>
 * @brief Automatically generate vectors of variables
 * @version 0.1
 * @date 2015-10-14
 *
 * @copyright Copyright (c) 2015-2022
 */

#ifndef VARSGENERATOR_H_
#define VARSGENERATOR_H_

#include <string>
#include <sstream>

#include "SymbolicAlgebra.h"

/**
 * @brief Generate a vector of symbols
 *
 * This method builds a vector of symbols whose
 * name have the form `basename_{vector_index}`.
 *
 * @tparam T is the domain of the symbol
 * @param basename is the name prefix of the symbols
 * @param size is the size of the output vector
 * @return a vector of `size` symbols
 */
template<typename T>
std::vector<SymbolicAlgebra::Symbol<T>>
get_symbol_vector(const std::string &basename, const size_t size)
{
  using namespace SymbolicAlgebra;

  std::vector<Symbol<T>> symbol_vector;
  symbol_vector.reserve(size);
  for (unsigned int i = 0; i < size; ++i) {
    std::ostringstream oss;
    oss << basename << i;

    symbol_vector.push_back(Symbol<T>(oss.str()));
  }

  return symbol_vector;
}

/**
 * @brief Generate a vector of symbols
 *
 * This method builds a vector of symbols whose
 * name have the form `basename_{vector_index}`.
 *
 * @tparam T is the domain of the symbols
 * @param basename is the name prefix of the symbols
 * @param size is the size of the output vector
 * @return a vector of `size` symbols
 */
template<typename T>
inline std::vector<SymbolicAlgebra::Symbol<T>>
get_symbol_vector(const char *basename, const size_t size)
{
  return get_symbol_vector<T>(std::string(basename), size);
}

#endif /* VARSGENERATOR_H_*/
