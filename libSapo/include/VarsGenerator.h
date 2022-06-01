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

#include "SymbolicAlgebra.h"

/**
 * @brief Generate a vector of symbols
 *
 * This method builds a vector of symbols whose
 * name have the form `basename_{vector_index}`.
 *
 * @param basename is the name prefix of the symbols
 * @param size is the size of the output vector
 * @return a vector of `size` symbols
 */
std::vector<SymbolicAlgebra::Symbol<>>
get_symbol_vector(const std::string &basename, const size_t size);

/**
 * @brief Generate a vector of symbols
 *
 * This method builds a vector of symbols whose
 * name have the form `basename_{vector_index}`.
 *
 * @param basename is the name prefix of the symbols
 * @param size is the size of the output vector
 * @return a vector of `size` symbols
 */
inline std::vector<SymbolicAlgebra::Symbol<>>
get_symbol_vector(const char *basename, const size_t size)
{
  return get_symbol_vector(std::string(basename), size);
}

#endif /* VARSGENERATOR_H_*/
