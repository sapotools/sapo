/**
 * @file VarsGenerator.h
 * Automatically generate variables.
 *
 * @author Alberto Casagrande <acasagrande@units.it>
 * @version 0.1
 */

#ifndef VARSGENERATOR_H_
#define VARSGENERATOR_H_

#include <string>

#include "SymbolicAlgebra.h"

std::vector<SymbolicAlgebra::Symbol<>>
get_symbol_vector(const std::string &basename,
                  const unsigned int number_of_symbols);

inline std::vector<SymbolicAlgebra::Symbol<>>
get_symbol_vector(const char *basename, const unsigned int number_of_symbols)
{
  return get_symbol_vector(std::string(basename), number_of_symbols);
}

#endif /* VARSGENERATOR_H_*/
