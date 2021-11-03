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

#include <ginac/ginac.h>

GiNaC::lst get_symbol_lst(const std::string &basename,
                          const unsigned int number_of_symbols);

inline GiNaC::lst get_symbol_lst(const char *basename,
                                 const unsigned int number_of_symbols)
{
  return get_symbol_lst(std::string(basename), number_of_symbols);
}

#endif /* VARSGENERATOR_H_*/
