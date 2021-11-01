/**
 * @file VarsGenerator.cpp
 * Automatically generate variables.
 *
 * @author Alberto Casagrande <acasagrande@units.it>
 * @version 0.1
 */

#include "VarsGenerator.h"

#include <sstream>

GiNaC::lst get_symbol_lst(const std::string &basename,
                          const unsigned int number_of_symbols)
{
  using namespace GiNaC;

  GiNaC::lst symbol_lst;
  for (unsigned int i = 0; i < number_of_symbols; ++i) {
    std::ostringstream oss;
    oss << basename << i;

    symbol_lst.append(symbol(oss.str()));
  }

  return symbol_lst;
}
