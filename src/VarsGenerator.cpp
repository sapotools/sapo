/**
 * @file VarsGenerator.cpp
 * Automatically generate variables for paralleltope generator functions.
 * For high dimensions declaring manually the variables can be tedious...
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
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

inline GiNaC::lst get_symbol_lst(const char *basename,
                                 const unsigned int number_of_symbols)
{
  return get_symbol_lst(std::string(basename), number_of_symbols);
}

/**
 * Constructor that instantiates the variable generator
 *
 * @param[in] dim dimension of the model/parallelotope
 */
VarsGenerator::VarsGenerator(const unsigned int dim):
    dim(dim), qs(get_symbol_lst("q", dim)), as(get_symbol_lst("q", dim)),
    bs(get_symbol_lst("q", dim)), ls(get_symbol_lst("q", dim))
{
  using namespace GiNaC;

  lst symbol_lst;

  for (unsigned int i = 0; i < dim; ++i) {
    std::ostringstream oss;
    oss << "u" << i;

    this->us.push_back(get_symbol_lst(oss.str(), 20));
  }
}

/**
 * Generate a box out of the variables
 *
 * @param[in] b box offsets
 * @returns generatred box
 */
LinearSystem VarsGenerator::genBox(const std::vector<double> &b) const
{
  using namespace std;

  unsigned int n = b.size() / 2;
  vector<double> Ai(n, 0);
  vector<vector<double>> A(2 * n, Ai);

  for (unsigned int i = 0; i < n; i++) {
    A[2 * i][i] = 1;
    A[2 * i + 1][i] = -1;
  }

  return LinearSystem(A, b);
}

VarsGenerator::~VarsGenerator()
{
  // TODO Auto-generated destructor stub
}
