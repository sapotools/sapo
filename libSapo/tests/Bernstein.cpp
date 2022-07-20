#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Bernstein

#include <boost/test/unit_test.hpp>
#include <gmpxx.h>

#include "Bundle.h"

#include "SymbolicAlgebra.h"
#include "Bernstein.h"

BOOST_AUTO_TEST_CASE(test_polynomial_Bernstein_double)
{
  using namespace SymbolicAlgebra;

  Symbol<> x("x"), y("y");

  std::vector<Symbol<>> vars{x, y};

  auto f1 = 2 * x*x;
  auto f2 = (2 * x * y * y + x + 1);

  Expression<>::replacement_type repl{{x, 1 + 3 * x}, {y, -2 + 2 * y}};

  f1.replace(repl);
  f2.replace(repl);

  auto coeff1 = get_Bernstein_coefficients(vars, f1);
  auto coeff2 = get_Bernstein_coefficients(vars, f2);

  std::vector<double> result1{2, 8, 32};
  std::vector<double> result2{10, 2, 2, 37, 5, 5};

  for (size_t i = 0; i < result1.size(); ++i) {
    BOOST_CHECK(coeff1[i] == result1[i]);
  }
  BOOST_CHECK(coeff1.size() == result1.size());

  for (size_t i = 0; i < result2.size(); ++i) {
    BOOST_CHECK(coeff2[i] == result2[i]);
  }
  BOOST_CHECK(coeff2.size() == result2.size());
}

BOOST_AUTO_TEST_CASE(test_polynomial_Bernstein_Q)
{
  using namespace SymbolicAlgebra;

  Symbol<mpq_class> x("x"), y("y");

  std::vector<Symbol<mpq_class>> vars{x, y};

  auto f1 = 2 * x*x;
  auto f2 = (2 * x * y * y + x + 1);

  Expression<mpq_class>::replacement_type repl{{x, 1 + 3 * x}, {y, -2 + 2 * y}};

  f1.replace(repl);
  f2.replace(repl);

  auto coeff1 = get_Bernstein_coefficients(vars, f1);
  auto coeff2 = get_Bernstein_coefficients(vars, f2);

  std::vector<mpq_class> result1{2, 8, 32};
  std::vector<mpq_class> result2{10, 2, 2, 37, 5, 5};

  for (size_t i = 0; i < result1.size(); ++i) {
    BOOST_CHECK(coeff1[i] == result1[i]);
  }
  BOOST_CHECK(coeff1.size() == result1.size());

  for (size_t i = 0; i < result2.size(); ++i) {
    BOOST_CHECK(coeff2[i] == result2[i]);
  }
  BOOST_CHECK(coeff2.size() == result2.size());
}

BOOST_AUTO_TEST_CASE(test_rational_Bernstein_double)
{
  using namespace SymbolicAlgebra;

  Symbol<> x("x"), y("y");

  std::vector<Symbol<>> vars{x, y};

  auto f = 2 * pow(x,3) / (2 * pow(x,2) * pow(y,2) + x + 1);

  Expression<>::replacement_type repl{{x, 1 + 3 * x}, {y, 1 + 3 * y}};

  f.replace(repl);

  auto coeff = get_Bernstein_coefficients(vars, f);

  std::vector<double> result{0.5,
                             0.2,
                             0.058823529411764705066012481893267,
                             0.8888888888888888,
                             0.2962962962962963,
                             0.08080808080808081,
                             1.6,
                             0.47058823529411764052809985514614,
                             0.12307692307692308375521861307789,
                             3.459459459459459651498036691919,
                             0.9624060150375939,
                             0.24758220502901354120872667863296};

  for (size_t i = 0; i < result.size(); ++i) {
    BOOST_CHECK(coeff[i] == result[i]);
  }
  BOOST_CHECK(coeff.size() == result.size());
}

BOOST_AUTO_TEST_CASE(test_rational_Bernstein_Q)
{
  using namespace SymbolicAlgebra;

  Symbol<mpq_class> x("x"), y("y");

  std::vector<Symbol<mpq_class>> vars{x, y};

  auto f = 2 * pow(x,3) / (2 * pow(x,2) * pow(y,2) + x + 1);

  Expression<mpq_class>::replacement_type repl{{x, 1 + 3 * x}, {y, 1 + 3 * y}};

  f.replace(repl);

  auto coeff = get_Bernstein_coefficients(vars, f);

  std::vector<mpq_class> result{{1,2}, {1,5}, {1,17}, {8,9}, {8,27},
                                {8,99}, {8,5}, {8,17}, {8,65},
                                {128,37}, {128,133}, {128,517}};

  for (size_t i = 0; i < result.size(); ++i) {
    BOOST_CHECK(coeff[i] == result[i]);
  }
  BOOST_CHECK(coeff.size() == result.size());
}