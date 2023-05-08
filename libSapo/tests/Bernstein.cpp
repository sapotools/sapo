#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Bernstein

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include "Bernstein.h"

#ifdef HAVE_GMP
#include <gmpxx.h>

typedef boost::mpl::list<double, mpq_class> test_types;
#else
typedef boost::mpl::list<double> test_types;
#endif

BOOST_AUTO_TEST_CASE_TEMPLATE(test_polynomial_Bernstein1, T, test_types)
{
  using namespace SymbolicAlgebra;

  Symbol<T> x("x"), y("y");

  auto coeffs = get_Bernstein_coefficients({x, y}, 2 * x*x);

  std::vector<T> results{0, 0, 2};

  BOOST_CHECK(coeffs.size() == results.size());
  for (size_t i = 0; i < results.size(); ++i) {
    BOOST_CHECK(coeffs[i] == results[i]);
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_polynomial_Bernstein2, T, test_types)
{
  using namespace SymbolicAlgebra;

  Symbol<T> x("x"), y("y");

  auto coeffs = get_Bernstein_coefficients({x, y}, 2 * x * y * y + x + 1);

  std::vector<T> results{1, 1, 1, 2, 2, 4};

  BOOST_CHECK(coeffs.size() == results.size());
  for (size_t i = 0; i < results.size(); ++i) {
    BOOST_CHECK(coeffs[i] == results[i]);
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_polynomial_Bernstein3, T, test_types)
{
  using namespace SymbolicAlgebra;

  Symbol<T> x("x"), y("y");

  auto coeffs = get_Bernstein_coefficients({x, y}, 2*(1+3*x)*(1+3*x));

  std::vector<T> results{2, 8, 32};

  BOOST_CHECK(coeffs.size() == results.size());
  for (size_t i = 0; i < results.size(); ++i) {
    BOOST_CHECK(coeffs[i] == results[i]);
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_polynomial_Bernstein4, T, test_types)
{
  using namespace SymbolicAlgebra;

  Symbol<T> x("x"), y("y");

  auto coeffs = get_Bernstein_coefficients({x, y}, 2*(3*x+1)*(2*y-2)*(2*y-2)+3*x+2);

  std::vector<T> results{10, 2, 2, 37, 5, 5};

  BOOST_CHECK(coeffs.size() == results.size());
  for (size_t i = 0; i < results.size(); ++i) {
    BOOST_CHECK(coeffs[i] == results[i]);
  }
}
 
BOOST_AUTO_TEST_CASE_TEMPLATE(test_rational_Bernstein, T, test_types)
{
  using namespace SymbolicAlgebra;

  Symbol<T> x("x"), y("y");

  std::vector<Symbol<T>> vars{x, y};

  auto f = 2 * pow(x,3) / (2 * pow(x,2) * pow(y,2) + x + 1);

  typename Expression<T>::replacement_type repl{{x, 1 + 3 * x}, {y, 1 + 3 * y}};

  f.replace(repl);

  auto coeff = get_Bernstein_coefficients(vars, f);

  std::vector<T> result{T(1)/2, T(1)/5, T(1)/17, T(8)/9, T(8)/27,
                        T(8)/99, T(8)/5, T(8)/17, T(8)/65, T(128)/37, 
                        T(128)/133, T(128)/517};

  for (size_t i = 0; i < result.size(); ++i) {
    BOOST_CHECK(coeff[i] == result[i]);
  }
  BOOST_CHECK(coeff.size() == result.size());
}
