#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE discrete_system

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#ifdef HAVE_GMP
#include <gmpxx.h>

template<>
struct std::is_arithmetic<mpq_class> : std::integral_constant<bool, true> {};

typedef boost::mpl::list<double, mpq_class> test_types;
#else
typedef boost::mpl::list<double> test_types;
#endif


#include "DiscreteSystem.h"

BOOST_AUTO_TEST_CASE_TEMPLATE(test_discrete_system, T, test_types)
{
    using namespace SymbolicAlgebra;
    using namespace LinearAlgebra;

    Symbol<T> s("s"), i("i"), r("r"), x("x");
    Symbol<T> alpha("alpha"), beta("beta");

    // variable vector
    std::vector<Symbol<T>> vars{s, i, r};

    // parameter vector
    std::vector<Symbol<T>> params{alpha, beta};

    // dynamic laws
    std::vector<Expression<T>> dyns{
        s-beta*s*i,          // s's law
        i+beta*s*i-alpha*i,  // i's law
        r+alpha*i            // r's law
    };

    DiscreteSystem<T>(vars, params, dyns);
    DiscreteSystem<T>({s, i, r}, {alpha, beta}, dyns);

    BOOST_REQUIRE_THROW(DiscreteSystem<T>({s, i}, {alpha, beta}, dyns), std::domain_error);
    BOOST_REQUIRE_THROW(DiscreteSystem<T>({s, i, r}, {alpha}, dyns), std::domain_error);
    BOOST_REQUIRE_THROW(DiscreteSystem<T>({s, i, i}, {alpha, beta}, dyns), std::domain_error);
    BOOST_REQUIRE_THROW(DiscreteSystem<T>({s, i, r}, {alpha, alpha}, dyns), std::domain_error);
    BOOST_REQUIRE_THROW(DiscreteSystem<T>({s, i, s}, {alpha, beta}, dyns), std::domain_error);
    BOOST_REQUIRE_THROW(DiscreteSystem<T>({s, i, r, x}, {alpha, beta}, dyns), std::domain_error);
    BOOST_REQUIRE_THROW(DiscreteSystem<T>({s, i, r}, {r, alpha, beta}, dyns), std::domain_error);

    std::map<Symbol<T>, Expression<T>> varDyn{
        {s, s-beta*s*i},
        {i, i+beta*s*i-alpha*i},
        {r, r+alpha*i}
    };

    DiscreteSystem<T>(varDyn, {alpha, beta});
    BOOST_REQUIRE_THROW(DiscreteSystem<T>(varDyn, {alpha}), std::domain_error);

    std::map<Symbol<T>, Expression<T>> varDyn2{
        {s, s-beta*s*i},
        {i, i+beta*s*i-alpha*i},
        {x, r+alpha*i}
    };
    BOOST_REQUIRE_THROW(DiscreteSystem<T>(varDyn2, {alpha, beta}), std::domain_error);
}
