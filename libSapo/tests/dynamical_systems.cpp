#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE dynamics

#include <boost/test/unit_test.hpp>

#include <sstream>

#include "DynamicalSystem.h"

BOOST_AUTO_TEST_CASE(test_dynamical_system)
{
    using namespace SymbolicAlgebra;
    using namespace LinearAlgebra;

    Symbol<> s("s"), i("i"), r("r"), x("x");
    Symbol<> alpha("alpha"), beta("beta");

    // variable vector
    std::vector<Symbol<>> vars{s, i, r};

    // parameter vector
    std::vector<Symbol<>> params{alpha, beta};

    // dynamic laws
    std::vector<Expression<>> dyns{
        s-beta*s*i,          // s's law
        i+beta*s*i-alpha*i,  // i's law
        r+alpha*i            // r's law
    };

    DynamicalSystem<double>(vars, params, dyns);
    DynamicalSystem<double>({s, i, r}, {alpha, beta}, dyns);

    BOOST_REQUIRE_THROW(DynamicalSystem<double>({s, i}, {alpha, beta}, dyns), std::domain_error);
    BOOST_REQUIRE_THROW(DynamicalSystem<double>({s, i, r}, {alpha}, dyns), std::domain_error);
    BOOST_REQUIRE_THROW(DynamicalSystem<double>({s, i, i}, {alpha, beta}, dyns), std::domain_error);
    BOOST_REQUIRE_THROW(DynamicalSystem<double>({s, i, r}, {alpha, alpha}, dyns), std::domain_error);
    BOOST_REQUIRE_THROW(DynamicalSystem<double>({s, i, s}, {alpha, beta}, dyns), std::domain_error);
    BOOST_REQUIRE_THROW(DynamicalSystem<double>({s, i, r, x}, {alpha, beta}, dyns), std::domain_error);
    BOOST_REQUIRE_THROW(DynamicalSystem<double>({s, i, r}, {r, alpha, beta}, dyns), std::domain_error);

    std::map<Symbol<>, Expression<>> varDyn{
        {s, s-beta*s*i},
        {i, i+beta*s*i-alpha*i},
        {r, r+alpha*i}
    };

    DynamicalSystem<double>(varDyn, {alpha, beta});
    BOOST_REQUIRE_THROW(DynamicalSystem<double>(varDyn, {alpha}), std::domain_error);
}
