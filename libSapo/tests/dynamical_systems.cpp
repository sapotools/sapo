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

BOOST_AUTO_TEST_CASE(test_runge_kutta)
{
    using namespace SymbolicAlgebra;
    using namespace LinearAlgebra;

    Symbol<> x("x"), y("y");
    Symbol<> timestep{"timestep"};

    // variable vector
    /*
    DynamicalSystem<double> ds({
        {x, 2*x*(x*x-3)*(4*(x*x)-3)*(x*x+21*(y*y)-12)},
        {y, y*(35*(x*x*x*x*x*x)+105*(x*x*x*x)*(y*y*y*y*y)-315*(x*x*x*x)-63*(x*x)*(y*y*y*y)+378*(x*x)+27*(y*y*y*y*y*y)-189*(y*y*y*y)+378*(y*y)-216)}
    });
    */
    
    DynamicalSystem<double> ODE({
        {x, -y},
        {y, x}
    });
    
    DynamicalSystem<double> ds({
        {x, -y},
        {y, x}
    });

    auto rk_system = runge_kutta4(ODE.variables(), ODE.dynamics(), timestep);

    

}