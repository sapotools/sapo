#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE ode

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

#include "DifferentialSystem.h"
#include "Integrator.h"
#include "Evolver.h"

#define APPROX_ERR 1e-14

inline bool epsilon_equivalent(const Polytope<double>& P1,
                               const Polytope<double>& P2, const double epsilon)
{
    return is_true(((expand(P1, epsilon).includes(P2)) && 
                   (expand(P2, epsilon).includes(P1))));    
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_ode, T, test_types)
{
    using namespace SymbolicAlgebra;
    using namespace LinearAlgebra;

    Symbol<T> x("x"), y("y");
    ODE<T> ode({x,y},{-y,x}, "time");
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_runge_kutta, T, test_types)
{
    using namespace SymbolicAlgebra;
    using namespace LinearAlgebra;

    Symbol<T> x("x"), y("y");
    ODE<T> ode({x,y},{-y,x}, "time");

    RungeKutta4Integrator rk4;

    auto rk_system = rk4(ode,  Symbol<T>("time step"));
    auto time_step = rk_system.time_variable();

    std::vector<Expression<T>> ds{
        time_step*(-2*(0.5*time_step*x + y) + -y + -2*(0.5*time_step*(-0.5*y*time_step + x) + y) + -(time_step*(-0.5*(0.5*time_step*x + y)*time_step + x) + y))/6 + x,
        time_step*(2*(-0.5*y*time_step + x) + x + 2*(-0.5*(0.5*time_step*x + y)*time_step + x) + -(0.5*time_step*(-0.5*y*time_step + x) + y)*time_step + x)/6 + y
    };

    BOOST_REQUIRE(rk_system.dynamics()[0] - ds[0] == 0);
    BOOST_REQUIRE(rk_system.dynamics()[1] - ds[1] == 0);
}

BOOST_AUTO_TEST_CASE(test_transform_dynamic_bundle)
{
    using namespace SymbolicAlgebra;
    using namespace LinearAlgebra;

    Symbol<double> x("x"), y("y");
    ODE<double> ode({x,y},{-y,x}, "time");

    RungeKutta4Integrator rk4;

    auto rk_system = rk4(ode, 0.011990811705);

    Evolver<double> T(rk_system);

    Dense::Matrix<double> rA{
        {1,0},
        {0,1}
    };

    Bundle rSet(rA, {-0.05,0.05}, {9.95,10.05}, {}, {0,1});

    Bundle next = T(rSet);
    for (size_t i=0; i<523; ++i) {
        next = T(next);
    }
    BOOST_CHECK(epsilon_equivalent(rSet, next, 3e-7));
}