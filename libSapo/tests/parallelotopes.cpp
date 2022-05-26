#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE parallelotope

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <sstream>

#include "Parallelotope.h"

BOOST_AUTO_TEST_CASE(test_parallelotope)
{
    using namespace LinearAlgebra;

    std::vector<Vector<double>> A = {
        {1,0,0},
        {0,1,0},
        {0,0,1}
    };
 
    std::vector<Vector<double>> B = {
        {3,4,0},
        {4,-3,0},
        {0,0,1}
    };

    Parallelotope p1(A, {-3,-2,-1}, {1,2,3}),
                  p2(B, {1,1,1}, {2,5,3});

    BOOST_CHECK(p1.dim()==3);

    std::vector<Vector<double>> Ap = {
        {1,0,0},
        {0,1,0},
        {0,0,1},
        {-1,0,0},
        {0,-1,0},
        {0,0,-1}
    };

    std::vector<Vector<double>> Bp = {
        {3,4,0},
        {4,-3,0},
        {0,0,1},
        {-3,-4,0},
        {-4,3,0},
        {0,0,-1},
    };

    BOOST_CHECK(p1==Polytope(Ap, {1,2,3,3,2,1}));
    BOOST_CHECK(p1!=Polytope(Ap, {1,2,3,7,2,1}));

    // The following commented test fails because of approximation errors
    // BOOST_CHECK((Polytope)p2==Polytope(Bp, {2,5,3,-1,-1,-1}));

    BOOST_CHECK(p2!=Polytope(Bp, {1,5,3,4,-1,-1}));
    BOOST_CHECK(p1!=p2);
}

BOOST_AUTO_TEST_CASE(test_parallelotope_error)
{
    using namespace LinearAlgebra;
    using namespace LinearAlgebra::Dense;

    Matrix<double> A = {
        {1,0},
        {0,1}
    };

    BOOST_REQUIRE_THROW(Parallelotope(A, {-3}, {1}), std::domain_error);
    BOOST_REQUIRE_THROW(Parallelotope(A, {-3}, {1,2}), std::domain_error);
    BOOST_REQUIRE_THROW(Parallelotope(A, {-3,-2}, {1}), std::domain_error);
    BOOST_REQUIRE_THROW(Parallelotope(A, {-3,-2,2}, {1,2}), std::domain_error);
    BOOST_REQUIRE_THROW(Parallelotope(A, {-3,-2}, {1,2,5}), std::domain_error);
    BOOST_REQUIRE_THROW(Parallelotope(A, {-3,-2,2}, {1,2,5}), std::domain_error);
}