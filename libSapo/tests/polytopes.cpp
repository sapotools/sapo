#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE polytope

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <sstream>

#include "Polytope.h"

BOOST_AUTO_TEST_CASE(test_polytope)
{
    using namespace LinearAlgebra;
    using namespace LinearAlgebra::Dense;

    Matrix<double> A = {
        {1,0,0},
        {0,1,0},
        {0,0,1},
        {-1,0,0},
        {0,-1,0},
        {0,0,-1}
    };

    Polytope p;

    p = Polytope(A,{1,2,3,3,2,1});
}

BOOST_AUTO_TEST_CASE(test_polytope_error)
{
    using namespace LinearAlgebra;
    using namespace LinearAlgebra::Dense;

    Matrix<double> A = {
        {1,0,0},
        {0,1,0}
    };

    BOOST_REQUIRE_THROW(Polytope(A, {1}), std::domain_error);
    BOOST_REQUIRE_THROW(Polytope(A, {1,2,3}), std::domain_error);
}

BOOST_AUTO_TEST_CASE(test_is_empty_polytope)
{
    using namespace LinearAlgebra;
    using namespace LinearAlgebra::Dense;

    Matrix<double> A = {
        {1,0,0},
        {0,1,0},
        {0,0,1},
        {-1,0,0},
        {0,-1,0},
        {0,0,-1}
    };

    Vector<double> b = {1,2,3,3,2,1};
 
    Polytope p(A,b);

    BOOST_CHECK(is_false(p.is_empty()));

    b[3]=-3;

    p = Polytope(A,b);
    BOOST_CHECK(is_true(p.is_empty()));
}

BOOST_AUTO_TEST_CASE(test_includes_polytope)
{
    using namespace LinearAlgebra;
    using namespace LinearAlgebra::Dense;

    Matrix<double> A = {
        {1,0,0},
        {0,1,0},
        {0,0,1},
        {-1,0,0},
        {0,-1,0},
        {0,0,-1}
    };

    Matrix<double> B = {
        {1,1,0},
        {1,-1,0},
        {0,1,1},
        {0,-1,1},
        {1,0,1},
        {1,0,-1},
        {-1,0,-1},
        {-1,0,1},
        {0,-1,-1},
        {0,1,-1},
        {-1,-1,0},
        {-1,1,0},
    };

    Polytope p1(A,{1,2,3,3,2,1}), p2(A,{1,2,3,0,0,0}), 
             p3(A,{3,4,5,0,0,0}), p4(A,{3,4,5,-2,-3,-4}), 
             p5(A,{1,2,3,-3,2,1}), 
             p6(B,{1,1,3,3,0,2,4,2,1,1,3,3});

    BOOST_CHECK(is_true(p5.is_empty()));

    BOOST_CHECK(is_true(p1.includes(p2)));
    BOOST_CHECK(is_true(p1.includes(p1)));
    BOOST_CHECK(is_true(p2.includes(p5)));
    BOOST_CHECK(is_true(p3.includes(p2)));
    BOOST_CHECK(is_true(p5.includes(p5)));
    BOOST_CHECK(is_true(p1.includes(p6)));
    BOOST_CHECK(is_false(p3.includes(p1)));
    BOOST_CHECK(is_false(p1.includes(p3)));
    BOOST_CHECK(is_false(p2.includes(p3)));
    BOOST_CHECK(is_false(p1.includes(p4)));
    BOOST_CHECK(is_false(p4.includes(p1)));
    BOOST_CHECK(is_false(p5.includes(p1)));
    BOOST_CHECK(is_false(p2.includes(p1)));
    BOOST_CHECK(is_false(p2.includes(p3)));
    BOOST_CHECK(is_false(p6.includes(p1)));
}

BOOST_AUTO_TEST_CASE(test_are_disjoint_polytope)
{
    using namespace LinearAlgebra;
    using namespace LinearAlgebra::Dense;

    Matrix<double> A = {
        {1,0,0},
        {0,1,0},
        {0,0,1},
        {-1,0,0},
        {0,-1,0},
        {0,0,-1}
    };

    Matrix<double> B = {
        {1,1},
        {1,-1},
    };

    Polytope p1(A,{1,2,3,3,2,1}), p2(A,{1,2,3,0,0,0}), 
             p3(A,{3,4,5,0,0,0}), p4(A,{3,4,5,-2,-3,-4}), 
             p5(A,{1,2,3,-3,2,1}), 
             p6(B,{1,1});

    BOOST_CHECK(is_false(are_disjoint(p1, p1)));
    BOOST_CHECK(is_false(are_disjoint(p1, p2)));
    BOOST_CHECK(is_false(are_disjoint(p1, p3)));
    BOOST_CHECK(is_true(are_disjoint(p1, p4)));
    BOOST_CHECK(is_true(are_disjoint(p1, p5)));
    BOOST_REQUIRE_THROW(are_disjoint(p1, p6), std::domain_error);

    BOOST_CHECK(is_false(are_disjoint(p2, p1)));
    BOOST_CHECK(is_false(are_disjoint(p2, p2)));
    BOOST_CHECK(is_false(are_disjoint(p2, p3)));
    BOOST_CHECK(is_true(are_disjoint(p2, p4)));
    BOOST_CHECK(is_true(are_disjoint(p2, p5)));
    BOOST_REQUIRE_THROW(are_disjoint(p2, p6), std::domain_error);

    BOOST_CHECK(is_false(are_disjoint(p3, p1)));
    BOOST_CHECK(is_false(are_disjoint(p3, p2)));
    BOOST_CHECK(is_false(are_disjoint(p3, p3)));
    BOOST_CHECK(is_false(are_disjoint(p3, p4)));
    BOOST_CHECK(is_true(are_disjoint(p3, p5)));
    BOOST_REQUIRE_THROW(are_disjoint(p3, p6), std::domain_error);

    BOOST_CHECK(is_true(are_disjoint(p4, p1)));
    BOOST_CHECK(is_true(are_disjoint(p4, p2)));
    BOOST_CHECK(is_false(are_disjoint(p4, p3)));
    BOOST_CHECK(is_false(are_disjoint(p4, p4)));
    BOOST_CHECK(is_true(are_disjoint(p4, p5)));
    BOOST_REQUIRE_THROW(are_disjoint(p4, p6), std::domain_error);

    BOOST_CHECK(is_true(are_disjoint(p5, p1)));
    BOOST_CHECK(is_true(are_disjoint(p5, p2)));
    BOOST_CHECK(is_true(are_disjoint(p5, p3)));
    BOOST_CHECK(is_true(are_disjoint(p5, p4)));
    BOOST_CHECK(is_true(are_disjoint(p5, p5)));
    BOOST_REQUIRE_THROW(are_disjoint(p5, p6), std::domain_error);

    BOOST_REQUIRE_THROW(are_disjoint(p6, p1), std::domain_error);
    BOOST_REQUIRE_THROW(are_disjoint(p6, p2), std::domain_error);
    BOOST_REQUIRE_THROW(are_disjoint(p6, p3), std::domain_error);
    BOOST_REQUIRE_THROW(are_disjoint(p6, p4), std::domain_error);
    BOOST_REQUIRE_THROW(are_disjoint(p6, p5), std::domain_error);
    BOOST_REQUIRE(is_false(are_disjoint(p6, p6)));
}

BOOST_AUTO_TEST_CASE(test_intersect_polytope)
{
    using namespace LinearAlgebra;
    using namespace LinearAlgebra::Dense;

    Matrix<double> A = {
        {1,0,0},
        {0,1,0},
        {0,0,1},
        {-1,0,0},
        {0,-1,0},
        {0,0,-1}
    };

    Matrix<double> B = {
        {1,1,0},
        {1,-1,0},
        {0,1,1},
        {0,-1,1},
        {1,0,1},
        {1,0,-1},
        {-1,0,-1},
        {-1,0,1},
        {0,-1,-1},
        {0,1,-1},
        {-1,-1,0},
        {-1,1,0},
    };

    Polytope p1(A,{1,2,3,3,2,1}), p2(A,{1,2,3,0,0,0}), 
             p3(A,{3,4,5,0,0,0}), p4(A,{3,4,5,-2,-3,-4}), 
             p5(A,{1,2,3,-3,2,1}), 
             p6(B,{1,1,3,3,0,2,4,2,1,1,3,3});

    Polytope pa = intersect(p1, p2),
             pb = intersect(p2, p1);
    
    BOOST_CHECK(is_true(pa == pb));
    BOOST_CHECK(is_true(pa == p2));
    BOOST_CHECK(is_true(pa != p1));

    pa = intersect(p1, p3);
    pb = intersect(p3, p1);

    BOOST_CHECK(is_true(pa == pb));
    BOOST_CHECK(is_true(pa == p2));
    BOOST_CHECK(is_true(pa != p3));

    pa = intersect(p3, p2);
    pb = intersect(p2, p3);

    BOOST_CHECK(is_true(pa == pb));
    BOOST_CHECK(is_true(pa == p2));
    BOOST_CHECK(is_true(pa != p3));

    pa = intersect(p3, p2);
    pb = intersect(p2, p3);

    BOOST_CHECK(is_true(pa == pb));
    BOOST_CHECK(is_true(pa == p2));
    BOOST_CHECK(is_true(pa != p3));

    pa = intersect(p1, p4);
    pb = intersect(p4, p1);

    BOOST_CHECK(is_true(pa == pb));
    BOOST_CHECK(is_false(p4.is_empty()));
    BOOST_CHECK(is_false(p1.is_empty()));
    BOOST_CHECK(is_true(pa.is_empty()));

    pa = intersect(p1, p5);
    pb = intersect(p5, p1);

    BOOST_CHECK(is_true(pa == pb));
    BOOST_CHECK(is_true(p5.is_empty()));
    BOOST_CHECK(is_true(pa.is_empty()));
    BOOST_CHECK(is_false(p1.is_empty()));

    pa = intersect(p1, p6);
    pb = intersect(p6, p1);

    BOOST_CHECK(is_true(pa == pb));
    BOOST_CHECK(is_true(pa == p6));
    BOOST_CHECK(is_true(pa != p1));

    pa = intersect(p2, p6);
    pb = intersect(p6, p2);

    BOOST_CHECK(is_true(pa == pb));
    BOOST_CHECK(is_true(p2.includes(pa)));
    BOOST_CHECK(is_true(p6.includes(pa)));
}

BOOST_AUTO_TEST_CASE(test_intersect_with_polytope)
{
    using namespace LinearAlgebra;
    using namespace LinearAlgebra::Dense;

    Matrix<double> A = {
        {1,0,0},
        {0,1,0},
        {0,0,1},
        {-1,0,0},
        {0,-1,0},
        {0,0,-1}
    };

    Matrix<double> B = {
        {1,1,0},
        {1,-1,0},
        {0,1,1},
        {0,-1,1},
        {1,0,1},
        {1,0,-1},
        {-1,0,-1},
        {-1,0,1},
        {0,-1,-1},
        {0,1,-1},
        {-1,-1,0},
        {-1,1,0},
    };
    Polytope p1(A,{1,2,3,3,2,1}), p2(A,{1,2,3,0,0,0}), 
             p3(A,{3,4,5,0,0,0}), p4(A,{3,4,5,-2,-3,-4}), 
             p5(A,{1,2,3,-3,2,1}),
             p6(B,{1,1,3,3,0,2,4,2,1,1,3,3});

    Polytope pa = p1, pb = p2;
    
    pa.intersect_with(p2);
    pb.intersect_with(p1);
    BOOST_CHECK(is_true(intersect(p1, p2) == pa));
    BOOST_CHECK(is_true(pa == pb));
   
    pa = p1;
    pb = p3;
    pa.intersect_with(p3);
    pb.intersect_with(p1);
    BOOST_CHECK(is_true(intersect(p1, p3) == pa));
    BOOST_CHECK(is_true(pa == pb));

    pa = p2;
    pb = p3;
    pa.intersect_with(p3);
    pb.intersect_with(p2);
    BOOST_CHECK(is_true(intersect(p2, p3) == pa));
    BOOST_CHECK(is_true(pa == pb));

    pa = p1;
    pb = p4;
    pa.intersect_with(p4);
    pb.intersect_with(p1);
    BOOST_CHECK(is_true(intersect(p1, p4) == pa));
    BOOST_CHECK(is_true(pa == pb));

    pa = p1;
    pb = p5;
    pa.intersect_with(p5);
    pb.intersect_with(p1);
    BOOST_CHECK(is_true(intersect(p1, p5) == pa));
    BOOST_CHECK(is_true(pa == pb));

    pa = p1;
    pb = p6;
    pa.intersect_with(p6);
    pb.intersect_with(p1);
    BOOST_CHECK(is_true(intersect(p1, p6) == pa));
    BOOST_CHECK(is_true(pa == pb));
}

BOOST_AUTO_TEST_CASE(test_over_approximate_union)
{
    using namespace LinearAlgebra;
    using namespace LinearAlgebra::Dense;

    Matrix<double> A = {
        {1,0},
        {0,1},
        {-1,0},
        {0,-1}
    };

    Matrix<double> B = {
        {1},
        {-1}
    };

    Polytope p1(A,{2,3,0,-1}), p2(A, {4,2,-1,0}), p3(A,{5,3,-3,-1.5}),
             p4(A,{12,3,-10,-1}), p5(A,{14,2,-11,0}), p6(A,{13.5,3.5,-11.5,-2.5}),
             p7(A,{1.5,2.5,-0.5,-1}), p8(B,{1,-1}), p9(A,{-3,-4,1,0});

    BOOST_CHECK(is_true(p1.includes(p7)));
    BOOST_CHECK(is_true(p9.is_empty()));

    Polytope res = over_approximate_union(p1, p2);
    
    BOOST_CHECK(is_true(res.includes(p1)));
    BOOST_CHECK(is_true(res.includes(p2)));
    BOOST_CHECK(is_false(res.includes(p3)));

    res = over_approximate_union(p1, p3);

    BOOST_CHECK(is_true(res.includes(p1)));
    BOOST_CHECK(is_false(res.includes(p2)));
    BOOST_CHECK(is_true(res.includes(p3)));

    res = over_approximate_union(p1, p9);

    BOOST_CHECK(is_true(res.includes(p1)));
    BOOST_CHECK(is_true(res.includes(p9)));

    BOOST_REQUIRE_THROW(over_approximate_union(p1, p8), std::domain_error);
}
