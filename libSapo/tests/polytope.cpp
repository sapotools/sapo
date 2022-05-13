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

    Vector<double> b = {1,2,3,3,2,1};
 
    Polytope p;

    p = Polytope(A,b);
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

    BOOST_CHECK(!p.is_empty());

    b[3]=-3;

    p = Polytope(A,b);
    BOOST_CHECK(p.is_empty());
}

BOOST_AUTO_TEST_CASE(test_contains_polytope)
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
        {0,1,1,0},
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

    BOOST_CHECK(p5.is_empty());

    BOOST_CHECK(p1.contains(p2));
    BOOST_CHECK(p1.contains(p1));
    BOOST_CHECK(p2.contains(p5));
    BOOST_CHECK(p3.contains(p2));
    BOOST_CHECK(p5.contains(p5));
    BOOST_CHECK(p1.contains(p6));
    BOOST_CHECK(!p3.contains(p1));
    BOOST_CHECK(!p1.contains(p3));
    BOOST_CHECK(!p2.contains(p3));
    BOOST_CHECK(!p1.contains(p4));
    BOOST_CHECK(!p4.contains(p1));
    BOOST_CHECK(!p5.contains(p1));
    BOOST_CHECK(!p2.contains(p1));
    BOOST_CHECK(!p2.contains(p3));
    BOOST_CHECK(!p6.contains(p1));
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
        {0,1,1,0},
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
    
    BOOST_CHECK(pa == pb);
    BOOST_CHECK(pa == p2);
    BOOST_CHECK(pa != p1);

    pa = intersect(p1, p3);
    pb = intersect(p3, p1);

    BOOST_CHECK(pa == pb);
    BOOST_CHECK(pa == p2);
    BOOST_CHECK(pa != p3);

    pa = intersect(p3, p2);
    pb = intersect(p2, p3);

    BOOST_CHECK(pa == pb);
    BOOST_CHECK(pa == p2);
    BOOST_CHECK(pa != p3);

    pa = intersect(p3, p2);
    pb = intersect(p2, p3);

    BOOST_CHECK(pa == pb);
    BOOST_CHECK(pa == p2);
    BOOST_CHECK(pa != p3);

    pa = intersect(p1, p4);
    pb = intersect(p4, p1);

    BOOST_CHECK(pa == pb);
    BOOST_CHECK(!p4.is_empty());
    BOOST_CHECK(!p1.is_empty());
    BOOST_CHECK(pa.is_empty());

    pa = intersect(p1, p5);
    pb = intersect(p5, p1);

    BOOST_CHECK(pa == pb);
    BOOST_CHECK(p5.is_empty());
    BOOST_CHECK(pa.is_empty());
    BOOST_CHECK(!p1.is_empty());

    pa = intersect(p1, p6);
    pb = intersect(p6, p1);

    BOOST_CHECK(pa == pb);
    BOOST_CHECK(pa == p6);
    BOOST_CHECK(pa != p1);
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
        {0,1,1,0},
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
    BOOST_CHECK(intersect(p1, p2) == pa);
    BOOST_CHECK(pa == pb);
   
    pa = p1;
    pb = p3;
    pa.intersect_with(p3);
    pb.intersect_with(p1);
    BOOST_CHECK(intersect(p1, p3) == pa);
    BOOST_CHECK(pa == pb);

    pa = p2;
    pb = p3;
    pa.intersect_with(p3);
    pb.intersect_with(p2);
    BOOST_CHECK(intersect(p2, p3) == pa);
    BOOST_CHECK(pa == pb);

    pa = p1;
    pb = p4;
    pa.intersect_with(p4);
    pb.intersect_with(p1);
    BOOST_CHECK(intersect(p1, p4) == pa);
    BOOST_CHECK(pa == pb);

    pa = p1;
    pb = p5;
    pa.intersect_with(p5);
    pb.intersect_with(p1);
    BOOST_CHECK(intersect(p1, p5) == pa);
    BOOST_CHECK(pa == pb);

    pa = p1;
    pb = p6;
    pa.intersect_with(p6);
    pb.intersect_with(p1);
    BOOST_CHECK(intersect(p1, p6) == pa);
    BOOST_CHECK(pa == pb);
}