#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE bundle

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <sstream>

#include "Bundle.h"

#define APPROX_ERR 1e-14

BOOST_AUTO_TEST_CASE(test_bundle)
{
    using namespace LinearAlgebra;

    std::vector<Vector<double>> A = {
        {1,0,0},
        {0,1,0},
        {0,0,1},
        {1,1,0},
        {0,1,1}
    };
 
    std::vector<Vector<unsigned int>> T = {
        {0,1,2},
        {0,3,4}
    };

    std::vector<Vector<double>> B = {
        {1,1,0},
        {0,1,0},
        {1,0,0},
        {0,1,1},
        {0,0,1}
    };

    Bundle b1(A,{0,0,0,0,0},{5,5,5,3,7},T),
           b2(B,{0,0,0,0,0},{3,5,5,7,5},{{2,1,4}}),
           b3(B,{0,0,0,0,0},{3,5,5,2,5},{{2,1,4}});


    BOOST_CHECK(b1.dim()==3);

    std::vector<Vector<double>> Ap = {
        {1,0,0},
        {0,1,0},
        {0,0,1},
        {1,1,0},
        {0,1,1},
        {-1,0,0},
        {0,-1,0},
        {0,0,-1},
        {-1,-1,0},
        {0,-1,-1}
    };

    BOOST_CHECK(b1==Polytope(Ap, {5,5,5,3,7,0,0,0,0,0}));
    BOOST_CHECK(b1==b2);
    BOOST_CHECK(b1!=b3);
    BOOST_CHECK(b1!=Polytope(Ap, {5,5,5,10,7,0,0,0,0,0}));
}

BOOST_AUTO_TEST_CASE(test_bundle_error)
{
    using namespace LinearAlgebra;

    std::vector<Vector<double>> A = {
        {1,0,0},
        {0,1,0},
        {0,0,1},
        {1,1,0},
        {0,1,1}
    };
 
    std::vector<Vector<unsigned int>> T = {
        {0,1,2},
        {0,3,4}
    };

    BOOST_REQUIRE_THROW(Bundle b1(A,{0,0,0,0}, {5,5,5,3,7},T);, std::domain_error);
    BOOST_REQUIRE_THROW(Bundle b1(A,{0,0,0,0,0}, {5,5,3,7},T);, std::domain_error);
    BOOST_REQUIRE_THROW(Bundle b1(A,{0,0,0,0,0,0}, {5,5,5,3,7},T);, std::domain_error);
    BOOST_REQUIRE_THROW(Bundle b1(A,{0,0,0,0,0}, {5,5,5,5,3,7},T);, std::domain_error);
    BOOST_REQUIRE_THROW(Bundle b1(A,{0,0,0,0,0}, {5,5,5,3,7}, {{7}});, std::domain_error);
    BOOST_REQUIRE_THROW(Bundle b1(A,{0,0,0,0,0}, {5,5,5,3,7}, {{7,0,0}}), std::domain_error);
    BOOST_REQUIRE_THROW(Bundle b1(A,{0,0,0,0,0}, {5,5,5,3,7}, {{0,0,0}}), std::domain_error);
}


BOOST_AUTO_TEST_CASE(test_is_empty_bundle)
{
    using namespace LinearAlgebra;

    std::vector<Vector<double>> A = {
        {1,0,0},
        {0,1,0},
        {0,0,1},
        {1,1,0},
        {0,1,1}
    };
 
    std::vector<Vector<unsigned int>> T = {
        {0,1,2},
        {0,3,4}
    };

    Bundle b1(A, {0,0,0,0,0}, {5,5,5,3,7}, T),
           b2(A, {0,0,0,0,0}, {-5,5,5,3,7}, T);

    BOOST_CHECK(((Polytope)b1).is_empty()==false);
    BOOST_CHECK(((Polytope)b2).is_empty()==true);
}

BOOST_AUTO_TEST_CASE(test_intersect_bundle)
{
    using namespace LinearAlgebra;

    std::vector<Vector<double>> A = {
        {1,0,0},
        {0,1,0},
        {0,0,1},
        {1,1,0},
        {0,1,1}
    };
 
    std::vector<Vector<unsigned int>> TA = {
        {0,1,2},
        {0,3,4}
    };


    std::vector<Vector<double>> B = {
        {1,1,0},
        {0,1,0},
        {1,0,0},
        {0,1,1},
        {0,0,1}
    };
 
    std::vector<Vector<unsigned int>> TB = {
        {2,1,4},
        {2,0,3}
    };

    Bundle b1(A,{0,0,0,0,0},{5,5,5,3,7}, TA),
           b2(B,{0,0,0,0,0},{3,5,5,7,5}, TB),
           b3(A,{-1,-1,-1,-7,-7},{5,5,5,3,7}, TA),
           b4(A,{1,1,1,1,1},{5,5,5,3,7}, TA),
           b5(A,{0,0,0,0,0},{-5,5,5,3,7}, TA);

    Bundle ba = intersect(b1,b1);
    BOOST_CHECK(ba==b1);
    
    ba = intersect(b1,b2);
    BOOST_CHECK(ba==b1);
    BOOST_CHECK(ba==b2);
    BOOST_CHECK(ba==intersect(b2,b1));

    ba = intersect(b1,b3);
    BOOST_CHECK(ba==b1);
    BOOST_CHECK(ba!=b3);
    BOOST_CHECK(ba==intersect(b3,b1));

    ba = intersect(b1,b4);
    BOOST_CHECK(ba==b4);
    BOOST_CHECK(ba!=b1);
    BOOST_CHECK(ba==intersect(b4,b1));

    ba = intersect(b1,b5);
    BOOST_CHECK(((Polytope)ba).is_empty());
    BOOST_CHECK(((Polytope)b5).is_empty());
    BOOST_CHECK(!((Polytope)b1).is_empty());
    BOOST_CHECK(ba==intersect(b5,b1));

    for (const auto& b: {b1,b2,b3,b4,b5}) {
        Bundle bi = intersect(b1,b);
        ba = b1;
        ba.intersect_with(b);
        BOOST_CHECK(ba==bi);

        ba = b;
        ba.intersect_with(b1);
        BOOST_CHECK(ba==bi);
    }
}

BOOST_AUTO_TEST_CASE(test_intersect_with_ls_bundle)
{
    using namespace LinearAlgebra;

    std::vector<Vector<double>> A = {
        {1,0,0},
        {0,1,0},
        {0,0,1},
        {1,1,0},
        {0,1,1}
    };
 
    std::vector<Vector<unsigned int>> TA = {
        {0,1,2},
        {0,3,4}
    };

    std::vector<Vector<double>> B = {{1,0,1},
                                     {-1,0,1}};

    Bundle b(A,{0,0,0,0,0},{5,5,5,3,7}, TA);

    LinearSystem ls1(B, {5, 0}),
                 ls2(B, {100, 100}),
                 ls3(B, {100, 0}),
                 ls4(B, {5, 100}),
                 ls5(B, {-1, 0}),
                 ls6(B, {5, -5}),
                 ls7(B, {-5, -5});

    std::vector<Vector<double>> C = {
        {1,0,0},
        {0,1,0},
        {0,0,1},
        {1,1,0},
        {0,1,1},
        {1,0,1},
        {1,0,-1},
        {-1,0,0},
        {0,-1,0},
        {0,0,-1},
        {-1,-1,0},
        {0,-1,-1},
        {-1,0,-1},
        {-1,0,1}
    };

    Bundle ba = b;
    ba.intersect_with(ls1);
    BOOST_CHECK(((Polytope)b).contains(ba));
    BOOST_CHECK(((Polytope)ba).satisfies(ls1));
    BOOST_CHECK(ba==Polytope(C, {5,5,5,3,7,5,7,0,0,0,0,0,0,0}));

    ba = b;
    ba.intersect_with(ls2);
    BOOST_CHECK(((Polytope)b).contains(ba));
    BOOST_CHECK(((Polytope)ba).satisfies(ls2));
    BOOST_CHECK(ba==b);

    ba = b;
    ba.intersect_with(ls3);
    BOOST_CHECK(((Polytope)b).contains(ba));
    BOOST_CHECK(((Polytope)ba).satisfies(ls3));
    BOOST_CHECK(ba==Polytope(C, {5,5,5,3,7,100,7,0,0,0,0,0,0,0}));

    ba = b;
    ba.intersect_with(ls4);
    BOOST_CHECK(((Polytope)b).contains(ba));
    BOOST_CHECK(((Polytope)ba).satisfies(ls4));
    BOOST_CHECK(ba==Polytope(C, {5,5,5,3,7,5,7,0,0,0,0,0,0,100}));


    std::vector<Vector<double>> D = {
        {1,0,0},
        {-1,0,0}
    };

    Polytope empty_poly(D, {0, -1});

    ba = b;
    ba.intersect_with(ls5);
    BOOST_CHECK(((Polytope)b).contains(ba));
    BOOST_CHECK(((Polytope)ba).satisfies(ls5));
    BOOST_CHECK(ba==empty_poly);

    ba = b;
    ba.intersect_with(ls6);
    BOOST_CHECK(((Polytope)b).contains(ba));
    BOOST_CHECK(((Polytope)ba).satisfies(ls6));
    BOOST_CHECK(ba==empty_poly);

    ba = b;
    ba.intersect_with(ls7);
    BOOST_CHECK(((Polytope)b).contains(ba));
    BOOST_CHECK(((Polytope)ba).satisfies(ls7));
    BOOST_CHECK(ba==empty_poly);
}

BOOST_AUTO_TEST_CASE(test_canonical_bundle)
{
    using namespace LinearAlgebra;

    std::vector<Vector<double>> A = {
        {1,0,0},
        {0,1,0},
        {0,0,1},
        {1,1,0},
        {0,1,1}
    };
 
    std::vector<Vector<unsigned int>> TA = {
        {0,1,2},
        {0,3,4}
    };

    Bundle b_orig(A,{0,0,0,-100,-100},{5,5,5,300,700}, TA);

    Bundle b(b_orig);

    Bundle bc = b.get_canonical();

    BOOST_CHECK(b_orig==b);
    BOOST_CHECK(bc==b);

    Vector<double> deltas = b_orig.upper_bounds()-b_orig.lower_bounds(),
                   deltas_c = bc.upper_bounds()-bc.lower_bounds();

    auto d_it = std::begin(deltas);
    auto dc_it = std::begin(deltas_c);
    for (; d_it != std::end(deltas); ++d_it, ++dc_it) {
        BOOST_CHECK(*dc_it <= *d_it);
    }

    BOOST_CHECK(norm_infinity(deltas_c) < norm_infinity(deltas));

    b.canonize();
    BOOST_CHECK(bc==b);

    Vector<double> deltas_b = b.upper_bounds()-b.lower_bounds();

    auto db_it = std::begin(deltas_b);

    dc_it = std::begin(deltas_c);
    for (; d_it != std::end(deltas); ++db_it, ++dc_it) {
        BOOST_CHECK(*dc_it == *db_it);
    }
}
