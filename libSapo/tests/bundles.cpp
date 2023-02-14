#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE bundle

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include "Bundle.h"

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
 
    std::set<Vector<unsigned int>> T = {
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

    Bundle b1(A,{0,0,0,0,0},{5,5,5,3,7}, T),
           b2(B,{0,0,0,0,0},{3,5,5,7,5}, 
              {{2,1,4}}), // all the directions are static
           b3(B,{0,0,0,0,0},{3,5,5,7,5}, 
              {{2,1,4}},{0,3}), // 0 and 3 are adaptive directions
           b4(B,{0,0,0,0,0},{3,5,5,2,5}, 
              {{2,1,4}}), // all the directions are static
           b5(B,{0,0,0,0,0},{3,5,5,2,5}, 
              {{2,1,4}},{0,3}), // 0 and 3 are adaptive directions
           b6(B,{0,0,0,0,0},{3,5,5,7,5},
              {{2,1,4}},true), // remove unused directions
           b7(B,{0,0,0,0,0},{3,5,5,7,5},{{2,1,4},
              {0,1,4}},{0,1,4},true); // 0, 1, and 4 are adaptive directions
                                      // remove  unused directions


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
    BOOST_CHECK(b1==b3);
    BOOST_CHECK(b1!=b4);
    BOOST_CHECK(b1!=b5);
    BOOST_CHECK(b1!=b6);
    BOOST_CHECK(b4==b5);
    BOOST_CHECK(b1!=Polytope(Ap, {5,5,5,10,7,0,0,0,0,0}));
    BOOST_CHECK(b4.size()==B.size());
    BOOST_CHECK(b4.templates().size()==2);

    const auto num_b5_templates = ceil(float(B.size())/b5.dim());
    BOOST_CHECK(b5.size()==b4.size());
    BOOST_CHECK(b5.templates().size()==num_b5_templates);

    const auto num_b6_templates = 1;
    BOOST_CHECK(b6.size()==num_b6_templates*b6.dim());
    BOOST_CHECK(b6.templates().size()==num_b6_templates);

    const auto num_b7_templates = 2;
    BOOST_CHECK(b7.size()==num_b7_templates*b7.dim());
    BOOST_CHECK(b7.templates().size()==num_b7_templates);
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
 
    std::set<Vector<unsigned int>> T = {
        {0,1,2},
        {0,3,4}
    };

    BOOST_REQUIRE_THROW(Bundle b1(A,{0,0,0,0}, {5,5,5,3,7}, T);, std::domain_error);
    BOOST_REQUIRE_THROW(Bundle b1(A,{0,0,0,0,0}, {5,5,3,7}, T);, std::domain_error);
    BOOST_REQUIRE_THROW(Bundle b1(A,{0,0,0,0,0,0}, {5,5,5,3,7}, T);, std::domain_error);
    BOOST_REQUIRE_THROW(Bundle b1(A,{0,0,0,0,0}, {5,5,5,5,3,7}, T);, std::domain_error);
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
 
    std::set<Vector<unsigned int>> T = {
        {0,1,2},
        {0,3,4}
    };

    Bundle b1(A, {0,0,0,0,0}, {5,5,5,3,7}, T),
           b2(A, {0,0,0,0,0}, {-5,5,5,3,7}, T);

    BOOST_CHECK(b1.is_empty()==false);
    BOOST_CHECK(b2.is_empty()==true);
}


BOOST_AUTO_TEST_CASE(test_are_disjoint_polytope)
{
    using namespace LinearAlgebra;
    using namespace LinearAlgebra::Dense;

    Matrix<double> A = {
        {1,0,0},
        {0,1,0},
        {0,0,1}
    };

    Matrix<double> B = {
        {1,1},
        {1,-1},
    };

    Bundle b1(A,{-1,-2,-3},{3,2,1}), b2(A,{-1,-2,-3},{0,0,0}), 
           b3(A,{-3,-4,-5},{0,0,0}), b4(A,{-3,-4,-5},{-2,-3,-4}), 
           b5(A,{-1,-2,-3},{-3,2,1}), 
           b6(B,{0,0},{10,10});

    BOOST_CHECK(!are_disjoint(b1, b1));
    BOOST_CHECK(!are_disjoint(b1, b2));
    BOOST_CHECK(!are_disjoint(b1, b3));
    BOOST_CHECK(are_disjoint(b1, b4));
    BOOST_CHECK(are_disjoint(b1, b5));
    BOOST_REQUIRE_THROW(are_disjoint(b1, b6), std::domain_error);

    BOOST_CHECK(!are_disjoint(b2, b1));
    BOOST_CHECK(!are_disjoint(b2, b2));
    BOOST_CHECK(!are_disjoint(b2, b3));
    BOOST_CHECK(are_disjoint(b2, b4));
    BOOST_CHECK(are_disjoint(b2, b5));
    BOOST_REQUIRE_THROW(are_disjoint(b2, b6), std::domain_error);

    BOOST_CHECK(!are_disjoint(b3, b1));
    BOOST_CHECK(!are_disjoint(b3, b2));
    BOOST_CHECK(!are_disjoint(b3, b3));
    BOOST_CHECK(!are_disjoint(b3, b4));
    BOOST_CHECK(are_disjoint(b3, b5));
    BOOST_REQUIRE_THROW(are_disjoint(b3, b6), std::domain_error);

    BOOST_CHECK(are_disjoint(b4, b1));
    BOOST_CHECK(are_disjoint(b4, b2));
    BOOST_CHECK(!are_disjoint(b4, b3));
    BOOST_CHECK(!are_disjoint(b4, b4));
    BOOST_CHECK(are_disjoint(b4, b5));
    BOOST_REQUIRE_THROW(are_disjoint(b4, b6), std::domain_error);

    BOOST_CHECK(are_disjoint(b5, b1));
    BOOST_CHECK(are_disjoint(b5, b2));
    BOOST_CHECK(are_disjoint(b5, b3));
    BOOST_CHECK(are_disjoint(b5, b4));
    BOOST_CHECK(are_disjoint(b5, b5));
    BOOST_REQUIRE_THROW(are_disjoint(b5, b6), std::domain_error);

    BOOST_REQUIRE_THROW(are_disjoint(b6, b1), std::domain_error);
    BOOST_REQUIRE_THROW(are_disjoint(b6, b2), std::domain_error);
    BOOST_REQUIRE_THROW(are_disjoint(b6, b3), std::domain_error);
    BOOST_REQUIRE_THROW(are_disjoint(b6, b4), std::domain_error);
    BOOST_REQUIRE_THROW(are_disjoint(b6, b5), std::domain_error);
    BOOST_REQUIRE(!are_disjoint(b6, b6));
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
 
    std::set<Vector<unsigned int>> TA = {
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
 
    std::set<Vector<unsigned int>> TB = {
        {2,1,4},
        {2,0,3}
    };

    Bundle b1(A,{0,0,0,0,0},{5,5,5,3,7}, TA),
           b2(B,{0,0,0,0,0},{3,5,5,7,5}, TB),
           b3(A,{-1,-1,-1,-7,-7},{5,5,5,3,7}, TA),
           b4(A,{1,1,1,1,1},{5,5,5,3,7}, TA),
           b5(A,{0,0,0,0,0},{-5,5,5,3,7}, TA),
           b_err({{1}}, {0}, {1});

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
    BOOST_CHECK(ba.is_empty());
    BOOST_CHECK(b5.is_empty());
    BOOST_CHECK(!b1.is_empty());
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

    BOOST_REQUIRE_THROW(intersect(b1,b_err), std::domain_error);
    BOOST_REQUIRE_THROW(intersect(b_err,b1), std::domain_error);
}

BOOST_AUTO_TEST_CASE(test_is_subset_of_bundle)
{
    using namespace LinearAlgebra;

    std::vector<Vector<double>> A = {
        {1,0,0},
        {0,1,0},
        {0,0,1},
        {1,1,0},
        {0,1,1}
    };

    Bundle b1(A,{0,0,0,0,0},{5,5,5,3,7}),
           b2(A,{-1,0,-1,-3,-3},{6,7,8,9,10}),
           b3(A,{-1,-1,-1,-7,-7},{5,5,5,3,7}),
           b4(A,{1,1,1,1,1},{5,5,5,3,7}),
           b5(A,{0,0,0,0,0},{1,5,5,3,7}),
           b6(A,{0,0,0,0,0},{-5,5,5,3,7});

    for (const Bundle& bundle: {b1, b2, b3, b4, b5, b6}) {
        BOOST_CHECK(bundle.is_subset_of(bundle));
    }

    BOOST_CHECK(b1.is_subset_of(b2));
    BOOST_CHECK(!b2.is_subset_of(b1));

    BOOST_CHECK(b1.is_subset_of(b3));
    BOOST_CHECK(!b3.is_subset_of(b1));

    BOOST_CHECK(b4.is_subset_of(b1));
    BOOST_CHECK(!b1.is_subset_of(b4));

    BOOST_CHECK(!b1.is_subset_of(b5));
    BOOST_CHECK(b5.is_subset_of(b1));

    BOOST_CHECK(!b4.is_subset_of(b5));
    BOOST_CHECK(!b5.is_subset_of(b4));

    for (const Bundle& bundle: {b1, b2, b3, b4, b5, b6}) {
        BOOST_CHECK(b6.is_subset_of(bundle));

        if (!bundle.is_empty()) {
            BOOST_CHECK(!bundle.is_subset_of(b6));
        }
    }
}


BOOST_AUTO_TEST_CASE(test_includes_bundle)
{
    using namespace LinearAlgebra;

    std::vector<Vector<double>> A = {
        {1,0,0},
        {0,1,0},
        {0,0,1},
        {1,1,0},
        {0,1,1}
    };

    Bundle b1(A,{0,0,0,0,0},{5,5,5,3,7}),
           b2(A,{-1,0,-1,-3,-3},{6,7,8,9,10}),
           b3(A,{-1,-1,-1,-7,-7},{5,5,5,3,7}),
           b4(A,{1,1,1,1,1},{5,5,5,3,7}),
           b5(A,{0,0,0,0,0},{1,5,5,3,7}),
           b6(A,{0,0,0,0,0},{-5,5,5,3,7});

    for (const Bundle& bundle: {b1, b2, b3, b4, b5, b6}) {
        BOOST_CHECK(bundle.includes(bundle));
    }

    BOOST_CHECK(b2.includes(b1));
    BOOST_CHECK(!b1.includes(b2));

    BOOST_CHECK(b3.includes(b1));
    BOOST_CHECK(!b1.includes(b3));

    BOOST_CHECK(b1.includes(b4));
    BOOST_CHECK(!b4.includes(b1));

    BOOST_CHECK(!b5.includes(b1));
    BOOST_CHECK(b1.includes(b5));

    BOOST_CHECK(!b4.includes(b5));
    BOOST_CHECK(!b5.includes(b4));

    for (const Bundle& bundle: {b1, b2, b3, b4, b5, b6}) {
        BOOST_CHECK(bundle.includes(b6));

        if (!bundle.is_empty()) {
            BOOST_CHECK(!b6.includes(bundle));
        }
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
 
    std::set<Vector<unsigned int>> TA = {
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
    BOOST_CHECK(b.includes(ba));
    BOOST_CHECK(ba.satisfies(ls1));
    BOOST_CHECK(ba==Polytope(C, {5,5,5,3,7,5,7,0,0,0,0,0,0,0}));

    ba = b;
    ba.intersect_with(ls2);
    BOOST_CHECK(b.includes(ba));
    BOOST_CHECK(ba.satisfies(ls2));
    BOOST_CHECK(ba==b);

    ba = b;
    ba.intersect_with(ls3);
    BOOST_CHECK(b.includes(ba));
    BOOST_CHECK(ba.satisfies(ls3));
    BOOST_CHECK(ba==Polytope(C, {5,5,5,3,7,100,7,0,0,0,0,0,0,0}));

    ba = b;
    ba.intersect_with(ls4);
    BOOST_CHECK(b.includes(ba));
    BOOST_CHECK(ba.satisfies(ls4));
    BOOST_CHECK(ba==Polytope(C, {5,5,5,3,7,5,7,0,0,0,0,0,0,100}));

    std::vector<Vector<double>> D = {
        {1,0,0},
        {-1,0,0}
    };

    Polytope empty_poly(D, {0, -1});

    ba = b;
    ba.intersect_with(ls5);

    BOOST_CHECK(b.includes(ba));
    BOOST_CHECK(ba.satisfies(ls5));
    BOOST_CHECK(ba==empty_poly);

    ba = b;
    ba.intersect_with(ls6);
    BOOST_CHECK(b.includes(ba));
    BOOST_CHECK(ba.satisfies(ls6));
    BOOST_CHECK(ba==empty_poly);

    ba = b;
    ba.intersect_with(ls7);
    BOOST_CHECK(b.includes(ba));
    BOOST_CHECK(ba.satisfies(ls7));
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
 
    std::set<Vector<unsigned int>> TA = {
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

BOOST_AUTO_TEST_CASE(test_approximate_union_bundle)
{
    using namespace LinearAlgebra;

    std::vector<Vector<double>> A = {
        {1,0,0},
        {0,1,0},
        {0,0,1},
        {1,1,0},
        {0,1,1}
    };

    std::vector<Vector<double>> B = {
        {1,0,0},
        {0,1,0},
        {0,0,1}
    };
 
    std::vector<Vector<double>> C = {
        {0,1,0},
        {1,1,0},
        {0,1,1}
    };

    std::set<Vector<unsigned int>> TA = {
        {0,1,2},
        {0,3,4}
    };

    Bundle b1(A,{6,6,6,2,2},{12,12,12,15,15}, TA),
           b2(B,{0,0,0},{5,5,5}),
           b3(C,{0,0,0},{5,5,5}),
           b_empty(B,{0,0,6},{5,5,5}),
           b_err({{1}},{0},{1});

    std::vector<std::pair<Bundle, Bundle>> b_pairs{
        {b1,b2}, {b1,b3}, {b2,b3}, {b1, b_empty}
    };

    for (const auto& b_pair: b_pairs) {
        Bundle res = over_approximate_union(b_pair.first, b_pair.second);
        BOOST_CHECK(res.includes(b_pair.first));
        BOOST_CHECK(res.includes(b_pair.second));
        
        res = over_approximate_union(b_pair.second, b_pair.first);
        BOOST_CHECK(res.includes(b_pair.first));
        BOOST_CHECK(res.includes(b_pair.second));
    }

    BOOST_REQUIRE_THROW(intersect(b_err,b1), std::domain_error);
    BOOST_REQUIRE_THROW(intersect(b1,b_err), std::domain_error);
}