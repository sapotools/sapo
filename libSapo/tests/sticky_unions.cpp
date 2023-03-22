#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE sticky_union

#include <boost/test/unit_test.hpp>

#include <sstream>

#include "Polytope.h"
#include "Bundle.h"
#include "StickyUnion.h"
#include "SetsUnion.h"

BOOST_AUTO_TEST_CASE(test_sticky_union)
{
    using namespace LinearAlgebra;
    using namespace LinearAlgebra::Dense;

    Matrix<double> A = {
        {1,0},
        {0,1}
    };

    Matrix<double> B = {
        {1}
    };

    Bundle b1(A,{3,1.5}, {5,3}), b2(A, {2,0}, {4,2}), b3(A,{0,1}, {2,3}),
           b4(A,{10,1}, {12,3}), b5(A,{11,0},{14,2}), b6(A,{13,2.5},{15,3.5}),
           b7(A,{0.5,1},{1.5,2.5}), b8(A,{4.5,1},{6,2}), b9(A,{-1,0},{-3,-4}),
           b10(B,{0},{1});

    Bundle c1(A,{0,1.5},{5.5,3}), c2(A,{10.5,1.5},{13,2});

    BOOST_CHECK(is_true(b9.is_empty()));

    BOOST_CHECK(is_true(b3.includes(b7)));

    BOOST_CHECK(is_false(are_disjoint(b1,b2)));
    BOOST_CHECK(is_true(are_disjoint(b1,b3)));
    BOOST_CHECK(is_true(are_disjoint(b1,b4)));
    BOOST_CHECK(is_true(are_disjoint(b1,b5)));
    BOOST_CHECK(is_true(are_disjoint(b1,b6)));

    BOOST_CHECK(is_false(are_disjoint(b2,b3)));
    BOOST_CHECK(is_true(are_disjoint(b2,b4)));
    BOOST_CHECK(is_true(are_disjoint(b2,b5)));
    BOOST_CHECK(is_true(are_disjoint(b2,b6)));

    BOOST_CHECK(is_true(are_disjoint(b3,b4)));
    BOOST_CHECK(is_true(are_disjoint(b3,b5)));
    BOOST_CHECK(is_true(are_disjoint(b3,b6)));

    BOOST_CHECK(is_false(are_disjoint(b4,b5)));
    BOOST_CHECK(is_true(are_disjoint(b4,b6)));

    BOOST_CHECK(is_true(are_disjoint(b5,b6)));

    std::list<Bundle> b_list{b1,b2,b3,b4,b5,b6,b7,b8,b9};

    StickyUnion<Bundle> Su(b_list);

    BOOST_CHECK(is_false(Su.any_includes(c1)));
    BOOST_CHECK(is_false(Su.any_includes(c2)));

    BOOST_CHECK(Su.size()==7);
    BOOST_CHECK(Su.number_of_classes()==3);

    for (auto l_it = std::begin(b_list); l_it != std::end(b_list); ++l_it) {
        BOOST_CHECK(is_true(Su.includes(*l_it)));
    }

    BOOST_CHECK(is_true(Su.includes(c1)));
    BOOST_CHECK(is_true(Su.includes(c2)));

    BOOST_REQUIRE_THROW(StickyUnion<Bundle>({b1, b10}), std::domain_error);
}

BOOST_AUTO_TEST_CASE(test_any_includes)
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

    Polytope p1(A,{3,4,5,-2,-3,-4}), p2(A,{1,2,3,0,0,0}),
             p3(A,{1,2,3,-3,2,1}), p4(A,{1,2,3,3,2,1}),
             p5(A,{0.5,1.5,2.5,-0.5,-0.5,-0.5}); 

    StickyUnion<Polytope> Su1, Su2(p1), Su3({p1, p2}), Su4(p3);

    BOOST_CHECK(is_false(p1.is_empty()));
    BOOST_CHECK(is_false(p2.is_empty()));
    BOOST_CHECK(is_true(p3.is_empty()));
    BOOST_CHECK(is_false(p4.is_empty()));
    BOOST_CHECK(is_false(p5.is_empty()));

    BOOST_CHECK(is_false(Su1.any_includes(p4)));
    BOOST_CHECK(is_false(Su2.any_includes(p4)));
    BOOST_CHECK(is_false(Su3.any_includes(p4)));
    BOOST_CHECK(is_false(Su4.any_includes(p4)));

    BOOST_CHECK(is_true(Su1.any_includes(p3)));
    BOOST_CHECK(is_true(Su2.any_includes(p3)));
    BOOST_CHECK(is_true(Su3.any_includes(p3)));
    BOOST_CHECK(is_true(Su4.any_includes(p3)));

    BOOST_CHECK(is_false(Su1.any_includes(p5)));
    BOOST_CHECK(is_false(Su2.any_includes(p5)));
    BOOST_CHECK(is_true(Su3.any_includes(p5)));
    BOOST_CHECK(is_false(Su4.any_includes(p5)));
}


BOOST_AUTO_TEST_CASE(test_copy_sticky_union)
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

    Matrix<double> C = {
        {1,0},
        {0,1},
        {-1,0},
        {0,-1}
    };

    Polytope p1(A,{1,2,3,0,0,0}), p2(A,{3,4,5,-2,-3,-4}),
             p3(B,{1,1,3,3,0,2,4,2,1,1,3,3}),
             p4(C,{1,1,0,0});

    BOOST_CHECK(is_false(p1.is_empty()));
    BOOST_CHECK(is_false(p2.is_empty()));
    BOOST_CHECK(is_false(p3.is_empty()));
    BOOST_CHECK(is_false(p4.is_empty()));

    StickyUnion<Polytope> Su1({p1,p2,p3}), Su2(p4), Su3, res;
    BOOST_CHECK(Su1.size()==3);

    res = Su3;

    BOOST_CHECK(is_true(res.is_empty()));

    for (const auto &Su:{Su1, Su2, Su3}) {
        res = Su;

        BOOST_CHECK(is_true(res.is_empty()==Su.is_empty()));
        BOOST_CHECK(res.size()==Su.size());
        BOOST_CHECK(res.number_of_classes()==Su.number_of_classes());

        SetsUnion<Polytope> Pu = Su;
        SetsUnion<Polytope> Pu_res = res;
        BOOST_CHECK(is_true(Pu_res.is_empty()==Pu.is_empty()));
        BOOST_CHECK(Pu_res.size()==Pu.size());

        auto res_it = std::end(Pu_res);
        auto p_it = std::end(Pu);

        for (; res_it != std::end(Pu_res); ++res_it, ++p_it) {
            BOOST_CHECK(is_true(*res_it == *p_it));
        }
    }
}

BOOST_AUTO_TEST_CASE(test_intersect_sets_union)
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

    Matrix<double> C = {
        {1,0},
        {0,1},
        {-1,0},
        {0,-1}
    };

    Polytope p1(A,{3,4,6,3,2,1}), p2(A,{1,2,3,0,0,0}), 
             p3(A,{4,4,5,0,0,0}), p4(A,{3,4,6,-2,-3,-4}), 
             p5(A,{1,2,3,-3,2,1}), 
             p6(B,{1,1,3,3,0,2,4,2,1,1,3,3}),
             p7(C,{1,1,0,0}), p8(A,{0,1,3,3,2,1}), 
             p9(A,{4,4,5,-1,-1,-1});

    BOOST_CHECK(is_false(p1.is_empty()));
    BOOST_CHECK(is_false(p2.is_empty()));
    BOOST_CHECK(is_false(p3.is_empty()));
    BOOST_CHECK(is_false(p4.is_empty()));
    BOOST_CHECK(is_true(p5.is_empty()));
    BOOST_CHECK(is_false(p6.is_empty()));
    BOOST_CHECK(is_false(p7.is_empty()));
    BOOST_CHECK(is_false(p8.is_empty()));
    BOOST_CHECK(is_false(p9.is_empty()));

    BOOST_CHECK(is_true(p1.includes(p2)));
    BOOST_CHECK(is_true(p1.includes(p4)));
    BOOST_CHECK(is_true(p1.includes(p6)));

    BOOST_CHECK(is_false(p2.includes(p1)));
    BOOST_CHECK(is_false(p4.includes(p1)));
    BOOST_CHECK(is_false(p6.includes(p1)));

    BOOST_CHECK(is_false(p3.includes(p1)));
    BOOST_CHECK(is_false(p1.includes(p3)));

    BOOST_CHECK(is_true(p3.includes(p2)));
    BOOST_CHECK(is_false(p3.includes(p4)));
    BOOST_CHECK(is_false(p3.includes(p6)));

    BOOST_CHECK(is_false(p2.includes(p3)));
    BOOST_CHECK(is_false(p4.includes(p3)));
    BOOST_CHECK(is_false(p6.includes(p3)));

    StickyUnion<Polytope> Su1({p1,p3}), Su2({p2,p4,p6}), Su3({p8, p9}),
                        Su4({intersect(p2, p8), p1}), Su5({p7});
    SetsUnion<Polytope> res;

    // intersect the over-approximation of the union {p1, p3}
    res = intersect(Su1, Su1);

    // res is an over-approximation of Su1, thus their sets 
    // may be not included by any of the sets in Su1
    BOOST_CHECK(is_false(Su1.any_includes(res)));

    res = intersect(Su1, Su2);
    BOOST_CHECK(is_true(Su1.any_includes(res)));
    BOOST_CHECK(is_false(Su2.any_includes(res)));

    res = intersect(Su1, Su3);
    BOOST_CHECK(is_true(Su1.any_includes(res)));
    BOOST_CHECK(is_true(Su3.any_includes(res)));

    res = intersect(Su1, Su4);
    BOOST_CHECK(is_true(Su1.any_includes(res)));
    BOOST_CHECK(is_true(Su4.any_includes(res)));

    res = intersect(Su2, Su2);
    BOOST_CHECK(is_false(Su2.any_includes(res)));

    res = intersect(Su2, Su3);
    BOOST_CHECK(is_false(Su2.any_includes(res)));
    BOOST_CHECK(is_true(Su3.any_includes(res)));

    res = intersect(Su2, Su4);
    BOOST_CHECK(is_false(Su2.any_includes(res)));
    BOOST_CHECK(is_true(Su4.any_includes(res)));

    res = intersect(Su3, Su3);
    BOOST_CHECK(is_true(Su3.any_includes(res)));

    res = intersect(Su3, Su4);
    BOOST_CHECK(is_true(Su3.any_includes(res)));
    BOOST_CHECK(is_true(Su4.any_includes(res)));

    res = intersect(Su4, Su4);
    BOOST_CHECK(is_true(Su4.any_includes(res)));
}

BOOST_AUTO_TEST_CASE(test_includes_sticky_union)
{
    using namespace LinearAlgebra;
    using namespace LinearAlgebra::Dense;

    Matrix<double> A = {
        {1,0},
        {0,1}
    };

    Matrix<double> B = {
        {1}
    };

    Bundle b1(A,{3,1.5}, {5,3}), b2(A, {2,0}, {4,2}), b3(A,{0,1}, {2,3}),
           b4(A,{10,1}, {12,3}), b5(A,{11,0},{14,2}), b6(A,{13,2.5},{15,3.5}),
           b7(A,{0.5,1},{1.5,2.5}), b8(A,{4.5,1},{6,2}), b9(A,{-1,0},{-3,-4});

    Bundle c1(A,{0,1.5},{5.5,3}), c2(A,{10.5,1.5},{13,2});

    BOOST_CHECK(is_true(b9.is_empty()));

    BOOST_CHECK(is_true(b3.includes(b7)));

    BOOST_CHECK(is_false(are_disjoint(b1,b2)));
    BOOST_CHECK(is_true(are_disjoint(b1,b3)));
    BOOST_CHECK(is_true(are_disjoint(b1,b4)));
    BOOST_CHECK(is_true(are_disjoint(b1,b5)));
    BOOST_CHECK(is_true(are_disjoint(b1,b6)));

    BOOST_CHECK(is_false(are_disjoint(b2,b3)));
    BOOST_CHECK(is_true(are_disjoint(b2,b4)));
    BOOST_CHECK(is_true(are_disjoint(b2,b5)));
    BOOST_CHECK(is_true(are_disjoint(b2,b6)));

    BOOST_CHECK(is_true(are_disjoint(b3,b4)));
    BOOST_CHECK(is_true(are_disjoint(b3,b5)));
    BOOST_CHECK(is_true(are_disjoint(b3,b6)));

    BOOST_CHECK(is_false(are_disjoint(b4,b5)));
    BOOST_CHECK(is_true(are_disjoint(b4,b6)));

    BOOST_CHECK(is_true(are_disjoint(b5,b6)));

    std::list<Bundle> blist{b1,b2,b3,b4,b5,b6,b7,b8,b9};

    SetsUnion<Bundle> bsu(blist);

    BOOST_CHECK(is_false(bsu.any_includes(c1)));
    BOOST_CHECK(is_false(bsu.any_includes(c2)));

    BOOST_CHECK(bsu.size()==7);

    StickyUnion<Bundle> sticky_union = bsu;

    BOOST_CHECK(sticky_union.size()==7);
    BOOST_CHECK(sticky_union.number_of_classes()==3);

    for (auto l_it = std::begin(blist); l_it != std::end(blist); ++l_it) {
        BOOST_CHECK(is_true(sticky_union.any_includes(*l_it)));
    }

    BOOST_CHECK(is_true(sticky_union.includes(c1)));
    BOOST_CHECK(is_true(sticky_union.includes(c2)));
}
