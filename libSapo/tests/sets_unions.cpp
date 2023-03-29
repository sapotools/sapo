#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE sets_union

#include <boost/test/unit_test.hpp>

#include <sstream>

#include "Polytope.h"
#include "Bundle.h"
#include "SetsUnion.h"

BOOST_AUTO_TEST_CASE(test_sets_union)
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

    Polytope<double> p1(A,{1,2,3,0,0,0}), p2(A,{3,4,5,-2,-3,-4}),
             p3(A,{1,2,3,-3,2,1}); 

    SetsUnion<Polytope<double>> Pu1, Pu2(p1), Pu3({p1, p2}), Pu4(p3);

    BOOST_CHECK(is_false(p1.is_empty()));
    BOOST_CHECK(is_false(p2.is_empty()));
    BOOST_CHECK(is_true(p3.is_empty()));

    BOOST_CHECK(is_true(Pu1.is_empty()));
    BOOST_CHECK(is_false(Pu2.is_empty()));
    BOOST_CHECK(is_false(Pu3.is_empty()));
    BOOST_CHECK(is_true(Pu4.is_empty()));

    BOOST_CHECK(Pu1.size()==0);
    BOOST_CHECK(Pu2.size()==1);
    BOOST_CHECK(Pu3.size()==2);
    BOOST_CHECK(Pu4.size()==0);

    BOOST_CHECK(Pu1.dim()==0);
    BOOST_CHECK(Pu2.dim()==3);
    BOOST_CHECK(Pu3.dim()==3);
    BOOST_CHECK(Pu4.dim()==0);
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

    Polytope<double> p1(A,{3,4,5,-2,-3,-4}), p2(A,{1,2,3,0,0,0}),
                     p3(A,{1,2,3,-3,2,1}), p4(A,{1,2,3,3,2,1}),
                     p5(A,{0.5,1.5,2.5,-0.5,-0.5,-0.5}); 

    SetsUnion<Polytope<double>> Pu1, Pu2(p1), Pu3({p1, p2}), Pu4(p3);

    BOOST_CHECK(is_false(p1.is_empty()));
    BOOST_CHECK(is_false(p2.is_empty()));
    BOOST_CHECK(is_true(p3.is_empty()));
    BOOST_CHECK(is_false(p4.is_empty()));
    BOOST_CHECK(is_false(p5.is_empty()));

    BOOST_CHECK(is_false(Pu1.any_includes(p4)));
    BOOST_CHECK(is_false(Pu2.any_includes(p4)));
    BOOST_CHECK(is_false(Pu3.any_includes(p4)));
    BOOST_CHECK(is_false(Pu4.any_includes(p4)));

    BOOST_CHECK(is_true(Pu1.any_includes(p3)));
    BOOST_CHECK(is_true(Pu2.any_includes(p3)));
    BOOST_CHECK(is_true(Pu3.any_includes(p3)));
    BOOST_CHECK(is_true(Pu4.any_includes(p3)));

    BOOST_CHECK(is_false(Pu1.any_includes(p5)));
    BOOST_CHECK(is_false(Pu2.any_includes(p5)));
    BOOST_CHECK(is_true(Pu3.any_includes(p5)));
    BOOST_CHECK(is_false(Pu4.any_includes(p5)));
}

BOOST_AUTO_TEST_CASE(test_copy_sets_union)
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

    Polytope<double> p1(A,{1,2,3,0,0,0}), p2(A,{3,4,5,-2,-3,-4}),
             p3(B,{1,1,3,3,0,2,4,2,1,1,3,3}),
             p4(C,{1,1,0,0});

    BOOST_CHECK(is_false(p1.is_empty()));
    BOOST_CHECK(is_false(p2.is_empty()));
    BOOST_CHECK(is_false(p3.is_empty()));
    BOOST_CHECK(is_false(p4.is_empty()));

    SetsUnion<Polytope<double>> Pu1({p1,p2,p3}), Pu2(p4), Pu3, res;
    BOOST_CHECK(Pu1.size()==3);

    res = Pu3;

    BOOST_CHECK(is_true(res.is_empty()));

    for (const auto &Pu:{Pu1, Pu2, Pu3}) {
        res = Pu;
        BOOST_CHECK(is_true(res.is_empty()==Pu.is_empty()));
        BOOST_CHECK(res.size()==Pu.size());

        auto res_it = std::end(res);
        auto p_it = std::end(Pu);

        for (; res_it != std::end(res); ++res_it, ++p_it) {
            BOOST_CHECK(is_true(*res_it == *p_it));
        }
    }
}

BOOST_AUTO_TEST_CASE(test_add_sets_union)
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

    Polytope<double> p1(A,{3,4,5,3,2,1}), p2(A,{1,2,3,0,0,0}), 
             p3(A,{4,4,5,0,0,0}), p4(A,{3,4,5,-2,-3,-4}), 
             p5(A,{1,2,3,-3,2,1}), 
             p6(B,{1,1,3,3,0,2,4,2,1,1,3,3}),
             p7(C,{1,1,0,0});

    BOOST_CHECK(is_false(p1.is_empty()));
    BOOST_CHECK(is_false(p2.is_empty()));
    BOOST_CHECK(is_false(p3.is_empty()));
    BOOST_CHECK(is_false(p4.is_empty()));
    BOOST_CHECK(is_true(p5.is_empty()));
    BOOST_CHECK(is_false(p6.is_empty()));
    BOOST_CHECK(is_false(p7.is_empty()));

    BOOST_CHECK(is_true(p1.includes(p2)));
    BOOST_CHECK(is_false(p3.includes(p1)));
    BOOST_CHECK(is_false(p1.includes(p3)));
    BOOST_CHECK(is_true(p1.includes(p4)));
    BOOST_CHECK(is_false(p4.includes(p1)));
    BOOST_CHECK(is_true(p3.includes(p2)));
    BOOST_CHECK(is_true(p3.includes(p4)));
    BOOST_CHECK(is_false(p3.includes(p6)));
    BOOST_CHECK(is_false(p2.includes(p3)));
    BOOST_CHECK(is_true(p1.includes(p6)));
    BOOST_CHECK(is_false(p6.includes(p1)));
    BOOST_CHECK(is_false(p6.includes(p2)));
    BOOST_CHECK(is_false(p6.includes(p3)));
    BOOST_CHECK(is_false(p6.includes(p4)));

    SetsUnion<Polytope<double>> Pu1({p2,p4});

    BOOST_CHECK(Pu1.size()==2);

    Pu1.add(p5);

    BOOST_CHECK(Pu1.size()==2);

    Pu1.add(p6);

    BOOST_CHECK(Pu1.size()==3);

    Pu1.add(p3);
    
    BOOST_CHECK(Pu1.size()==2);

    Pu1 = SetsUnion<Polytope<double>>({p2,p4,p6});

    BOOST_CHECK(Pu1.size()==3);

    Pu1.add(p1);

    BOOST_CHECK(Pu1.size()==1);

    BOOST_REQUIRE_THROW(Pu1.add(p7), std::domain_error);
}

BOOST_AUTO_TEST_CASE(test_make_union_sets_union)
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

    Polytope<double> p1(A,{3,4,6,3,2,1}), p2(A,{1,2,3,0,0,0}), 
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

    SetsUnion<Polytope<double>> Pu1({p1,p3}), Pu2({p2,p4,p6}), Pu3({p8, p9}),
                        Pu4({intersect(p2, p8), p1}), Pu5({p7}), res1, res2;


    for (auto it=std::begin(Pu2); it != std::end(Pu2); ++it) {
        BOOST_CHECK(is_true(Pu1.any_includes(*it)));
    }
    
    for (auto it=std::begin(Pu3); it != std::end(Pu3); ++it) {
        BOOST_CHECK(is_true(Pu1.any_includes(*it)));
    }

    for (auto it=std::begin(Pu2); it != std::end(Pu2); ++it) {
        BOOST_CHECK(is_false(Pu3.any_includes(*it)));
    }

    BOOST_CHECK(Pu2.size()==3);

    res1 = make_union(Pu2, Pu3);
    res2 = make_union(Pu3, Pu2);

    for (auto it=std::begin(Pu2); it != std::end(Pu2); ++it) {
        BOOST_CHECK(is_true(res1.any_includes(*it)));
    }

    for (auto it=std::begin(Pu3); it != std::end(Pu3); ++it) {
        BOOST_CHECK(is_true(res1.any_includes(*it)));
    }

    for (auto it=std::begin(res2); it != std::end(res2); ++it) {
        BOOST_CHECK(is_true(res1.any_includes(*it)));
    }

    for (auto it=std::begin(res1); it != std::end(res1); ++it) {
        BOOST_CHECK(is_true(res2.any_includes(*it)));
    }

    BOOST_CHECK(res1.size()==5);

    res1 = make_union(res1, Pu4);
    res2 = make_union(Pu4, res2);

    BOOST_CHECK(res1.size()==2);
    BOOST_CHECK(res1.size()==res2.size());

    for (auto it=std::begin(Pu2); it != std::end(Pu2); ++it) {
        BOOST_CHECK(is_true(res1.any_includes(*it)));
    }

    for (auto it=std::begin(Pu3); it != std::end(Pu3); ++it) {
        BOOST_CHECK(is_true(res1.any_includes(*it)));
    }

    for (auto it=std::begin(Pu4); it != std::end(Pu4); ++it) {
        BOOST_CHECK(is_true(res1.any_includes(*it)));
    }

    for (auto it=std::begin(res2); it != std::end(res2); ++it) {
        BOOST_CHECK(is_true(res1.any_includes(*it)));
    }

    for (auto it=std::begin(res1); it != std::end(res1); ++it) {
        BOOST_CHECK(is_true(res2.any_includes(*it)));
    }

    res1 = make_union(Pu1, Pu2);
    res2 = make_union(Pu2, Pu1);

    BOOST_CHECK(res1.size()==2);
    BOOST_CHECK(res1.size()==res2.size());

    for (auto it=std::begin(Pu2); it != std::end(Pu2); ++it) {
        BOOST_CHECK(is_true(res1.any_includes(*it)));
    }

    for (auto it=std::begin(Pu1); it != std::end(Pu1); ++it) {
        BOOST_CHECK(is_true(res1.any_includes(*it)));
    }

    for (auto it=std::begin(res2); it != std::end(res2); ++it) {
        BOOST_CHECK(is_true(res1.any_includes(*it)));
    }

    for (auto it=std::begin(res1); it != std::end(res1); ++it) {
        BOOST_CHECK(is_true(res2.any_includes(*it)));
    }

    BOOST_REQUIRE_THROW(res1.update(Pu5), std::domain_error);
}

BOOST_AUTO_TEST_CASE(test_update_sets_union)
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

    Polytope<double> p1(A,{3,4,6,3,2,1}), p2(A,{1,2,3,0,0,0}), 
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

    SetsUnion<Polytope<double>> Pu1({p1,p3}), Pu2({p2,p4,p6}), Pu3({p8, p9}),
                        Pu4({intersect(p2, p8), p1}), Pu5({p7}), res;


    for (auto it=std::begin(Pu2); it != std::end(Pu2); ++it) {
        BOOST_CHECK(is_true(Pu1.any_includes(*it)));
    }
    
    for (auto it=std::begin(Pu3); it != std::end(Pu3); ++it) {
        BOOST_CHECK(is_true(Pu1.any_includes(*it)));
    }

    for (auto it=std::begin(Pu2); it != std::end(Pu2); ++it) {
        BOOST_CHECK(is_false(Pu3.any_includes(*it)));
    }

    BOOST_CHECK(Pu2.size()==3);

    res = Pu2;

    res.update(Pu3);

    for (auto it=std::begin(Pu2); it != std::end(Pu2); ++it) {
        BOOST_CHECK(is_true(res.any_includes(*it)));
    }

    for (auto it=std::begin(Pu3); it != std::end(Pu3); ++it) {
        BOOST_CHECK(is_true(res.any_includes(*it)));
    }

    BOOST_CHECK(res.size()==5);

    res.update(Pu4);

    BOOST_CHECK(res.size()==2);

    for (auto it=std::begin(Pu2); it != std::end(Pu2); ++it) {
        BOOST_CHECK(is_true(res.any_includes(*it)));
    }

    for (auto it=std::begin(Pu3); it != std::end(Pu3); ++it) {
        BOOST_CHECK(is_true(res.any_includes(*it)));
    }

    for (auto it=std::begin(Pu4); it != std::end(Pu4); ++it) {
        BOOST_CHECK(is_true(res.any_includes(*it)));
    }

    res = Pu2;

    res.update(Pu1);

     for (auto it=std::begin(Pu2); it != std::end(Pu2); ++it) {
        BOOST_CHECK(is_true(res.any_includes(*it)));
    }

    for (auto it=std::begin(Pu1); it != std::end(Pu1); ++it) {
        BOOST_CHECK(is_true(res.any_includes(*it)));
    }
    
    BOOST_CHECK(res.size()==2);

    BOOST_REQUIRE_THROW(res.update(Pu5), std::domain_error);
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

    Polytope<double> p1(A,{3,4,6,3,2,1}), p2(A,{1,2,3,0,0,0}), 
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

    SetsUnion<Polytope<double>> Pu1({p1,p3}), Pu2({p2,p4,p6}), Pu3({p8, p9}),
                        Pu4({intersect(p2, p8), p1}), Pu5({p7}), res;

    res = intersect(Pu1, Pu1);
    for (auto it=std::begin(res); it != std::end(res); ++it) {
        BOOST_CHECK(is_true(Pu1.any_includes(*it)));
    }
    for (auto it=std::begin(Pu1); it != std::end(Pu1); ++it) {
        BOOST_CHECK(is_true(res.any_includes(*it)));
    }

    res = intersect(Pu1, Pu2);

    for (auto it=std::begin(res); it != std::end(res); ++it) {
        BOOST_CHECK(is_true(Pu1.any_includes(*it))&&
                    is_true(Pu2.any_includes(*it)));
    }

    res = intersect(Pu1, Pu3);

    for (auto it=std::begin(res); it != std::end(res); ++it) {
        BOOST_CHECK(is_true(Pu1.any_includes(*it))&&
                    is_true(Pu3.any_includes(*it)));
    }

    res = intersect(Pu1, Pu4);

    for (auto it=std::begin(res); it != std::end(res); ++it) {
        BOOST_CHECK(is_true(Pu1.any_includes(*it))&&
                    is_true(Pu4.any_includes(*it)));
    }

    res = intersect(Pu2, Pu2);
    for (auto it=std::begin(res); it != std::end(res); ++it) {
        BOOST_CHECK(is_true(Pu2.any_includes(*it)));
    }
    for (auto it=std::begin(Pu2); it != std::end(Pu2); ++it) {
        BOOST_CHECK(is_true(res.any_includes(*it)));
    }

    res = intersect(Pu2, Pu3);

    for (auto it=std::begin(res); it != std::end(res); ++it) {
        BOOST_CHECK(is_true(Pu2.any_includes(*it)));
        BOOST_CHECK(is_true(Pu3.any_includes(*it)));
    }

    res = intersect(Pu2, Pu4);

    for (auto it=std::begin(res); it != std::end(res); ++it) {
        BOOST_CHECK(is_true(Pu2.any_includes(*it)));
        BOOST_CHECK(is_true(Pu4.any_includes(*it)));
    }

    res = intersect(Pu3, Pu3);
    for (auto it=std::begin(res); it != std::end(res); ++it) {
        BOOST_CHECK(is_true(Pu3.any_includes(*it)));
    }
    for (auto it=std::begin(Pu3); it != std::end(Pu3); ++it) {
        BOOST_CHECK(is_true(res.any_includes(*it)));
    }

    res = intersect(Pu3, Pu4);

    for (auto it=std::begin(res); it != std::end(res); ++it) {
        BOOST_CHECK(is_true(Pu3.any_includes(*it)));
        BOOST_CHECK(is_true(Pu4.any_includes(*it)));
    }

    res = intersect(Pu4, Pu4);
    for (auto it=std::begin(res); it != std::end(res); ++it) {
        BOOST_CHECK(is_true(Pu4.any_includes(*it)));
    }
    for (auto it=std::begin(Pu4); it != std::end(Pu4); ++it) {
        BOOST_CHECK(is_true(res.any_includes(*it)));
    }

    BOOST_CHECK(is_true(are_disjoint(Pu2, SetsUnion<Polytope<double>>())));
    BOOST_CHECK(is_true(are_disjoint(SetsUnion<Polytope<double>>(), Pu2)));
}

BOOST_AUTO_TEST_CASE(test_includes_sets_union)
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

    Bundle b1(A,{3,1.5}, {5,3}), b2(A, {2,0}, {4,2.5}), b3(A,{0,1}, {2,3}),
           b4(A,{10,1}, {12,3}), b5(A,{11,0},{14,2}), b6(A,{13,2.5},{15,3.5}),
           b7(A,{0.5,1},{1.5,2.5}), b8(A,{4.5,1},{6,2}), b9(A,{-1,0},{-3,-4});

    Bundle c1(A,{0,1.5},{5.5,2}), c2(A,{10.5,1},{13,1.5}), c3(A, {0,1},{5.5,3});

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

    std::list<Bundle<double>> blist{b1,b2,b3,b4,b5,b6,b7,b8,b9};

    SetsUnion<Bundle<double>> bsu(blist);
    SetsUnion<Bundle<double>> bsu2({c1, c2});

    BOOST_CHECK(is_false(bsu.any_includes(c1)));
    BOOST_CHECK(is_false(bsu.any_includes(c2)));
    BOOST_CHECK(is_false(bsu.any_includes(c3)));

    BOOST_CHECK(is_true(bsu.includes(c1)));
    /*
    BOOST_CHECK(is_true(bsu.includes(c2)));
    BOOST_CHECK(is_false(bsu.includes(c3)));

    BOOST_CHECK(is_true(bsu.includes(bsu2)));
    BOOST_CHECK(is_false(bsu2.includes(bsu)));
    */
}
