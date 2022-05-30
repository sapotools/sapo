#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE sets_union

#include <boost/test/unit_test.hpp>

#include <sstream>

#include "Polytope.h"
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

    Polytope p1(A,{1,2,3,0,0,0}), p2(A,{3,4,5,-2,-3,-4}),
             p3(A,{1,2,3,-3,2,1}); 

    SetsUnion<Polytope> Pu1, Pu2(p1), Pu3({p1, p2}), Pu4(p3);

    BOOST_CHECK(!p1.is_empty());
    BOOST_CHECK(!p2.is_empty());
    BOOST_CHECK(p3.is_empty());

    BOOST_CHECK(Pu1.is_empty());
    BOOST_CHECK(!Pu2.is_empty());
    BOOST_CHECK(!Pu3.is_empty());
    BOOST_CHECK(Pu4.is_empty());

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

    Polytope p1(A,{3,4,5,-2,-3,-4}), p2(A,{1,2,3,0,0,0}),
             p3(A,{1,2,3,-3,2,1}), p4(A,{1,2,3,3,2,1}),
             p5(A,{0.5,1.5,2.5,-0.5,-0.5,-0.5}); 

    SetsUnion<Polytope> Pu1, Pu2(p1), Pu3({p1, p2}), Pu4(p3);

    BOOST_CHECK(!p1.is_empty());
    BOOST_CHECK(!p2.is_empty());
    BOOST_CHECK(p3.is_empty());
    BOOST_CHECK(!p4.is_empty());
    BOOST_CHECK(!p5.is_empty());

    BOOST_CHECK(!Pu1.any_includes(p4));
    BOOST_CHECK(!Pu2.any_includes(p4));
    BOOST_CHECK(!Pu3.any_includes(p4));
    BOOST_CHECK(!Pu4.any_includes(p4));

    BOOST_CHECK(Pu1.any_includes(p3));
    BOOST_CHECK(Pu2.any_includes(p3));
    BOOST_CHECK(Pu3.any_includes(p3));
    BOOST_CHECK(Pu4.any_includes(p3));

    BOOST_CHECK(!Pu1.any_includes(p5));
    BOOST_CHECK(!Pu2.any_includes(p5));
    BOOST_CHECK(Pu3.any_includes(p5));
    BOOST_CHECK(!Pu4.any_includes(p5));
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

    Matrix<double> C = {
        {1,0},
        {0,1},
        {-1,0},
        {0,-1}
    };

    Polytope p1(A,{1,2,3,0,0,0}), p2(A,{3,4,5,-2,-3,-4}),
             p3(B,{1,1,3,3,0,2,4,2,1,1,3,3}),
             p4(C,{1,1,0,0});

    BOOST_CHECK(!p1.is_empty());
    BOOST_CHECK(!p2.is_empty());
    BOOST_CHECK(!p3.is_empty());
    BOOST_CHECK(!p4.is_empty());

    SetsUnion<Polytope> Pu1({p1,p2,p3}), Pu2(p4), Pu3, res;
    BOOST_CHECK(Pu1.size()==3);

    res = Pu3;

    BOOST_CHECK(res.is_empty());

    for (const auto &Pu:{Pu1, Pu2, Pu3}) {
        res = Pu;
        BOOST_CHECK(res.is_empty()==Pu.is_empty());
        BOOST_CHECK(res.size()==Pu.size());

        auto res_it = std::end(res);
        auto p_it = std::end(Pu);

        for (; res_it != std::end(res); ++res_it, ++p_it) {
            BOOST_CHECK(*res_it == *p_it);
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

    Matrix<double> C = {
        {1,0},
        {0,1},
        {-1,0},
        {0,-1}
    };

    Polytope p1(A,{3,4,5,3,2,1}), p2(A,{1,2,3,0,0,0}), 
             p3(A,{4,4,5,0,0,0}), p4(A,{3,4,5,-2,-3,-4}), 
             p5(A,{1,2,3,-3,2,1}), 
             p6(B,{1,1,3,3,0,2,4,2,1,1,3,3}),
             p7(C,{1,1,0,0});

    BOOST_CHECK(!p1.is_empty());
    BOOST_CHECK(!p2.is_empty());
    BOOST_CHECK(!p3.is_empty());
    BOOST_CHECK(!p4.is_empty());
    BOOST_CHECK(p5.is_empty());
    BOOST_CHECK(!p6.is_empty());
    BOOST_CHECK(!p7.is_empty());

    BOOST_CHECK(p1.includes(p2));
    BOOST_CHECK(!p3.includes(p1));
    BOOST_CHECK(!p1.includes(p3));
    BOOST_CHECK(p1.includes(p4));
    BOOST_CHECK(!p4.includes(p1));
    BOOST_CHECK(p3.includes(p2));
    BOOST_CHECK(p3.includes(p4));
    BOOST_CHECK(!p3.includes(p6));
    BOOST_CHECK(!p2.includes(p3));
    BOOST_CHECK(p1.includes(p6));
    BOOST_CHECK(!p6.includes(p1));
    BOOST_CHECK(!p6.includes(p2));
    BOOST_CHECK(!p6.includes(p3));
    BOOST_CHECK(!p6.includes(p4));

    SetsUnion<Polytope> Pu1({p2,p4});

    BOOST_CHECK(Pu1.size()==2);

    Pu1.add(p5);

    BOOST_CHECK(Pu1.size()==2);

    Pu1.add(p6);

    BOOST_CHECK(Pu1.size()==3);

    Pu1.add(p3);
    
    BOOST_CHECK(Pu1.size()==2);

    Pu1 = SetsUnion<Polytope>({p2,p4,p6});

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

    BOOST_CHECK(!p1.is_empty());
    BOOST_CHECK(!p2.is_empty());
    BOOST_CHECK(!p3.is_empty());
    BOOST_CHECK(!p4.is_empty());
    BOOST_CHECK(p5.is_empty());
    BOOST_CHECK(!p6.is_empty());
    BOOST_CHECK(!p7.is_empty());
    BOOST_CHECK(!p8.is_empty());
    BOOST_CHECK(!p9.is_empty());

    BOOST_CHECK(p1.includes(p2));
    BOOST_CHECK(p1.includes(p4));
    BOOST_CHECK(p1.includes(p6));

    BOOST_CHECK(!p2.includes(p1));
    BOOST_CHECK(!p4.includes(p1));
    BOOST_CHECK(!p6.includes(p1));

    BOOST_CHECK(!p3.includes(p1));
    BOOST_CHECK(!p1.includes(p3));

    BOOST_CHECK(p3.includes(p2));
    BOOST_CHECK(!p3.includes(p4));
    BOOST_CHECK(!p3.includes(p6));

    BOOST_CHECK(!p2.includes(p3));
    BOOST_CHECK(!p4.includes(p3));
    BOOST_CHECK(!p6.includes(p3));

    SetsUnion<Polytope> Pu1({p1,p3}), Pu2({p2,p4,p6}), Pu3({p8, p9}),
                        Pu4({intersect(p2, p8), p1}), Pu5({p7}), res1, res2;


    for (auto it=std::begin(Pu2); it != std::end(Pu2); ++it) {
        BOOST_CHECK(Pu1.any_includes(*it));
    }
    
    for (auto it=std::begin(Pu3); it != std::end(Pu3); ++it) {
        BOOST_CHECK(Pu1.any_includes(*it));
    }

    for (auto it=std::begin(Pu2); it != std::end(Pu2); ++it) {
        BOOST_CHECK(!Pu3.any_includes(*it));
    }

    BOOST_CHECK(Pu2.size()==3);

    res1 = make_union(Pu2, Pu3);
    res2 = make_union(Pu3, Pu2);

    for (auto it=std::begin(Pu2); it != std::end(Pu2); ++it) {
        BOOST_CHECK(res1.any_includes(*it));
    }

    for (auto it=std::begin(Pu3); it != std::end(Pu3); ++it) {
        BOOST_CHECK(res1.any_includes(*it));
    }

    for (auto it=std::begin(res2); it != std::end(res2); ++it) {
        BOOST_CHECK(res1.any_includes(*it));
    }

    for (auto it=std::begin(res1); it != std::end(res1); ++it) {
        BOOST_CHECK(res2.any_includes(*it));
    }

    BOOST_CHECK(res1.size()==5);

    res1 = make_union(res1, Pu4);
    res2 = make_union(Pu4, res2);

    BOOST_CHECK(res1.size()==2);
    BOOST_CHECK(res1.size()==res2.size());

    for (auto it=std::begin(Pu2); it != std::end(Pu2); ++it) {
        BOOST_CHECK(res1.any_includes(*it));
    }

    for (auto it=std::begin(Pu3); it != std::end(Pu3); ++it) {
        BOOST_CHECK(res1.any_includes(*it));
    }

    for (auto it=std::begin(Pu4); it != std::end(Pu4); ++it) {
        BOOST_CHECK(res1.any_includes(*it));
    }

    for (auto it=std::begin(res2); it != std::end(res2); ++it) {
        BOOST_CHECK(res1.any_includes(*it));
    }

    for (auto it=std::begin(res1); it != std::end(res1); ++it) {
        BOOST_CHECK(res2.any_includes(*it));
    }

    res1 = make_union(Pu1, Pu2);
    res2 = make_union(Pu2, Pu1);

    BOOST_CHECK(res1.size()==2);
    BOOST_CHECK(res1.size()==res2.size());

    for (auto it=std::begin(Pu2); it != std::end(Pu2); ++it) {
        BOOST_CHECK(res1.any_includes(*it));
    }

    for (auto it=std::begin(Pu1); it != std::end(Pu1); ++it) {
        BOOST_CHECK(res1.any_includes(*it));
    }

    for (auto it=std::begin(res2); it != std::end(res2); ++it) {
        BOOST_CHECK(res1.any_includes(*it));
    }

    for (auto it=std::begin(res1); it != std::end(res1); ++it) {
        BOOST_CHECK(res2.any_includes(*it));
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

    BOOST_CHECK(!p1.is_empty());
    BOOST_CHECK(!p2.is_empty());
    BOOST_CHECK(!p3.is_empty());
    BOOST_CHECK(!p4.is_empty());
    BOOST_CHECK(p5.is_empty());
    BOOST_CHECK(!p6.is_empty());
    BOOST_CHECK(!p7.is_empty());
    BOOST_CHECK(!p8.is_empty());
    BOOST_CHECK(!p9.is_empty());

    BOOST_CHECK(p1.includes(p2));
    BOOST_CHECK(p1.includes(p4));
    BOOST_CHECK(p1.includes(p6));

    BOOST_CHECK(!p2.includes(p1));
    BOOST_CHECK(!p4.includes(p1));
    BOOST_CHECK(!p6.includes(p1));

    BOOST_CHECK(!p3.includes(p1));
    BOOST_CHECK(!p1.includes(p3));

    BOOST_CHECK(p3.includes(p2));
    BOOST_CHECK(!p3.includes(p4));
    BOOST_CHECK(!p3.includes(p6));

    BOOST_CHECK(!p2.includes(p3));
    BOOST_CHECK(!p4.includes(p3));
    BOOST_CHECK(!p6.includes(p3));

    SetsUnion<Polytope> Pu1({p1,p3}), Pu2({p2,p4,p6}), Pu3({p8, p9}),
                        Pu4({intersect(p2, p8), p1}), Pu5({p7}), res;


    for (auto it=std::begin(Pu2); it != std::end(Pu2); ++it) {
        BOOST_CHECK(Pu1.any_includes(*it));
    }
    
    for (auto it=std::begin(Pu3); it != std::end(Pu3); ++it) {
        BOOST_CHECK(Pu1.any_includes(*it));
    }

    for (auto it=std::begin(Pu2); it != std::end(Pu2); ++it) {
        BOOST_CHECK(!Pu3.any_includes(*it));
    }

    BOOST_CHECK(Pu2.size()==3);

    res = Pu2;

    res.update(Pu3);

    for (auto it=std::begin(Pu2); it != std::end(Pu2); ++it) {
        BOOST_CHECK(res.any_includes(*it));
    }

    for (auto it=std::begin(Pu3); it != std::end(Pu3); ++it) {
        BOOST_CHECK(res.any_includes(*it));
    }

    BOOST_CHECK(res.size()==5);

    res.update(Pu4);

    BOOST_CHECK(res.size()==2);

    for (auto it=std::begin(Pu2); it != std::end(Pu2); ++it) {
        BOOST_CHECK(res.any_includes(*it));
    }

    for (auto it=std::begin(Pu3); it != std::end(Pu3); ++it) {
        BOOST_CHECK(res.any_includes(*it));
    }

    for (auto it=std::begin(Pu4); it != std::end(Pu4); ++it) {
        BOOST_CHECK(res.any_includes(*it));
    }

    res = Pu2;

    res.update(Pu1);

     for (auto it=std::begin(Pu2); it != std::end(Pu2); ++it) {
        BOOST_CHECK(res.any_includes(*it));
    }

    for (auto it=std::begin(Pu1); it != std::end(Pu1); ++it) {
        BOOST_CHECK(res.any_includes(*it));
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

    BOOST_CHECK(!p1.is_empty());
    BOOST_CHECK(!p2.is_empty());
    BOOST_CHECK(!p3.is_empty());
    BOOST_CHECK(!p4.is_empty());
    BOOST_CHECK(p5.is_empty());
    BOOST_CHECK(!p6.is_empty());
    BOOST_CHECK(!p7.is_empty());
    BOOST_CHECK(!p8.is_empty());
    BOOST_CHECK(!p9.is_empty());

    BOOST_CHECK(p1.includes(p2));
    BOOST_CHECK(p1.includes(p4));
    BOOST_CHECK(p1.includes(p6));

    BOOST_CHECK(!p2.includes(p1));
    BOOST_CHECK(!p4.includes(p1));
    BOOST_CHECK(!p6.includes(p1));

    BOOST_CHECK(!p3.includes(p1));
    BOOST_CHECK(!p1.includes(p3));

    BOOST_CHECK(p3.includes(p2));
    BOOST_CHECK(!p3.includes(p4));
    BOOST_CHECK(!p3.includes(p6));

    BOOST_CHECK(!p2.includes(p3));
    BOOST_CHECK(!p4.includes(p3));
    BOOST_CHECK(!p6.includes(p3));

    SetsUnion<Polytope> Pu1({p1,p3}), Pu2({p2,p4,p6}), Pu3({p8, p9}),
                        Pu4({intersect(p2, p8), p1}), Pu5({p7}), res;

    res = intersect(Pu1, Pu1);
    for (auto it=std::begin(res); it != std::end(res); ++it) {
        BOOST_CHECK(Pu1.any_includes(*it));
    }
    for (auto it=std::begin(Pu1); it != std::end(Pu1); ++it) {
        BOOST_CHECK(res.any_includes(*it));
    }

    res = intersect(Pu1, Pu2);

    for (auto it=std::begin(res); it != std::end(res); ++it) {
        BOOST_CHECK(Pu1.any_includes(*it)&&Pu2.any_includes(*it));
    }

    res = intersect(Pu1, Pu3);

    for (auto it=std::begin(res); it != std::end(res); ++it) {
        BOOST_CHECK(Pu1.any_includes(*it)&&Pu3.any_includes(*it));
    }

    res = intersect(Pu1, Pu4);

    for (auto it=std::begin(res); it != std::end(res); ++it) {
        BOOST_CHECK(Pu1.any_includes(*it)&&Pu4.any_includes(*it));
    }

    res = intersect(Pu2, Pu2);
    for (auto it=std::begin(res); it != std::end(res); ++it) {
        BOOST_CHECK(Pu2.any_includes(*it));
    }
    for (auto it=std::begin(Pu2); it != std::end(Pu2); ++it) {
        BOOST_CHECK(res.any_includes(*it));
    }

    res = intersect(Pu2, Pu3);

    for (auto it=std::begin(res); it != std::end(res); ++it) {
        BOOST_CHECK(Pu2.any_includes(*it)&&Pu3.any_includes(*it));
    }

    res = intersect(Pu2, Pu4);

    for (auto it=std::begin(res); it != std::end(res); ++it) {
        BOOST_CHECK(Pu2.any_includes(*it)&&Pu4.any_includes(*it));
    }

    res = intersect(Pu3, Pu3);
    for (auto it=std::begin(res); it != std::end(res); ++it) {
        BOOST_CHECK(Pu3.any_includes(*it));
    }
    for (auto it=std::begin(Pu3); it != std::end(Pu3); ++it) {
        BOOST_CHECK(res.any_includes(*it));
    }

    res = intersect(Pu3, Pu4);

    for (auto it=std::begin(res); it != std::end(res); ++it) {
        BOOST_CHECK(Pu3.any_includes(*it)&&
                    Pu4.any_includes(*it));
    }

    res = intersect(Pu4, Pu4);
    for (auto it=std::begin(res); it != std::end(res); ++it) {
        BOOST_CHECK(Pu4.any_includes(*it));
    }
    for (auto it=std::begin(Pu4); it != std::end(Pu4); ++it) {
        BOOST_CHECK(res.any_includes(*it));
    }

    BOOST_CHECK(intersect(Pu2, SetsUnion<Polytope>()).is_empty());
}