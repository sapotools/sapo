#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE evolver

#include <boost/test/unit_test.hpp>

#include <sstream>

#include "Evolver.h"

#define APPROX_ERR 1e-14


inline bool epsilon_equivalent(const Polytope& P1,
                               const Polytope& P2, const double epsilon)
{
    return is_true((expand(P1, epsilon).includes(P2)) && 
                   (expand(P2, epsilon).includes(P1)));    
}

BOOST_AUTO_TEST_CASE(test_parametric_transform_bundle)
{
    using namespace SymbolicAlgebra;
    using namespace LinearAlgebra;

    Symbol<> s("s"), i("i"), r("r");
    Symbol<> alpha("alpha"), beta("beta");

    Dense::Matrix<double> pA{
        {1,0},
        {0,1}
    };

    std::map<Symbol<>, Expression<>> varDyn{
        {s, s-beta*s*i},
        {i, i+beta*s*i-alpha*i},
        {r, r+alpha*i}
    };

    Bundle pSet(pA, {0.5,0.1}, {0.6,0.2});

    Evolver<double> T(DiscreteSystem<double>(varDyn, {alpha, beta}));

    Dense::Matrix<double> rA{
        {1,0,0},
        {0,1,0},
        {0,0,1}
    };

    Bundle rSet(rA, {0,0,0}, {1,1,1});

    Bundle next = T(rSet, pSet);

    BOOST_CHECK(is_true(next==Bundle(rA, {0,0,0}, {1,0.7,1.6})));

    next = T(next, pSet);
    BOOST_CHECK(is_true(next==Bundle(rA, {0,0,0}, {1,0.49,2.02})));

    Dense::Matrix<double> errA1{
        {1,0,0,0},
        {0,1,0,0},
        {0,0,1,0},
        {0,0,0,1}
    };

    Bundle errSet1(errA1, {0,0,0,0}, {1,1,1,1});

    BOOST_REQUIRE_THROW(T(errSet1, pSet), std::domain_error);

    Dense::Matrix<double> errA2{
        {1,0},
        {0,1}
    };

    Bundle errSet2(errA2, {0,0}, {1,1});

    BOOST_REQUIRE_THROW(T(errSet2, pSet), std::domain_error);

    Dense::Matrix<double> errPA1{
        {1,0,0,0},
        {0,1,0,0},
        {0,0,1,0},
        {0,0,0,1}
    };

    Bundle errParam1(errPA1, {0,0,0,0}, {1,1,1,1});

    BOOST_REQUIRE_THROW(T(rSet, errParam1), std::domain_error);

    Dense::Matrix<double> errPA2{
        {1}
    };

    Bundle errParam2(errPA2, {0}, {1});

    BOOST_REQUIRE_THROW(T(rSet, errParam2), std::domain_error);
}

BOOST_AUTO_TEST_CASE(test_synthesis_bundle)
{
    using namespace SymbolicAlgebra;
    using namespace LinearAlgebra;

    Symbol<> s("s"), i("i"), r("r");
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

    Dense::Matrix<double> pA{
        {1,0},
        {0,1},
        {-1,0},
        {0,-1}
    };

    DiscreteSystem<double> f(vars, params, dyns); 

    std::shared_ptr<STL::Atom> atom1=std::make_shared<STL::Atom>(i-0.365),
                               atom2=std::make_shared<STL::Atom>(r-2);

    SetsUnion<Polytope> pSet(Polytope(pA, {0.6,0.2,-0.5,-0.1}));


    Dense::Matrix<double> rA{
        {1,0,0},
        {0,1,0},
        {0,0,1}
    };

    Bundle rSet(rA, {0,0,0}, {1,0.7,1.6});

    SetsUnion<Polytope> synthesized = synthesize(f, rSet, pSet, atom1);

    Polytope expected({{1,0},{0,-1},{-140,140}},{0.6,-0.1,-67.0});

    BOOST_CHECK(synthesized.size()==1);
    
    if (synthesized.size()==1) {
        // The following test does not work because of the approximation error
        // due to the double arithmetic
        // BOOST_CHECK(synthesized.begin()->simplify()==expected);
        BOOST_CHECK(epsilon_equivalent(*(synthesized.begin()), expected, APPROX_ERR));
    }

    synthesized = synthesize(f, rSet, pSet, atom2);
    expected = Polytope({{7,0},{0,1},{-1,0},{0,-1}},{4,0.2,-0.5,-0.1});

    BOOST_CHECK(synthesized.size()==1);
    
    if (synthesized.size()==1) {
        // The following test does not work because of the approximation error
        // due to the double arithmetic
        // BOOST_CHECK(synthesized.begin()->simplify()==expected);
        BOOST_CHECK(epsilon_equivalent(*(synthesized.begin()), expected, APPROX_ERR));
    }

    rSet = Bundle(rA, {0,0,0}, {1,0.65,1.55});
    pSet = Polytope(pA, {0.55,0.15,-0.5,-0.1});
    expected = Polytope({{1,0},{0,-1},{-130,130}},{0.55,-0.1,-57});

    synthesized = synthesize(f, rSet, pSet, atom1);

    BOOST_CHECK(synthesized.size()==1);
    
    if (synthesized.size()==1) {
        // The following test does not work because of the approximation error
        // due to the double arithmetic
        // BOOST_CHECK(synthesized.begin()->simplify()==expected);
        BOOST_CHECK(epsilon_equivalent(*(synthesized.begin()), expected, APPROX_ERR));
    }

    synthesized = synthesize(f, rSet, pSet, atom2);
    expected = Polytope({{1,0},{0,1},{-1,0},{0,-1}},{0.55,0.15,-0.5,-0.1});

    BOOST_CHECK(synthesized.size()==1);
    
    // The following test does not work because of the approximation error
    // due to the double arithmetic
    // BOOST_CHECK(synthesized.begin()->simplify()==expected);
    BOOST_CHECK(epsilon_equivalent(*(synthesized.begin()), expected, APPROX_ERR));
}