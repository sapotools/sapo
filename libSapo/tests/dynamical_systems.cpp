#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE dynamics

#include <boost/test/unit_test.hpp>

#include <sstream>

#include "DynamicalSystem.h"

#define APPROX_ERR 1e-14

BOOST_AUTO_TEST_CASE(test_transform_bundle)
{
    using namespace SymbolicAlgebra;
    using namespace LinearAlgebra;

    Symbol<> s("s"), i("i"), r("r");

    // variable vector
    std::vector<Symbol<>> vars{s, i, r};

    DynamicalSystem<double> f(vars, {
                            s-0.1*s*i,          // s's law
                            i+0.1*s*i-0.5*i,    // i's law
                            r+0.5*i             // r's law
                        });

    Dense::Matrix<double> rA{
        {1,0,0},
        {0,1,0},
        {0,0,1}
    };

    Bundle rSet(rA, {0,0,0}, {1,1,1});

    Bundle next = f.transform(rSet);
    BOOST_CHECK(next==Bundle(rA, {0,0,0}, {1,0.6,1.5}));

    next = f.transform(next);
    BOOST_CHECK(next==Bundle(rA, {0,0,0}, {1,0.36,1.8}));

    f = DynamicalSystem<double>(vars, {
                            s-0.2*s*i,          // s's law
                            i+0.2*s*i-0.6*i,    // i's law
                            r+0.6*i             // r's law
                        });

    next = f.transform(rSet);
    BOOST_CHECK(next==Bundle(rA, {0,0,0}, {1,0.6,1.6}));

    next = f.transform(next);
    BOOST_CHECK(next==Bundle(rA, {0,0,0}, {1,0.36,1.96}));
}

BOOST_AUTO_TEST_CASE(test_parametric_transform_bundle)
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
        {0,1}
    };

    Bundle pSet(pA, {0.5,0.1}, {0.6,0.2});

    DynamicalSystem<double> f(vars, params, dyns);

    Dense::Matrix<double> rA{
        {1,0,0},
        {0,1,0},
        {0,0,1}
    };

    Bundle rSet(rA, {0,0,0}, {1,1,1});

    Bundle next = f.transform(rSet, pSet);

    BOOST_CHECK(next==Bundle(rA, {0,0,0}, {1,0.7,1.6}));

    next = f.transform(next, pSet);
    BOOST_CHECK(next==Bundle(rA, {0,0,0}, {1,0.49,2.02}));
}

inline bool epsilon_equivalent(const Polytope& synthesized,
                               const Polytope& expected, const double epsilon)
{
    return ((expand(synthesized, epsilon).contains(expected)) && 
            (expand(expected, epsilon).contains(synthesized)));    
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

    DynamicalSystem<double> f(vars, params, dyns); 

    std::shared_ptr<STL::Atom> atom1=std::make_shared<STL::Atom>(i-0.365),
                               atom2=std::make_shared<STL::Atom>(r-2);

    PolytopesUnion pSet(Polytope(pA, {0.6,0.2,-0.5,-0.1}));


    Dense::Matrix<double> rA{
        {1,0,0},
        {0,1,0},
        {0,0,1}
    };

    Bundle rSet(rA, {0,0,0}, {1,0.7,1.6});

    PolytopesUnion synthesized = synthesize(f, rSet, pSet, atom1);

    Polytope expected({{1,0},{0,-1},{-140,140}},{0.6,-0.1,-67.0});

    BOOST_CHECK(synthesized.size()==1);
    
    // The following test does not work because of the approximation error
    // due to the double arithmetic
    // BOOST_CHECK(synthesized.begin()->simplify()==expected);
    BOOST_CHECK(epsilon_equivalent(*(synthesized.begin()), expected, APPROX_ERR));

    synthesized = synthesize(f, rSet, pSet, atom2);
    expected = Polytope({{7,0},{0,1},{-1,0},{0,-1}},{4,0.2,-0.5,-0.1});

    BOOST_CHECK(synthesized.size()==1);
    
    // The following test does not work because of the approximation error
    // due to the double arithmetic
    // BOOST_CHECK(synthesized.begin()->simplify()==expected);
    BOOST_CHECK(epsilon_equivalent(*(synthesized.begin()), expected, APPROX_ERR));

    rSet = Bundle(rA, {0,0,0}, {1,0.65,1.55});
    pSet = Polytope(pA, {0.55,0.15,-0.5,-0.1});
    expected = Polytope({{1,0},{0,-1},{-130,130}},{0.55,-0.1,-57});

    synthesized = synthesize(f, rSet, pSet, atom1);

    BOOST_CHECK(synthesized.size()==1);
    
    // The following test does not work because of the approximation error
    // due to the double arithmetic
    // BOOST_CHECK(synthesized.begin()->simplify()==expected);
    BOOST_CHECK(epsilon_equivalent(*(synthesized.begin()), expected, APPROX_ERR));

    synthesized = synthesize(f, rSet, pSet, atom2);
    expected = Polytope({{1,0},{0,1},{-1,0},{0,-1}},{0.55,0.15,-0.5,-0.1});

    BOOST_CHECK(synthesized.size()==1);
    
    // The following test does not work because of the approximation error
    // due to the double arithmetic
    // BOOST_CHECK(synthesized.begin()->simplify()==expected);
    BOOST_CHECK(epsilon_equivalent(*(synthesized.begin()), expected, APPROX_ERR));
}