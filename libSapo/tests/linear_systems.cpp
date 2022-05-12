#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE linear_system

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <sstream>

#include <glpk.h>

#include "LinearSystem.h"
#include "LinearAlgebra.h"
#include "LinearAlgebraIO.h"

template<typename T>
inline bool same_result_value(const OptimizationResult<T>& result, const T value)
{
    return (result.status()==GLP_OPT && result.optimum()==value);
}


template<typename T>
inline bool same_result_status(const OptimizationResult<T>& result, int status)
{
    return result.status()==status;
}

BOOST_AUTO_TEST_CASE(test_linear_systems)
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
 
    LinearSystem ls(A,b);

    std::vector<std::pair<std::pair<Vector<double>, bool>, double>> min_probs = {
        {{{1,0,0}, true}, 1},
        {{{0,1,0}, true}, 2},
        {{{0,0,1}, true}, 3},
        {{{25,0,0}, true}, 25},
        {{{-1,0,0}, true}, 3},
        {{{0,-1,0}, true}, 2},
        {{{0,0,-1}, true}, 1},
        {{{1,0,0}, false}, -3},
        {{{25,0,0}, false}, -75},
        {{{0,1,0}, false}, -2},
        {{{0,0,1}, false}, -1},
        {{{-1,0,0}, false}, -1},
        {{{0,-1,0}, false}, -2},
        {{{0,0,-1}, false}, -3}
    };

    for (auto &prob: min_probs) {
        auto input = prob.first;
        auto output = prob.second;
        auto result = ls.optimize(input.first, input.second);
        if (input.second) {
            auto result2 = ls.maximize(input.first);
            BOOST_CHECK(result.optimum()==result2.optimum() && result.status()==result2.status());
            BOOST_CHECK_MESSAGE(same_result_value(result, output), 
                                "maximizing " << input.first << " on " << ls 
                                << " produces " << result.optimum() << ": "
                                << output << " was expected.");
        } else {
            auto result2 = ls.minimize(input.first);
            BOOST_CHECK(result.optimum()==result2.optimum() && result.status()==result2.status());
            BOOST_CHECK_MESSAGE(same_result_value(result, output), 
                                "minimizing " << input.first << " on " << ls 
                                << " produces " << result.optimum() << ": "
                                << output << " was expected.");
        }
    }
}

BOOST_AUTO_TEST_CASE(test_unbounded_linear_systems)
{
    using namespace LinearAlgebra;
    using namespace LinearAlgebra::Dense;

    Matrix<double> A = {
        {0,1,0},
        {0,0,1},
        {-1,0,0},
        {0,-1,0}
    };

    Vector<double> b = {2,3,3,2};

    LinearSystem ls(A,b);

    Vector<double> obj{1,0,0};
    OptimizationResult<double> result = ls.maximize(obj);
    BOOST_CHECK_MESSAGE(same_result_status(result, GLP_UNBND),
                        "maximizing " << obj << " on " << ls 
                        << " produces " << result.optimum() 
                        << ": GLP_UNBND was expected.");

    obj = {0,0,1};
    result = ls.minimize(obj);
    BOOST_CHECK_MESSAGE(same_result_status(result, GLP_UNBND), 
                        "minimizing " << obj << " on " << ls 
                        << " produces " << result.optimum() << ": "
                        << "GLP_UNBND was expected.");
}

BOOST_AUTO_TEST_CASE(test_unfeasable_linear_systems)
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

    Vector<double> b = {1,2,3,-3,2,1};

    LinearSystem ls(A,b);

    Vector<double> obj{1,0,0};
    OptimizationResult<double> result = ls.maximize(obj);
    BOOST_CHECK_MESSAGE(same_result_status(result, GLP_NOFEAS),
                        "maximizing " << obj << " on " << ls 
                        << " produces " << result.status() 
                        << ": GLP_NOFEAS was expected.");

    obj = {0,0,1};
    result = ls.minimize(obj);
    BOOST_CHECK_MESSAGE(same_result_status(result, GLP_NOFEAS), 
                        "minimizing " << obj << " on " << ls 
                        << " produces " << result.status() << ": "
                        << "GLP_NOFEAS was expected.");
}

BOOST_AUTO_TEST_CASE(test_has_solutions_linear_systems)
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

    Vector<double> b = {1,2,3,-3,2,1};

    LinearSystem ls(A,b);

    // Test non-empty solution set
    BOOST_CHECK(!ls.has_solutions());

    b[3]=-1;
    ls = LinearSystem(A,b);

    // Test non-empty solution set
    BOOST_CHECK(ls.has_solutions());

    // Test non-empty interior
    BOOST_CHECK(!ls.has_solutions(true));

    b[3]=1;
    ls = LinearSystem(A,b);

    // Test non-empty interior
    BOOST_CHECK(ls.has_solutions());
}
