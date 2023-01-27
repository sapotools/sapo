#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE linear_system

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include "Simplex.h"
#include "LinearAlgebraIO.h"


BOOST_AUTO_TEST_CASE(test_simplex_optimum)
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

    std::vector<std::pair<std::pair<Vector<double>, OptimizationGoal>, 
                          double>> min_probs = {
        {{{1,0,0}, OptimizationGoal::MAXIMIZE}, 1},
        {{{0,1,0}, OptimizationGoal::MAXIMIZE}, 2},
        {{{0,0,1}, OptimizationGoal::MAXIMIZE}, 3},
        {{{25,0,0}, OptimizationGoal::MAXIMIZE}, 25},
        {{{-1,0,0}, OptimizationGoal::MAXIMIZE}, 3},
        {{{0,-1,0}, OptimizationGoal::MAXIMIZE}, 2},
        {{{0,0,-1}, OptimizationGoal::MAXIMIZE}, 1},
        {{{1,0,0}, OptimizationGoal::MINIMIZE}, -3},
        {{{25,0,0}, OptimizationGoal::MINIMIZE}, -75},
        {{{0,1,0}, OptimizationGoal::MINIMIZE}, -2},
        {{{0,0,1}, OptimizationGoal::MINIMIZE}, -1},
        {{{-1,0,0}, OptimizationGoal::MINIMIZE}, -1},
        {{{0,-1,0}, OptimizationGoal::MINIMIZE}, -2},
        {{{0,0,-1}, OptimizationGoal::MINIMIZE}, -3}
    };

    for (auto &prob: min_probs) {
        const auto& input = prob.first;
        const auto& output = prob.second;
        SimplexMethodOptimizer optimizer;
        auto result = optimizer(A, b, input.first, input.second);
        auto value = input.first*result.optimum();
        BOOST_CHECK_MESSAGE(output==value, 
                            "maximizing " << input.first << " on " << 
                            "A x<= b where A: " << 
                            A << " and b: " << b <<
                            " produces " << value << ": " <<
                            output << " was expected.");
    }
}

BOOST_AUTO_TEST_CASE(test_simplex_unbounded)
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

    SimplexMethodOptimizer optimizer;

    Vector<double> obj{1,0,0};
    auto result = optimizer(A, b, obj, OptimizationGoal::MAXIMIZE);
    BOOST_CHECK_MESSAGE(result.status()==OptimizationResult<double>::UNBOUNDED,
                        "maximizing " << obj << " on " << 
                        "A x<= b where A: " << 
                        A << " and b: " << b <<
                        " produces " << result.status() <<
                        ": UNBOUNDED was expected.");

    obj = {0,0,1};
    result = optimizer(A, b, obj, OptimizationGoal::MINIMIZE);
    BOOST_CHECK_MESSAGE(result.status()==OptimizationResult<double>::UNBOUNDED,
                        "minimizing " << obj << " on " << 
                        "A x<= b where A: " << 
                        A << " and b: " << b <<
                        " produces " << result.status() << 
                        ": UNBOUNDED was expected.");
}


BOOST_AUTO_TEST_CASE(test_simplex_unbounded2)
{
    using namespace LinearAlgebra;
    using namespace LinearAlgebra::Dense;

    Matrix<double> A = {
        {1,0},
        {0,1}
    };

    Vector<double> b = {2,3};

    SimplexMethodOptimizer optimizer;

    Vector<double> obj{1,1};
    auto result = optimizer(A, b, obj, OptimizationGoal::MAXIMIZE);
    BOOST_CHECK_MESSAGE(result.status()==OptimizationResult<double>::OPTIMUM_AVAILABLE,
                        "maximizing " << obj << " on " << 
                        "A x<= b where A: " << 
                        A << " and b: " << b <<
                        " produces " << result.status() <<
                        ": OPTIMUM_AVAILABLE was expected.");

    result = optimizer(A, b, obj, OptimizationGoal::MINIMIZE);
    BOOST_CHECK_MESSAGE(result.status()==OptimizationResult<double>::UNBOUNDED,
                        "minimizing " << obj << " on " << 
                        "A x<= b where A: " << 
                        A << " and b: " << b <<
                        " produces " << result.status() << 
                        ": UNBOUNDED was expected.");
}

BOOST_AUTO_TEST_CASE(test_simplex_unfeasible)
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

    Vector<double> obj{1,0,0};
    
    SimplexMethodOptimizer optimizer;

    auto result = optimizer(A,b,obj,OptimizationGoal::MAXIMIZE);
    BOOST_CHECK_MESSAGE(result.status()==OptimizationResult<double>::INFEASIBLE,
                        "maximizing " << obj << " on "
                        "A x<= b where A: " << 
                        A << " and b: " << b <<
                        " produces " << result.status() <<
                        ": INFEASIBLE was expected.");

    obj = {0,0,1};
    result = optimizer(A,b,obj,OptimizationGoal::MINIMIZE);
    BOOST_CHECK_MESSAGE(result.status()==OptimizationResult<double>::INFEASIBLE,
                        "minimizing " << obj << " on "
                        "A x<= b where A: " << 
                        A << " and b: " << b << 
                        " produces " << result.status() << 
                        ": INFEASIBLE was expected.");
}
