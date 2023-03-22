#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE simplex

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include "Simplex.h"
#include "LinearAlgebraIO.h"
#include "Approximation.h"

#ifdef HAVE_GMP
#include <gmpxx.h>

typedef boost::mpl::list<double, mpq_class> test_types;

template<>
struct is_punctual<mpq_class> : std::integral_constant<bool, true> {};

#else
typedef boost::mpl::list<double> test_types;
#endif

template<typename T>
inline T admitted_error();

template<>
inline double admitted_error<double>()
{
    return 1e-14;
}

template<>
inline mpq_class admitted_error<mpq_class>()
{
    return 0;
}

template<typename T>
inline bool approximate(const OptimizationResult<T>& result, const T exact)
{
    return (result.status()==result.OPTIMUM_AVAILABLE &&
            abs(result.objective_value()-exact) <= admitted_error<T>());
}

template<typename T>
inline bool is_exactly(const OptimizationResult<T>& result, const T exact)
{
    return (result.status()==result.OPTIMUM_AVAILABLE &&
            result.objective_value()==exact);
}


BOOST_AUTO_TEST_CASE_TEMPLATE(test_simplex_optimum, T, test_types)
{
    using namespace LinearAlgebra;
    using namespace LinearAlgebra::Dense;

    Matrix<T> A = {
        {1,0,0},
        {0,1,0},
        {0,0,1},
        {-1,0,0},
        {0,-1,0},
        {0,0,-1}
    };

    Vector<T> b = {1,2,3,3,2,1};

    std::vector<std::pair<std::pair<Vector<T>, OptimizationGoal>, 
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
        SimplexMethodOptimizer<T> optimizer;
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

BOOST_AUTO_TEST_CASE_TEMPLATE(test_simplex_unbounded, T, test_types)
{
    using namespace LinearAlgebra;
    using namespace LinearAlgebra::Dense;

    Matrix<T> A = {
        {0,1,0},
        {0,0,1},
        {-1,0,0},
        {0,-1,0}
    };

    Vector<T> b = {2,3,3,2};

    SimplexMethodOptimizer<T> optimizer;

    Vector<T> obj{1,0,0};
    auto result = optimizer(A, b, obj, OptimizationGoal::MAXIMIZE);
    BOOST_CHECK_MESSAGE(result.status()==OptimizationResult<T>::UNBOUNDED,
                        "maximizing " << obj << " on " << 
                        "A x<= b where A: " << 
                        A << " and b: " << b <<
                        " produces " << result.status() <<
                        ": UNBOUNDED was expected.");

    obj = {0,0,1};
    result = optimizer(A, b, obj, OptimizationGoal::MINIMIZE);
    BOOST_CHECK_MESSAGE(result.status()==OptimizationResult<T>::UNBOUNDED,
                        "minimizing " << obj << " on " << 
                        "A x<= b where A: " << 
                        A << " and b: " << b <<
                        " produces " << result.status() << 
                        ": UNBOUNDED was expected.");
}


BOOST_AUTO_TEST_CASE_TEMPLATE(test_simplex_unbounded2, T, test_types)
{
    using namespace LinearAlgebra;
    using namespace LinearAlgebra::Dense;

    Matrix<T> A = {
        {1,0},
        {0,1}
    };

    Vector<T> b = {2,3};

    SimplexMethodOptimizer<T> optimizer;

    Vector<T> obj{1,1};
    auto result = optimizer(A, b, obj, OptimizationGoal::MAXIMIZE);
    BOOST_CHECK_MESSAGE(result.status()==OptimizationResult<T>::OPTIMUM_AVAILABLE,
                        "maximizing " << obj << " on " << 
                        "A x<= b where A: " << 
                        A << " and b: " << b <<
                        " produces " << result.status() <<
                        ": OPTIMUM_AVAILABLE was expected.");

    result = optimizer(A, b, obj, OptimizationGoal::MINIMIZE);
    BOOST_CHECK_MESSAGE(result.status()==OptimizationResult<T>::UNBOUNDED,
                        "minimizing " << obj << " on " << 
                        "A x<= b where A: " << 
                        A << " and b: " << b <<
                        " produces " << result.status() << 
                        ": UNBOUNDED was expected.");
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_simplex_unfeasible, T, test_types)
{
    using namespace LinearAlgebra;
    using namespace LinearAlgebra::Dense;

    Matrix<T> A = {
        {1,0,0},
        {0,1,0},
        {0,0,1},
        {-1,0,0},
        {0,-1,0},
        {0,0,-1}
    };

    Vector<T> b = {1,2,3,-3,2,1};

    Vector<T> obj{1,0,0};
    
    SimplexMethodOptimizer<T> optimizer;

    auto result = optimizer(A,b,obj,OptimizationGoal::MAXIMIZE);
    BOOST_CHECK_MESSAGE(result.status()==OptimizationResult<T>::INFEASIBLE,
                        "maximizing " << obj << " on "
                        "A x<= b where A: " << 
                        A << " and b: " << b <<
                        " produces " << result.status() <<
                        ": INFEASIBLE was expected.");

    obj = {0,0,1};
    result = optimizer(A,b,obj,OptimizationGoal::MINIMIZE);
    BOOST_CHECK_MESSAGE(result.status()==OptimizationResult<T>::INFEASIBLE,
                        "minimizing " << obj << " on "
                        "A x<= b where A: " << 
                        A << " and b: " << b << 
                        " produces " << result.status() << 
                        ": INFEASIBLE was expected.");
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_simplex_3D_square, T, test_types)
{
    using namespace LinearAlgebra;
    using namespace LinearAlgebra::Dense;

    Matrix<T> A = {
        {1,0,0},
        {0,1,0},
        {0,0,1},
        {-1,0,0},
        {0,-1,0},
        {0,0,-1},
    };

    Vector<T> b = {1,2,3,-1,2,1};
    /*
        x in [1,1]
        y in [-2, 2]
        z in [-1, 3]
    */

    Vector<T> obj{1,0,0};
    
    SimplexMethodOptimizer<T> optimizer;
    
    auto result = optimizer(A,b,obj,OptimizationGoal::MAXIMIZE);
    BOOST_CHECK(is_exactly<T>(result,1));

    result = optimizer(A,b,obj,OptimizationGoal::MINIMIZE);
    BOOST_CHECK(is_exactly<T>(result,1));

    obj = {1,-1, 0};

    result = optimizer(A,b,obj,OptimizationGoal::MAXIMIZE);
    BOOST_CHECK(is_exactly<T>(result,3));

    result = optimizer(A,b,obj,OptimizationGoal::MINIMIZE);
    BOOST_CHECK(is_exactly<T>(result,-1));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_simplex_0_constraints, T, test_types)
{
    using namespace LinearAlgebra;
    using namespace LinearAlgebra::Dense;

    Matrix<T> A = {
        {1,0,0},
        {0,1,0},
        {0,0,1},
        {-1,0,0},
        {0,-1,0},
        {0,0,-1},
    };

    Vector<T> b = {10,7,16,0,0,0};
    b /= T(10);
    /*
        x in [0,1]
        y in [0,0.7]
        z in [0,1.6]
    */

    Vector<T> obj{1,0,0};
    
    SimplexMethodOptimizer<T> optimizer;
    
    auto result = optimizer(A,b,obj,OptimizationGoal::MAXIMIZE);
    BOOST_CHECK(is_exactly<T>(result,1));

    result = optimizer(A,b,obj,OptimizationGoal::MINIMIZE);
    BOOST_CHECK(is_exactly<T>(result,0));

    obj = {1,-1,0};

    result = optimizer(A,b,obj,OptimizationGoal::MAXIMIZE);
    BOOST_CHECK(is_exactly<T>(result,1));

    result = optimizer(A,b,obj,OptimizationGoal::MINIMIZE);
    BOOST_CHECK(is_exactly<T>(result,T(-7)/10));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_simplex_negative_constraints, T, test_types)
{
    using namespace LinearAlgebra;
    using namespace LinearAlgebra::Dense;

    Matrix<T> A = {
        {1,0},
        {0,1},
        {-1,0},
        {0,-1},
    };

    Vector<T> b = {995,1005,5,-5};
    b /= T(100);
    /*
        x in [-0.05, 9.95]
        y in [0.05, 10.05]
    */

    SimplexMethodOptimizer<T> optimizer;

    Vector<T> obj{1,-1};
    
    auto result = optimizer(A,b,obj,OptimizationGoal::MAXIMIZE);

    BOOST_CHECK(approximate<T>(result,T(99)/10));

    result = optimizer(A,b,obj,OptimizationGoal::MINIMIZE);
    BOOST_CHECK(approximate<T>(result,T(-101)/10));

    obj = {1,0};
    
    result = optimizer(A,b,obj,OptimizationGoal::MAXIMIZE);
    BOOST_CHECK(approximate<T>(result,T(995)/100));

    result = optimizer(A,b,obj,OptimizationGoal::MINIMIZE);
    BOOST_CHECK(approximate<T>(result,T(-5)/100));
}

BOOST_AUTO_TEST_CASE(test_simplex_GPLK_bug)
{
    using namespace LinearAlgebra;
    using namespace LinearAlgebra::Dense;

    Matrix<double> A = {
        {1,2.51581e-08},
        {-2.51581e-08,1},
        {-1,-2.51581e-08},
        {2.51581e-08,-1}
    };

    Vector<double> b = {9.95,10.05,0.05,-0.05};
    /*
       1*x + 2.51581e-08*y <= 9.95
       -2.51581e-08*x + 1*y <= 10.05
       -1*x + -2.51581e-08*y <= 0.05
       2.51581e-08*x + -1*y <= -0.05
    */

    SimplexMethodOptimizer<double> optimizer;

    Vector<double> obj{1,0};

    auto result = optimizer(A,b,obj,OptimizationGoal::MINIMIZE);
    auto objective_value = -0.05-2.52838e-07-9.04964e-13;

    for (const auto& value: b-A*result.optimum()) {
        BOOST_CHECK(value>=0);
    }
    BOOST_CHECK(approximate<double>(result,objective_value));
}

BOOST_AUTO_TEST_CASE(test_approximated_simplex)
{
    using namespace LinearAlgebra;
    using namespace LinearAlgebra::Dense;

    Matrix<double> A = {
        {1,2.51581e-08},
        {-2.51581e-08,1},
        {-1,-2.51581e-08},
        {2.51581e-08,-1}
    };

    Vector<double> b = {9.95,10.05,0.05,-0.05};
    /*
       1*x + 2.51581e-08*y <= 9.95
       -2.51581e-08*x + 1*y <= 10.05
       -1*x + -2.51581e-08*y <= 0.05
       2.51581e-08*x + -1*y <= -0.05
    */

    SimplexMethodOptimizer<double, Approximation<double>> optimizer;

    Vector<double> obj{1,0};

    auto result = optimizer(A,b,obj,OptimizationGoal::MINIMIZE);

    for (const auto& value: b-A*result.optimum()) {
        BOOST_CHECK(!((value<0).is_true()));
    }
}


#ifdef HAVE_GMP


BOOST_AUTO_TEST_CASE(test_Q_approximated_simplex_infeasible_despite_errors)
{
    using namespace LinearAlgebra;
    using namespace LinearAlgebra::Dense;

    mpq_class large_value{mpq_class(mpq_class(1) << (IEEE754Rounding<double>::mantissa_size+1))};

    Matrix<mpq_class> A = {
        {1,1},
        {(mpq_class(1))/large_value,1},
        {-1,0},
        {0,-1}
    };

    Vector<mpq_class> b = {2,1,-1,-1};

    Vector<mpq_class> obj{1,0};
    
    SimplexMethodOptimizer<mpq_class, Approximation<mpq_class>> optimizer;

    auto result = optimizer(A,b,obj,OptimizationGoal::MAXIMIZE);
    BOOST_CHECK_MESSAGE(result.status()==result.INFEASIBLE,
                        "maximizing " << obj << " on "
                        "A x<= b where A: " << 
                        A << " and b: " << b <<
                        " produces " << result.status() <<
                        ": INFEASIBLE was expected.");

    obj = {0,1};
    result = optimizer(A,b,obj,OptimizationGoal::MINIMIZE);
    BOOST_CHECK_MESSAGE(result.status()==result.INFEASIBLE,
                        "minimizing " << obj << " on "
                        "A x<= b where A: " << 
                        A << " and b: " << b << 
                        " produces " << result.status() << 
                        ": INFEASIBLE was expected.");
}
#endif