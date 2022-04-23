#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE linear_algebra

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include "LinearAlgebra.h"
#include "LinearAlgebraIO.h"

typedef boost::mpl::list<double> test_types;

BOOST_AUTO_TEST_CASE_TEMPLATE(test_norm_1, T, test_types)
{
    std::vector<std::pair<std::vector<T>, T>> tests{
        {{}, 0},
        {{1,7,5,0}, 13},
        {{-1,-7,-5,0}, 13},
        {{0,0,0,0}, 0},
        {{1,0,0,0}, 1},
        {{0,0,1,0}, 1}
    };

    for (auto test_it = std::begin(tests); test_it != std::end(tests); ++test_it) {
        T value = norm_1(test_it->first);
        bool beval = (value == test_it->second);
        BOOST_REQUIRE_MESSAGE(beval, "norm_1(" << test_it->first << ") == "
                                          << value  << " != " << test_it->second);   
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_norm_2, T, test_types)
{
    std::vector<std::pair<std::vector<T>, T>> tests{
        {{}, 0},
        {{1,7,5,5}, 10},
        {{-1,-7,-5,-5}, 10},
        {{-1,-7,-5,0}, std::sqrt(75)},
        {{0,0,0,0}, 0},
        {{1,0,0,0}, 1},
        {{0,0,-1,0}, 1}
    };

    for (auto test_it = std::begin(tests); test_it != std::end(tests); ++test_it) {
        T value = norm_2(test_it->first);
        bool beval = (value == test_it->second);
        BOOST_REQUIRE_MESSAGE(beval, "norm_2(" << test_it->first << ") == "
                                          << value  << " != " << test_it->second);  
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_norm_inf, T, test_types)
{
    std::vector<std::pair<std::vector<T>, T>> tests{
        {{}, 0},
        {{1,7,5,5}, 7},
        {{-1,7,-5,0}, 7},
        {{-1,-7,-5,0}, 7},
        {{0,0,0,0}, 0},
        {{1,0,0,0}, 1},
        {{0,0,-1,0}, 1}
    };

    for (auto test_it = std::begin(tests); test_it != std::end(tests); ++test_it) {
        T value = norm_infinity(test_it->first);
        bool beval = (value == test_it->second);
        BOOST_REQUIRE_MESSAGE(beval, "norm_inf(" << test_it->first << ") == "
                                          << value  << " != " << test_it->second);  
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_eq, T, test_types)
{
    std::vector<std::pair<std::pair<std::vector<T>, std::vector<T>>, bool>> tests{
        {{{},{}}, true},
        {{{1,7,5,5}, {1,7,5,5}}, true},
        {{{-1,-7,-5,-5}, {-1,-7,-5,-5}}, true},
        {{{0,3,-5,1}, {0,0,0,0}}, false},
        {{{0,0,0,0}, {0,3,-5,1}}, false},
        {{{1,7,5,5}, {-1,-7,-5,-5}}, false},
        {{{-1,-7,-5,-5}, {1,7,5,5}}, false},
        {{{0,-1,-7,-5,-5}, {0,1,7,5,5}}, false},
        {{{0,1,7,5,5}, {0,1,7,5,5}}, true},
        {{{-1,-7,-5,-5,0}, {0,1,7,5,5}}, false},
        {{{1,7,5,5,0}, {1,7,5,5,0}}, true},
        {{{1,7,5,5,0}, {1,7,5,-5,0}}, false},
    };

    for (auto test_it = std::begin(tests); test_it != std::end(tests); ++test_it) {
        auto eq = (test_it->first.first == test_it->first.second);
        bool beval = (eq == test_it->second);
        BOOST_REQUIRE_MESSAGE(beval, "(" << test_it->first.first << " == "
                                          << test_it->first.second  << ") == " 
                                          << eq << " != " << test_it->second);  
    }

    std::vector<T> v1{1,7,5,5}, v2{0,1,7,5,5};
    BOOST_REQUIRE_THROW(operator==(v1,v2), std::domain_error);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_neq, T, test_types)
{
    std::vector<std::pair<std::pair<std::vector<T>, std::vector<T>>, bool>> tests{
        {{{},{}}, false},
        {{{1,7,5,5}, {1,7,5,5}}, false},
        {{{-1,-7,-5,-5}, {-1,-7,-5,-5}}, false},
        {{{0,3,-5,1}, {0,0,0,0}}, true},
        {{{0,0,0,0}, {0,3,-5,1}}, true},
        {{{1,7,5,5}, {-1,-7,-5,-5}}, true},
        {{{-1,-7,-5,-5}, {1,7,5,5}}, true},
        {{{0,-1,-7,-5,-5}, {0,1,7,5,5}}, true},
        {{{0,1,7,5,5}, {0,1,7,5,5}}, false},
        {{{-1,-7,-5,-5,0}, {0,1,7,5,5}}, true},
        {{{1,7,5,5,0}, {1,7,5,5,0}}, false},
        {{{1,7,5,5,0}, {1,7,5,-5,0}}, true},
    };

    for (auto test_it = std::begin(tests); test_it != std::end(tests); ++test_it) {
        auto eq = (test_it->first.first != test_it->first.second);
        bool beval = (eq == test_it->second);
        BOOST_REQUIRE_MESSAGE(beval, "(" << test_it->first.first << " != "
                                          << test_it->first.second  << ") == " 
                                          << eq << " != " << test_it->second);  
    }

    std::vector<T> v1{1,7,5,5}, v2{0,1,7,5,5};
    BOOST_REQUIRE_THROW(operator!=(v1,v2), std::domain_error);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_subtraction, T, test_types)
{
    std::vector<std::pair<std::pair<std::vector<T>, std::vector<T>>, std::vector<T>>> tests{
        {{{},{}},{}},
        {{{1,7,5,5}, {1,7,5,5}}, {0,0,0,0}},
        {{{-1,-7,-5,-5}, {-1,-7,-5,-5}}, {0,0,0,0}},
        {{{0,3,-5,1}, {0,0,0,0}}, {0,3,-5,1}},
        {{{0,0,0,0}, {0,3,-5,1}}, {0,-3,5,-1}},
        {{{1,7,5,5}, {-1,-7,-5,-5}}, {2,14,10,10}},
        {{{-1,-7,-5,-5}, {1,7,5,5}}, {-2,-14,-10,-10}},
        {{{0,-1,-7,-5,-5}, {0,1,7,5,5}}, {0,-2,-14,-10,-10}},
        {{{0,1,7,5,5}, {0,1,7,5,5}}, {0,0,0,0,0}},
        {{{-1,7,-5,5,0}, {0,1,-7,5,5}}, {-1,6,2,0,-5}},
        {{{1,7,5,5,0}, {1,7,5,5,0}}, {0,0,0,0,0}},
        {{{1,7,5,5,0}, {1,7,5,-5,0}}, {0,0,0,10,0}},
    };

    for (auto test_it = std::begin(tests); test_it != std::end(tests); ++test_it) {
        auto difference = test_it->first.first - test_it->first.second;
        bool beval = (difference == test_it->second);
        BOOST_REQUIRE_MESSAGE(beval, "(" << test_it->first.first << " - "
                                          << test_it->first.second  << ") == " 
                                          << difference << " != " << test_it->second);   
    }

    std::vector<T> v1{1,7,5,5}, v2{0,1,7,5,5};
    BOOST_REQUIRE_THROW(v1-v2, std::domain_error);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_sum, T, test_types)
{
    std::vector<std::pair<std::pair<std::vector<T>, std::vector<T>>, std::vector<T>>> tests{
        {{{},{}},{}},
        {{{1,7,5,5}, {1,7,5,5}}, {2,14,10,10}},
        {{{-1,-7,-5,-5}, {-1,-7,-5,-5}}, {-2,-14,-10,-10}},
        {{{0,3,-5,1}, {0,0,0,0}}, {0,3,-5,1}},
        {{{0,0,0,0}, {0,3,-5,1}}, {0,3,-5,1}},
        {{{1,7,5,5}, {-1,-7,-5,-5}}, {0,0,0,0}},
        {{{-1,-7,-5,-5}, {1,7,5,5}}, {0,0,0,0}},
        {{{0,-1,-7,-5,-5}, {0,1,7,5,5}}, {0,0,0,0,0}},
        {{{0,1,7,5,5}, {0,1,7,5,5}}, {0,2,14,10,10}},
        {{{-1,7,-5,5,0}, {0,1,-7,5,5}}, {-1,8,-12,10,5}},
        {{{1,7,5,5,0}, {1,7,5,5,0}}, {2,14,10,10,0}},
        {{{1,7,5,5,0}, {1,7,5,-5,0}}, {2,14,10,0,0}},
    };

    for (auto test_it = std::begin(tests); test_it != std::end(tests); ++test_it) {
        auto sum = test_it->first.first + test_it->first.second;
        bool beval = (sum == test_it->second);
        BOOST_REQUIRE_MESSAGE(beval, "(" << test_it->first.first << " + "
                                          << test_it->first.second  << ") == " 
                                          << sum << " != " << test_it->second);   
    }

    std::vector<T> v1{1,7,5,5}, v2{0,1,7,5,5};
    BOOST_REQUIRE_THROW(v1+v2, std::domain_error);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_mult_pre, T, test_types)
{
    std::vector<std::pair<std::pair<T, std::vector<T>>, std::vector<T>>> tests{
        {{0,{}},{}},
        {{8,{}},{}},
        {{-4,{}},{}},
        {{1,{}},{}},
        {{0, {1,0,-5,7}}, {0,0,0,0}},
        {{8, {1,0,-5,7}}, {8,0,-40,56}},
        {{-4, {1,0,-5,7}}, {-4,0,20,-28}},
        {{1, {1,0,-5,7}}, {1,0,-5,7}}
    };

    for (auto test_it = std::begin(tests); test_it != std::end(tests); ++test_it) {
        auto mult = test_it->first.first * test_it->first.second;
        bool beval = (mult == test_it->second);
        BOOST_REQUIRE_MESSAGE(beval, "(" << test_it->first.first << " * "
                                          << test_it->first.second  << ") == " 
                                          << mult << " != " << test_it->second);   
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_mult_post, T, test_types)
{
    std::vector<std::pair<std::pair<std::vector<T>, T>, std::vector<T>>> tests{
        {{{},0},{}},
        {{{},8},{}},
        {{{},-4},{}},
        {{{},1},{}},
        {{{1,0,-5,7},0}, {0,0,0,0}},
        {{{1,0,-5,7},8}, {8,0,-40,56}},
        {{{1,0,-5,7},-4}, {-4,0,20,-28}},
        {{{1,0,-5,7},1}, {1,0,-5,7}}
    };

    for (auto test_it = std::begin(tests); test_it != std::end(tests); ++test_it) {
        auto mult = test_it->first.first * test_it->first.second;
        bool beval = (mult == test_it->second);
        BOOST_REQUIRE_MESSAGE(beval, "(" << test_it->first.first << " * "
                                          << test_it->first.second  << ") == " 
                                          << mult << " != " << test_it->second);   
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_div, T, test_types)
{
    std::vector<std::pair<std::pair<std::vector<T>, T>, std::vector<T>>> tests{
        {{{},8},{}},
        {{{},-4},{}},
        {{{},1},{}},
        {{{8,0,-40,56},8}, {1,0,-5,7}},
        {{{-4,0,20,-28},-4}, {1,0,-5,7}},
        {{{1,0,-5,7},1}, {1,0,-5,7}}
    };

    for (auto test_it = std::begin(tests); test_it != std::end(tests); ++test_it) {
        auto div = test_it->first.first / test_it->first.second;
        bool beval = (div == test_it->second);
        BOOST_REQUIRE_MESSAGE(beval, "(" << test_it->first.first << " / "
                                          << test_it->first.second  << ") == " 
                                          << div << " != " << test_it->second);   
    }

    std::vector<std::pair<std::vector<T>, T>> err_test{
        {{},0},
        {{1,0,-5,7},0},
    };
    for (auto test_it = std::begin(err_test); test_it != std::end(err_test); ++test_it) {
        BOOST_REQUIRE_THROW(operator/(test_it->first,test_it->second), std::domain_error); 
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_neg, T, test_types)
{
    std::vector<std::pair<std::vector<T>, std::vector<T>>> tests{
        {{},{}},
        {{1,7,5,5}, {-1,-7,-5,-5}},
        {{-1,-7,-5,-5}, {1,7,5,5}},
        {{-1,-7,0,-5}, {1,7,0,5}},
        {{0,0,0,0}, {0,0,0,0}},
    };

    for (auto test_it = std::begin(tests); test_it != std::end(tests); ++test_it) {
        auto neg = -test_it->first;
        bool beval = (neg == test_it->second);
        BOOST_REQUIRE_MESSAGE(beval, "-(" << test_it->first << ") == " 
                                          << neg << " != " << test_it->second);    
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_prod, T, test_types)
{
    std::vector<std::pair<std::pair<std::vector<T>, std::vector<T>>, T>> tests{
        {{{},{}},0},
        {{{1,7,5,5}, {1,7,5,5}}, 100},
        {{{-1,-7,-5,-5}, {-1,-7,-5,-5}}, 100},
        {{{0,3,-5,1}, {0,0,0,0}}, 0},
        {{{0,0,0,0}, {0,3,-5,1}}, 0},
        {{{1,7,5,5}, {-1,-7,-5,-5}}, -100},
        {{{-1,-7,-5,-5}, {1,7,5,5}}, -100},
        {{{0,-1,-7,-5,-5}, {0,1,7,5,5}}, -100},
        {{{0,1,7,5,5}, {0,1,7,5,5}}, 100},
        {{{-1,7,-5,5,0}, {0,1,-7,5,5}}, 67},
        {{{1,7,5,5,0}, {1,7,5,5,0}}, 100},
        {{{1,7,5,5,0}, {1,7,5,-5,0}}, 50},
    };

    for (auto test_it = std::begin(tests); test_it != std::end(tests); ++test_it) {
        auto prod = test_it->first.first * test_it->first.second;
        bool beval = (prod == test_it->second);
        BOOST_REQUIRE_MESSAGE(beval, "(" << test_it->first.first << " * "
                                          << test_it->first.second  << ") == " 
                                          << prod << " != " << test_it->second);   
    }

    std::vector<T> v1{1,7,5,5}, v2{0,1,7,5,5};
    BOOST_REQUIRE_THROW(v1*v2, std::domain_error);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_H_prod, T, test_types)
{
    std::vector<std::pair<std::pair<std::vector<T>, std::vector<T>>, std::vector<T>>> tests{
        {{{},{}},{}},
        {{{1,7,5,5}, {1,7,5,5}}, {1,49,25,25}},
        {{{-1,-7,-5,-5}, {-1,-7,-5,-5}}, {1,49,25,25}},
        {{{0,3,-5,1}, {0,0,0,0}}, {0,0,0,0}},
        {{{0,0,0,0}, {0,3,-5,1}}, {0,0,0,0}},
        {{{1,7,5,5}, {-1,-7,-5,-5}}, {-1,-49,-25,-25}},
        {{{-1,-7,-5,-5}, {1,7,5,5}}, {-1,-49,-25,-25}},
        {{{0,-1,-7,-5,-5}, {0,1,7,5,5}}, {0,-1,-49,-25,-25}},
        {{{0,1,7,5,5}, {0,1,7,5,5}}, {0,1,49,25,25}},
        {{{-1,7,-5,5,0}, {0,1,-7,5,5}}, {0,7,35,25,0}},
        {{{1,7,5,5,0}, {1,7,5,5,0}}, {1,49,25,25,0}},
        {{{1,7,5,5,0}, {1,7,5,-5,0}}, {1,49,25,-25,0}},
    };

    for (auto test_it = std::begin(tests); test_it != std::end(tests); ++test_it) {
        auto prod = H_prod(test_it->first.first, test_it->first.second);
        bool beval = (prod == test_it->second);
        BOOST_REQUIRE_MESSAGE(beval, "H_prod(" << test_it->first.first << ","
                                          << test_it->first.second  << ") == " 
                                          << prod << " != " << test_it->second);   
    }

    std::vector<T> v1{1,7,5,5}, v2{0,1,7,5,5};
    BOOST_REQUIRE_THROW(H_prod(v1,v2), std::domain_error);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_linear_dep, T, test_types)
{
    std::vector<std::pair<std::pair<std::vector<T>, std::vector<T>>, bool>> tests{
        {{{1,7,5,5}, {1,7,5,5}}, true},
        {{{2,14,10,10}, {1,7,5,5}}, true},
        {{{1,7,5,5}, {2,14,10,10}}, true},
        {{{-1,-7,-5,-5}, {-1,-7,-5,-5}}, true},
        {{{0,3,-5,1}, {0,0,0,0}}, false},
        {{{0,0,0,0}, {0,3,-5,1}}, false},
        {{{1,7,5,5}, {-1,-7,-5,-5}}, true},
        {{{-1,-7,-5,-5}, {1,7,5,5}}, true},
        {{{0,-1,-7,-5,-5}, {0,1,7,5,5}}, true},
        {{{0,1,7,5,5}, {0,1,49,25,25}}, false},
        {{{-1,7,-5,5,0}, {0,1,-7,5,5}}, false},
        {{{1,7,5,5,0}, {1,7,5,5,0}}, true},
        {{{1,7,5,5,0}, {1,7,5,-5,0}}, false},
    };

    for (auto test_it = std::begin(tests); test_it != std::end(tests); ++test_it) {
        auto ldep = are_linearly_dependent(test_it->first.first, test_it->first.second);
        bool beval = (ldep == test_it->second);
        BOOST_REQUIRE_MESSAGE(beval, "are_linearly_dependent(" << test_it->first.first << ","
                                          << test_it->first.second  << ") == " 
                                          << ldep << " != " << test_it->second);   
    }

    std::vector<std::pair<std::vector<T>, std::vector<T>>> err_test{
        {{},{}},
        {{1,0,-5,7},{0,1,7,5,5}},
    };
    for (auto test_it = std::begin(err_test); test_it != std::end(err_test); ++test_it) {
        BOOST_REQUIRE_THROW(are_linearly_dependent(test_it->first,test_it->second), std::domain_error); 
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_linear_dep_ratio, T, test_types)
{
    std::vector<std::pair<std::pair<std::vector<T>, std::vector<T>>, T>> tests{
        {{{1,7,5,5}, {1,7,5,5}}, 1},
        {{{2,14,10,10}, {1,7,5,5}}, 2},
        {{{1,7,5,5}, {2,14,10,10}}, ((T)1)/2},
        {{{-1,-7,-5,-5}, {-1,-7,-5,-5}}, 1},
        {{{1,7,5,5}, {-1,-7,-5,-5}}, -1},
        {{{-1,-7,-5,-5}, {1,7,5,5}}, -1},
        {{{0,-1,-7,-5,-5}, {0,1,7,5,5}}, -1},
        {{{1,7,5,5,0}, {1,7,5,5,0}}, 1},
    };

    for (auto test_it = std::begin(tests); test_it != std::end(tests); ++test_it) {
        auto ratio = test_it->first.first / test_it->first.second;
        bool beval = (ratio == test_it->second);
        BOOST_REQUIRE_MESSAGE(beval, "(" << test_it->first.first << " / "
                                          << test_it->first.second  << ") == " 
                                          << ratio << " != " << test_it->second);   
    }

    std::vector<std::pair<std::vector<T>, std::vector<T>>> err_test{
        {{},{}},
        {{1,0,-5,7},{0,1,7,5,5}},
        {{0,3,-5,1}, {0,0,0,0}},
        {{0,0,0,0}, {0,3,-5,1}},
        {{0,1,7,5,5}, {0,1,49,25,25}},
        {{-1,7,-5,5,0}, {0,1,-7,5,5}},
        {{1,7,5,5,0}, {1,7,5,-5,0}}
    };
    for (auto test_it = std::begin(err_test); test_it != std::end(err_test); ++test_it) {
        BOOST_REQUIRE_THROW(test_it->first / test_it->second, std::domain_error); 
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_dense_matrix_eq, T, test_types)
{
    using namespace DenseLinearAlgebra;
    std::vector<std::pair<std::pair<Matrix<T>, Matrix<T>>, bool>> tests{
        {{{{}},{{}}}, true},
        {{{{0,1,2},{3,4,5},{6,7,8}},{{0,1,2},{3,4,5},{6,7,8}}}, true},
        {{{{0,1,2},{3,4,5},{6,7,8}},{{0,1,2},{3,-4,5},{6,7,8}}}, false},
        {{{{0,1,2},{3,4,5},{6,7,8}},{{0,3,6},{1,4,7},{2,5,8}}}, false},
        {{{{0,1,2},{3,4,5},{6,7,8}},{{0,1,2},{3,4,5},{0,7,8}}}, false}
    };

    for (auto test_it = std::begin(tests); test_it != std::end(tests); ++test_it) {
        auto meq = test_it->first.first == test_it->first.second;
        bool beval = (meq == test_it->second);
        BOOST_REQUIRE_MESSAGE(beval, "(" << test_it->first.first << " == "
                                          << test_it->first.second  << ") == " 
                                          << meq << " != " << test_it->second);   
    }

    std::vector<std::pair<Matrix<T>, Matrix<T>>> err_test{
        {{{0,1,2},{3,4,5},{6,7,8}},
         {{0,0,1,2},{1,3,4,5},{2,6,7,8},{3,0,7,1}}},
    };
    for (auto test_it = std::begin(err_test); test_it != std::end(err_test); ++test_it) {
        BOOST_REQUIRE_THROW(operator==(test_it->first,test_it->second), std::domain_error); 
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_dense_matrix_neq, T, test_types)
{
    using namespace DenseLinearAlgebra;
    std::vector<std::pair<std::pair<Matrix<T>, Matrix<T>>, bool>> tests{
        {{{{}},{{}}}, false},
        {{{{0,1,2},{3,4,5},{6,7,8}},{{0,1,2},{3,4,5},{6,7,8}}}, false},
        {{{{0,1,2},{3,4,5},{6,7,8}},{{0,1,2},{3,-4,5},{6,7,8}}}, true},
        {{{{0,1,2},{3,4,5},{6,7,8}},{{0,3,6},{1,4,7},{2,5,8}}}, true},
        {{{{0,1,2},{3,4,5},{6,7,8}},{{0,1,2},{3,4,5},{0,7,8}}}, true}
    };

    for (auto test_it = std::begin(tests); test_it != std::end(tests); ++test_it) {
        auto mneq = test_it->first.first != test_it->first.second;
        bool beval = (mneq == test_it->second);
        BOOST_REQUIRE_MESSAGE(beval, "(" << test_it->first.first << " != "
                                          << test_it->first.second  << ") == " 
                                          << mneq << " != " << test_it->second);   
    }

    std::vector<std::pair<Matrix<T>, Matrix<T>>> err_test{
        {{{0,1,2},{3,4,5},{6,7,8}},
         {{0,0,1,2},{1,3,4,5},{2,6,7,8},{3,0,7,1}}},
    };
    for (auto test_it = std::begin(err_test); test_it != std::end(err_test); ++test_it) {
        BOOST_REQUIRE_THROW(operator!=(test_it->first,test_it->second), std::domain_error); 
    }
}


BOOST_AUTO_TEST_CASE_TEMPLATE(test_dense_matrix_sum, T, test_types)
{
    using namespace DenseLinearAlgebra;
    std::vector<std::pair<std::pair<Matrix<T>, Matrix<T>>, Matrix<T>>> tests{
        {{{{}},{{}}}, {{}}},
        {{{{0,1,2},{3,4,5},{6,7,8}},{{0,1,2},{3,4,5},{6,7,8}}}, {{0,2,4},{6,8,10},{12,14,16}}},
        {{{{0,1,2},{3,4,5},{6,7,8}},{{0,1,2},{3,-4,5},{6,7,8}}}, {{0,2,4},{6,0,10},{12,14,16}}},
        {{{{0,1,2},{3,4,5},{6,7,8}},{{0,3,6},{1,4,7},{2,5,8}}}, {{0,4,8},{4,8,12},{8,12,16}}},
        {{{{0,1,2},{3,4,5},{6,7,8}},{{0,1,2},{3,4,5},{0,7,8}}}, {{0,2,4},{6,8,10},{6,14,16}}}
    };

    for (auto test_it = std::begin(tests); test_it != std::end(tests); ++test_it) {
        auto sum = test_it->first.first + test_it->first.second;
        bool beval = (sum == test_it->second);
        BOOST_REQUIRE_MESSAGE(beval, "(" << test_it->first.first << " + "
                                          << test_it->first.second  << ") == " 
                                          << sum << " != " << test_it->second);   
    }

    std::vector<std::pair<Matrix<T>, Matrix<T>>> err_test{
        {{{0,1,2},{3,4,5},{6,7,8}},
         {{0,0,1,2},{1,3,4,5},{2,6,7,8},{3,0,7,1}}},
    };
    for (auto test_it = std::begin(err_test); test_it != std::end(err_test); ++test_it) {
        BOOST_REQUIRE_THROW(operator+(test_it->first,test_it->second), std::domain_error); 
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_dense_transpose, T, test_types)
{
    using namespace DenseLinearAlgebra;
    std::vector<std::pair<Matrix<T>, Matrix<T>>> tests{
        {{},{}},
        {{{0,1,2},{3,4,5},{6,7,8}},{{0,3,6},{1,4,7},{2,5,8}}}
    };

    for (auto test_it = std::begin(tests); test_it != std::end(tests); ++test_it) {
        auto transp = transpose(test_it->first);
        bool beval = (transp == test_it->second);
        BOOST_REQUIRE_MESSAGE(beval, "transpose(" << test_it->first << ") == "
                                          << transp << " != " << test_it->second);   
    }
}