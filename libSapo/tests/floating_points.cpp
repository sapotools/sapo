#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE floating_points

#include <vector>
#include <functional>

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include "FloatingPoints.h"

typedef boost::mpl::list<float, double> test_types;

template<typename T>
inline T apply(const std::function<T (const T&, const T&, int)> f, const T& init_value, 
               const T& value, int rounding, const size_t num_of_terms)
{
    T a{init_value};
    for (size_t i=0; i<num_of_terms; ++i) {
        a = f(a,value,rounding);
    }

    return a;
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_sum_downward, T, test_types)
{
    using mantissa_type = typename IEEE754Rounding<T>::mantissa_type;

    T large_value{T(mantissa_type(1) << (IEEE754Rounding<T>::mantissa_size))};

    const std::vector<T> larges = {large_value, 1, -large_value, -1, 0, 
                                   T(1)/3000, -T(1)/3000};
    const std::vector<T> smalls = {1, T(1)/large_value, -1, -T(1)/large_value, 
                                   0, T(1)/3000, -T(1)/3000};

    constexpr T medium{1<<7};
    for (const auto& large : larges) {
        for (const auto& small : smalls) {
            T a{large};
            for (size_t i=0; i<medium; ++i) {
                a = add(a,small,FE_DOWNWARD);
            }
            BOOST_CHECK_MESSAGE(a<=large+medium*small, "sum("<<large << ","
                                                       << small<<"," << medium 
                                                       << ", FE_DOWNWARD)="
                                                       << a << " > " 
                                                       << large << "+"<< medium
                                                       <<"*"<< small);
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_sum_upward, T, test_types)
{
    using mantissa_type = typename IEEE754Rounding<T>::mantissa_type;

    T large_value{T(mantissa_type(1) << (IEEE754Rounding<T>::mantissa_size))};

    const std::vector<T> larges = {large_value, 1, -large_value, -1, 0, 
                                   T(1)/3000, -T(1)/3000};
    const std::vector<T> smalls = {1, T(1)/large_value, -1, -T(1)/large_value, 
                                   0, T(1)/3000, -T(1)/3000};

    constexpr T medium{1<<7};
    for (const auto& large : larges) {
        for (const auto& small : smalls) {
            T a{large};
            for (size_t i=0; i<medium; ++i) {
                a = add(a,small,FE_UPWARD);
            }
            BOOST_CHECK_MESSAGE(a>=large+medium*small, "sum("<<large << ","
                                                       << small<<"," << medium 
                                                       << ", FE_UPWARD)="
                                                       << a << " < " 
                                                       << large << "+"<< medium
                                                       <<"*"<< small);
        }
    }
}


BOOST_AUTO_TEST_CASE_TEMPLATE(test_difference_downward, T, test_types)
{
    using mantissa_type = typename IEEE754Rounding<T>::mantissa_type;

    T large_value{T(mantissa_type(1) << (IEEE754Rounding<T>::mantissa_size))};

    const std::vector<T> larges = {large_value, 1, -large_value, -1, 0, 
                                   T(1)/3000, -T(1)/3000};
    const std::vector<T> smalls = {1, T(1)/large_value, -1, -T(1)/large_value, 
                                   0, T(1)/3000, -T(1)/3000};

    constexpr T medium{1<<7};
    for (const auto& large : larges) {
        for (const auto& small : smalls) {
            T a{large};
            for (size_t i=0; i<medium; ++i) {
                a = subtract(a,small,FE_DOWNWARD);
            }
            BOOST_CHECK_MESSAGE(a<=large-medium*small, "difference("<<large << ","
                                                       << small<<"," << medium 
                                                       << ", FE_DOWNWARD)="
                                                       << a << " > " << large << "-"
                                                       << medium<<"*"<< small << "=" 
                                                       << large-medium*small);
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_difference_upward, T, test_types)
{
    using mantissa_type = typename IEEE754Rounding<T>::mantissa_type;

    T large_value{T(mantissa_type(1) << (IEEE754Rounding<T>::mantissa_size))};

    const std::vector<T> larges = {large_value, 1, -large_value, -1, 0, 
                                   T(1)/3000, -T(1)/3000};
    const std::vector<T> smalls = {1, T(1)/large_value, -1, -T(1)/large_value, 
                                   0, T(1)/3000, -T(1)/3000};

    constexpr T medium{1<<7};
    for (const auto& large : larges) {
        for (const auto& small : smalls) {
            T a{large};
            for (size_t i=0; i<medium; ++i) {
                a = subtract(a,small,FE_UPWARD);
            }
            BOOST_CHECK_MESSAGE(a>=large-medium*small, "difference("<<large << ","<< small
                                                        <<"," << medium 
                                                        << ", FE_UPWARD)=" << a << " < " 
                                                        << large << "-" << medium<<"*"
                                                        << small << "=" << large-medium*small);
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_product_downward, T, test_types)
{
    using mantissa_type = typename IEEE754Rounding<T>::mantissa_type;

    T large_value{T(mantissa_type(1) << (IEEE754Rounding<T>::mantissa_size))};

    const std::vector<T> larges = {large_value, 1, -large_value, -1, 7, -7, 0};
    const std::vector<T> smalls = {1, T(1)/large_value, -1, -T(1)/large_value, 
                                   2047, -2047, 0, T(1)/3000};

    for (const auto& large : larges) {
        for (const auto& small : smalls) {
            const auto a = multiply(large,small,FE_DOWNWARD);
            BOOST_CHECK_MESSAGE(a<=large*small, "product("<< large << "," << small
                                                << ", FE_DOWNWARD)=" 
                                                << a << " > " << large << "*"<< small 
                                                << "=" << large*small);
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_product_upward, T, test_types)
{
    using mantissa_type = typename IEEE754Rounding<T>::mantissa_type;

    T large_value{T(mantissa_type(1) << (IEEE754Rounding<T>::mantissa_size))};

    const std::vector<T> larges = {large_value, 1, -large_value, -1, 7, -7, 0};
    const std::vector<T> smalls = {1, T(1)/large_value, -1, -T(1)/large_value, 
                                   2047, -2047, 0, T(1)/3000};

    for (const auto& large : larges) {
        for (const auto& small : smalls) {
            const auto a = multiply(large,small,FE_UPWARD);
            BOOST_CHECK_MESSAGE(a>=large*small, "product("<< large << "," << small
                                                << ", FE_UPWARD)=" 
                                                << a << " < " << large << "*"<< small 
                                                << "=" << large*small);
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_quotient_downward, T, test_types)
{
    using mantissa_type = typename IEEE754Rounding<T>::mantissa_type;

    T large_value{T(mantissa_type(1) << (IEEE754Rounding<T>::mantissa_size))};

    const std::vector<T> larges = {large_value, 1, -large_value, -1, 7, -7, 0};
    const std::vector<T> smalls = {1, T(1)/large_value, -1, -T(1)/large_value, 
                                   2047, -2047, T(1)/3000};

    for (const auto& large : larges) {
        for (const auto& small : smalls) {
            const auto a = divide(large,small,FE_DOWNWARD);
            bool test_res;
            char wrong_relation;
            if (small>0) {
                test_res = a*small <= large;
                wrong_relation = '>'; 
            } else {
                test_res = a*small >= large;
                wrong_relation = '<'; 
            }

            BOOST_CHECK_MESSAGE(test_res, "divide("<< large << "," << small
                                           << ", FE_DOWNWARD) * " << small 
                                           << "=" << (a*small) << " "
                                           << wrong_relation << " " 
                                           << large );
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_quotient_upward, T, test_types)
{
    using mantissa_type = typename IEEE754Rounding<T>::mantissa_type;

    T large_value{T(mantissa_type(1) << (IEEE754Rounding<T>::mantissa_size))};

    const std::vector<T> larges = {large_value, 1, -large_value, -1, 7, -7, 0};
    const std::vector<T> smalls = {1, T(1)/large_value, -1, -T(1)/large_value, 
                                   2047, -2047, T(1)/3000};

    for (const auto& large : larges) {
        for (const auto& small : smalls) {
            const auto a = divide(large,small,FE_UPWARD);
            bool test_res;
            char wrong_relation;
            if (small>0) {
                test_res = a*small >= large;
                wrong_relation = '<'; 
            } else {
                test_res = a*small <= large;
                wrong_relation = '>'; 
            }

            BOOST_CHECK_MESSAGE(test_res, "divide("<< large << "," << small
                                           << ", FE_UPWARD) * " << small 
                                           << "=" << (a*small) << " "
                                           << wrong_relation << " " 
                                           << large );
        }
    }
}