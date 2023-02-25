#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE floating_points

#include <vector>
#include <functional>

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include "FloatingPoints.h"

typedef boost::mpl::list<float, double> test_types;

template<typename T>
inline T apply(const std::function<T (const T&, const T&, int)> f, const T& init_value, const T& value, int rounding, const size_t num_of_terms)
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

    const std::vector<T> larges = {large_value, 1, -large_value, -1, 0};
    const std::vector<T> smalls = {1, T(1)/large_value, -1, -T(1)/large_value, 0};

    constexpr T medium{1<<7};
    for (const auto& large : larges) {
        for (const auto& small : smalls) {
            T a{large};
            for (size_t i=0; i<medium; ++i) {
                a = add(a,small,FE_DOWNWARD);
            }
            BOOST_CHECK_MESSAGE(a<=large+medium*small, "sum("<<large << ","<< small<<"," << medium << ")="<< a << " > " << large << "+"<< medium<<"*"<< small);
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_sum_upward, T, test_types)
{
    using mantissa_type = typename IEEE754Rounding<T>::mantissa_type;

    T large_value{T(mantissa_type(1) << (IEEE754Rounding<T>::mantissa_size))};

    const std::vector<T> larges = {large_value, 1, -large_value, -1, 0};
    const std::vector<T> smalls = {1, T(1)/large_value, -1, -T(1)/large_value, 0};

    constexpr T medium{1<<7};
    for (const auto& large : larges) {
        for (const auto& small : smalls) {
            T a{large};
            for (size_t i=0; i<medium; ++i) {
                a = add(a,small,FE_UPWARD);
            }
            BOOST_CHECK_MESSAGE(a>=large+medium*small, "sum("<<large << ","<< small<<"," << medium << ")="<< a << " < " << large << "+"<< medium<<"*"<< small);
        }
    }
}


BOOST_AUTO_TEST_CASE_TEMPLATE(test_difference_downward, T, test_types)
{
    using mantissa_type = typename IEEE754Rounding<T>::mantissa_type;

    T large_value{T(mantissa_type(1) << (IEEE754Rounding<T>::mantissa_size))};

    const std::vector<T> larges = {large_value, 1, -large_value, -1, 0};
    const std::vector<T> smalls = {1, T(1)/large_value, -1, -T(1)/large_value, 0};

    constexpr T medium{1<<7};
    for (const auto& large : larges) {
        for (const auto& small : smalls) {
            T a{large};
            for (size_t i=0; i<medium; ++i) {
                a = subtract(a,small,FE_DOWNWARD);
            }
            BOOST_CHECK_MESSAGE(a<=large-medium*small, "difference("<<large << ","<< small<<"," << medium << ")="<< a << " > " << large << "-"<< medium<<"*"<< small << "=" << large-medium*small);
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_difference_upward, T, test_types)
{
    using mantissa_type = typename IEEE754Rounding<T>::mantissa_type;

    T large_value{T(mantissa_type(1) << (IEEE754Rounding<T>::mantissa_size))};

    const std::vector<T> larges = {large_value, 1, -large_value, -1, 0};
    const std::vector<T> smalls = {1, T(1)/large_value, -1, -T(1)/large_value, 0};

    constexpr T medium{1<<7};
    for (const auto& large : larges) {
        for (const auto& small : smalls) {
            T a{large};
            for (size_t i=0; i<medium; ++i) {
                a = subtract(a,small,FE_UPWARD);
            }
            BOOST_CHECK_MESSAGE(a>=large-medium*small, "v("<<large << ","<< small<<"," << medium << ")="<< a << " < " << large << "-"<< medium<<"*"<< small << "=" << large-medium*small);
        }
    }
}
