#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE approximations

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include "Approximation.h"

typedef boost::mpl::list<double> test_types;

BOOST_AUTO_TEST_CASE_TEMPLATE(test_approximation_init, T, test_types)
{
    BOOST_REQUIRE_NO_THROW(Approximation<T>(0,0));
    BOOST_REQUIRE_NO_THROW(Approximation<T>(-1,0));
    BOOST_REQUIRE_NO_THROW(Approximation<T>(-2,-1));
    BOOST_REQUIRE_NO_THROW(Approximation<T>(1,2));

    BOOST_REQUIRE_THROW(Approximation<T>(0,-1), std::domain_error);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_approximation_contains, T, test_types)
{
    Approximation<T> a(0,1), b(0,0), c(1,1), d(-1,2), e(0, 1);

    BOOST_CHECK(a.contains(a));
    BOOST_CHECK(a.contains(b));
    BOOST_CHECK(a.contains(c));
    BOOST_CHECK(!a.contains(d));
    BOOST_CHECK(a.contains(e));

    BOOST_CHECK(!b.contains(a));
    BOOST_CHECK(b.contains(b));
    BOOST_CHECK(!b.contains(c));
    BOOST_CHECK(!b.contains(d));
    BOOST_CHECK(!b.contains(e));

    BOOST_CHECK(!c.contains(a));
    BOOST_CHECK(!c.contains(b));
    BOOST_CHECK(c.contains(c));
    BOOST_CHECK(!c.contains(d));
    BOOST_CHECK(!c.contains(e));

    BOOST_CHECK(d.contains(a));
    BOOST_CHECK(d.contains(b));
    BOOST_CHECK(d.contains(c));
    BOOST_CHECK(d.contains(d));
    BOOST_CHECK(d.contains(e));

    BOOST_CHECK(e.contains(a));
    BOOST_CHECK(e.contains(b));
    BOOST_CHECK(e.contains(c));
    BOOST_CHECK(!e.contains(d));
    BOOST_CHECK(e.contains(e));
}


BOOST_AUTO_TEST_CASE_TEMPLATE(test_approximation_interior_contains, T, test_types)
{
    Approximation<T> a(0,1), b(0,0), c(1,1), d(-1,2), e(0, 1);

    BOOST_CHECK(!a.interior_contains(a));
    BOOST_CHECK(!a.interior_contains(b));
    BOOST_CHECK(!a.interior_contains(c));
    BOOST_CHECK(!a.interior_contains(d));
    BOOST_CHECK(!a.interior_contains(e));

    BOOST_CHECK(!b.interior_contains(a));
    BOOST_CHECK(!b.interior_contains(b));
    BOOST_CHECK(!b.interior_contains(c));
    BOOST_CHECK(!b.interior_contains(d));
    BOOST_CHECK(!b.interior_contains(e));

    BOOST_CHECK(!c.interior_contains(a));
    BOOST_CHECK(!c.interior_contains(b));
    BOOST_CHECK(!c.interior_contains(c));
    BOOST_CHECK(!c.interior_contains(d));
    BOOST_CHECK(!c.interior_contains(e));

    BOOST_CHECK(d.interior_contains(a));
    BOOST_CHECK(d.interior_contains(b));
    BOOST_CHECK(d.interior_contains(c));
    BOOST_CHECK(!d.interior_contains(d));
    BOOST_CHECK(d.interior_contains(e));

    BOOST_CHECK(!e.interior_contains(a));
    BOOST_CHECK(!e.interior_contains(b));
    BOOST_CHECK(!e.interior_contains(c));
    BOOST_CHECK(!e.interior_contains(d));
    BOOST_CHECK(!e.interior_contains(e));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_approximation_equality_inequality, T, test_types)
{
    Approximation<T> a(0,1), b(0,0), c(1,1), d(-1,2), e(0, 1);

    BOOST_CHECK((a==a).is_true());
    BOOST_CHECK((a!=b).is_uncertain());
    BOOST_CHECK((a!=c).is_uncertain());
    BOOST_CHECK((a!=d).is_uncertain());
    BOOST_CHECK((a==e).is_uncertain());

    BOOST_CHECK((b!=a).is_uncertain());
    BOOST_CHECK((b==b).is_true());
    BOOST_CHECK((b!=c).is_true());
    BOOST_CHECK((b!=d).is_uncertain());
    BOOST_CHECK((b!=e).is_uncertain());

    BOOST_CHECK((c!=a).is_uncertain());
    BOOST_CHECK((c!=b).is_true());
    BOOST_CHECK((c==c).is_true());
    BOOST_CHECK((c!=d).is_uncertain());
    BOOST_CHECK((c!=e).is_uncertain());

    BOOST_CHECK((d!=a).is_uncertain());
    BOOST_CHECK((d!=b).is_uncertain());
    BOOST_CHECK((d!=c).is_uncertain());
    BOOST_CHECK((d==d).is_true());
    BOOST_CHECK((d!=e).is_uncertain());

    BOOST_CHECK((e==a).is_uncertain());
    BOOST_CHECK((e!=b).is_uncertain());
    BOOST_CHECK((e!=c).is_uncertain());
    BOOST_CHECK((e!=d).is_uncertain());
    BOOST_CHECK((e==e).is_true());
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_approximation_lesser, T, test_types)
{
    Approximation<T> a(0,1), b(0,0), c(1,1), d(-1,2), e(0, 1);
    T v0{0}, v1{1}, v2{2}, v3{-1};

    BOOST_CHECK((a<a).is_false());
    BOOST_CHECK((a<b).is_false());
    BOOST_CHECK((a<c).is_uncertain());
    BOOST_CHECK((a<d).is_uncertain());
    BOOST_CHECK((a<e).is_uncertain());

    BOOST_CHECK((a<v0).is_false());
    BOOST_CHECK((a<v1).is_uncertain());
    BOOST_CHECK((a<v2).is_true());
    BOOST_CHECK((a<v3).is_false());

    BOOST_CHECK((v0<a).is_uncertain());
    BOOST_CHECK((v1<a).is_false());
    BOOST_CHECK((v2<a).is_false());
    BOOST_CHECK((v3<a).is_true());

    BOOST_CHECK((b<a).is_uncertain());
    BOOST_CHECK((b<b).is_false());
    BOOST_CHECK((b<c).is_true());
    BOOST_CHECK((b<d).is_uncertain());
    BOOST_CHECK((b<e).is_uncertain());

    BOOST_CHECK((b<v0).is_false());
    BOOST_CHECK((b<v1).is_true());
    BOOST_CHECK((b<v2).is_true());
    BOOST_CHECK((b<v3).is_false());

    BOOST_CHECK((v0<b).is_false());
    BOOST_CHECK((v1<b).is_false());
    BOOST_CHECK((v2<b).is_false());
    BOOST_CHECK((v3<b).is_true());

    BOOST_CHECK((c<a).is_false());
    BOOST_CHECK((c<b).is_false());
    BOOST_CHECK((c<c).is_false());
    BOOST_CHECK((c<d).is_uncertain());
    BOOST_CHECK((c<e).is_false());

    BOOST_CHECK((c<v0).is_false());
    BOOST_CHECK((c<v1).is_false());
    BOOST_CHECK((c<v2).is_true());
    BOOST_CHECK((c<v3).is_false());

    BOOST_CHECK((v0<c).is_true());
    BOOST_CHECK((v1<c).is_false());
    BOOST_CHECK((v2<c).is_false());
    BOOST_CHECK((v3<c).is_true());

    BOOST_CHECK((d<a).is_uncertain());
    BOOST_CHECK((d<b).is_uncertain());
    BOOST_CHECK((d<c).is_uncertain());
    BOOST_CHECK((d<d).is_false());
    BOOST_CHECK((d<e).is_uncertain());

    BOOST_CHECK((d<v0).is_uncertain());
    BOOST_CHECK((d<v1).is_uncertain());
    BOOST_CHECK((d<v2).is_uncertain());
    BOOST_CHECK((d<v3).is_false());

    BOOST_CHECK((v0<d).is_uncertain());
    BOOST_CHECK((v1<d).is_uncertain());
    BOOST_CHECK((v2<d).is_false());
    BOOST_CHECK((v3<d).is_uncertain());

    BOOST_CHECK((e<a).is_uncertain());
    BOOST_CHECK((e<b).is_false());
    BOOST_CHECK((e<c).is_uncertain());
    BOOST_CHECK((e<d).is_uncertain());
    BOOST_CHECK((e<e).is_false());

    BOOST_CHECK((e<v0).is_false());
    BOOST_CHECK((e<v1).is_uncertain());
    BOOST_CHECK((e<v2).is_true());
    BOOST_CHECK((e<v3).is_false());

    BOOST_CHECK((v0<e).is_uncertain());
    BOOST_CHECK((v1<e).is_false());
    BOOST_CHECK((v2<e).is_false());
    BOOST_CHECK((v3<e).is_true());
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_approximation_lesser_or_equal, T, test_types)
{
    Approximation<T> a(0,1), b(0,0), c(1,1), d(-1,2), e(0, 1);
    T v0{0}, v1{1}, v2{2}, v3{-1};

    BOOST_CHECK((a<=a).is_true());
    BOOST_CHECK((a<=b).is_uncertain());
    BOOST_CHECK((a<=c).is_true());
    BOOST_CHECK((a<=d).is_uncertain());
    BOOST_CHECK((a<=e).is_uncertain());

    BOOST_CHECK((a<=v0).is_uncertain());
    BOOST_CHECK((a<=v1).is_true());
    BOOST_CHECK((a<=v2).is_true());
    BOOST_CHECK((a<=v3).is_false());

    BOOST_CHECK((v0<=a).is_true());
    BOOST_CHECK((v1<=a).is_uncertain());
    BOOST_CHECK((v2<=a).is_false());
    BOOST_CHECK((v3<=a).is_true());

    BOOST_CHECK((b<=a).is_true());
    BOOST_CHECK((b<=b).is_true());
    BOOST_CHECK((b<=c).is_true());
    BOOST_CHECK((b<=d).is_uncertain());
    BOOST_CHECK((b<=e).is_true());

    BOOST_CHECK((b<=v0).is_true());
    BOOST_CHECK((b<=v1).is_true());
    BOOST_CHECK((b<=v2).is_true());
    BOOST_CHECK((b<=v3).is_false());

    BOOST_CHECK((v0<=b).is_true());
    BOOST_CHECK((v1<=b).is_false());
    BOOST_CHECK((v2<=b).is_false());
    BOOST_CHECK((v3<=b).is_true());

    BOOST_CHECK((c<=a).is_uncertain());
    BOOST_CHECK((c<=b).is_false());
    BOOST_CHECK((c<=c).is_true());
    BOOST_CHECK((c<=d).is_uncertain());
    BOOST_CHECK((c<=e).is_uncertain());

    BOOST_CHECK((c<=v0).is_false());
    BOOST_CHECK((c<=v1).is_true());
    BOOST_CHECK((c<=v2).is_true());
    BOOST_CHECK((c<=v3).is_false());

    BOOST_CHECK((v0<=c).is_true());
    BOOST_CHECK((v1<=c).is_true());
    BOOST_CHECK((v2<=c).is_false());
    BOOST_CHECK((v3<=c).is_true());

    BOOST_CHECK((d<=a).is_uncertain());
    BOOST_CHECK((d<=b).is_uncertain());
    BOOST_CHECK((d<=c).is_uncertain());
    BOOST_CHECK((d<=d).is_true());
    BOOST_CHECK((d<=e).is_uncertain());

    BOOST_CHECK((d<=v0).is_uncertain());
    BOOST_CHECK((d<=v1).is_uncertain());
    BOOST_CHECK((d<=v2).is_true());
    BOOST_CHECK((d<=v3).is_uncertain());

    BOOST_CHECK((v0<=d).is_uncertain());
    BOOST_CHECK((v1<=d).is_uncertain());
    BOOST_CHECK((v2<=d).is_uncertain());
    BOOST_CHECK((v3<=d).is_true());

    BOOST_CHECK((e<=a).is_uncertain());
    BOOST_CHECK((e<=b).is_uncertain());
    BOOST_CHECK((e<=c).is_true());
    BOOST_CHECK((e<=d).is_uncertain());
    BOOST_CHECK((e<=e).is_true());

    BOOST_CHECK((e<=v0).is_uncertain());
    BOOST_CHECK((e<=v1).is_true());
    BOOST_CHECK((e<=v2).is_true());
    BOOST_CHECK((e<=v3).is_false());

    BOOST_CHECK((v0<=e).is_true());
    BOOST_CHECK((v1<=e).is_uncertain());
    BOOST_CHECK((v2<=e).is_false());
    BOOST_CHECK((v3<=e).is_true());
}


BOOST_AUTO_TEST_CASE_TEMPLATE(test_approximation_greater, T, test_types)
{
    Approximation<T> a(0,1), b(0,0), c(1,1), d(-1,2), e(0, 1);
    T v0{0}, v1{1}, v2{2}, v3{-1};

    BOOST_CHECK((a>a).is_false());
    BOOST_CHECK((b>a).is_false());
    BOOST_CHECK((c>a).is_uncertain());
    BOOST_CHECK((d>a).is_uncertain());
    BOOST_CHECK((e>a).is_uncertain());

    BOOST_CHECK((a>v0).is_uncertain());
    BOOST_CHECK((a>v1).is_false());
    BOOST_CHECK((a>v2).is_false());
    BOOST_CHECK((a>v3).is_true());

    BOOST_CHECK((v0>a).is_false());
    BOOST_CHECK((v1>a).is_uncertain());
    BOOST_CHECK((v2>a).is_true());
    BOOST_CHECK((v3>a).is_false());

    BOOST_CHECK((a>b).is_uncertain());
    BOOST_CHECK((b>b).is_false());
    BOOST_CHECK((c>b).is_true());
    BOOST_CHECK((d>b).is_uncertain());
    BOOST_CHECK((e>b).is_uncertain());

    BOOST_CHECK((b>v0).is_false());
    BOOST_CHECK((b>v1).is_false());
    BOOST_CHECK((b>v2).is_false());
    BOOST_CHECK((b>v3).is_true());

    BOOST_CHECK((v0>b).is_false());
    BOOST_CHECK((v1>b).is_true());
    BOOST_CHECK((v2>b).is_true());
    BOOST_CHECK((v3>b).is_false());

    BOOST_CHECK((c>a).is_uncertain());
    BOOST_CHECK((c>b).is_true());
    BOOST_CHECK((c>c).is_false());
    BOOST_CHECK((c>d).is_uncertain());
    BOOST_CHECK((c>e).is_uncertain());

    BOOST_CHECK((c>v0).is_true());
    BOOST_CHECK((c>v1).is_false());
    BOOST_CHECK((c>v2).is_false());
    BOOST_CHECK((c>v3).is_true());

    BOOST_CHECK((v0>c).is_false());
    BOOST_CHECK((v1>c).is_false());
    BOOST_CHECK((v2>c).is_true());
    BOOST_CHECK((v3>c).is_false());

    BOOST_CHECK((d>a).is_uncertain());
    BOOST_CHECK((d>b).is_uncertain());
    BOOST_CHECK((d>c).is_uncertain());
    BOOST_CHECK((d>d).is_false());
    BOOST_CHECK((d>e).is_uncertain());

    BOOST_CHECK((d>v0).is_uncertain());
    BOOST_CHECK((d>v1).is_uncertain());
    BOOST_CHECK((d>v2).is_false());
    BOOST_CHECK((d>v3).is_uncertain());

    BOOST_CHECK((v0>d).is_uncertain());
    BOOST_CHECK((v1>d).is_uncertain());
    BOOST_CHECK((v2>d).is_uncertain());
    BOOST_CHECK((v3>d).is_false());

    BOOST_CHECK((a>e).is_uncertain());
    BOOST_CHECK((b>e).is_false());
    BOOST_CHECK((c>e).is_uncertain());
    BOOST_CHECK((d>e).is_uncertain());
    BOOST_CHECK((e>e).is_false());

    BOOST_CHECK((e>v0).is_uncertain());
    BOOST_CHECK((e>v1).is_false());
    BOOST_CHECK((e>v2).is_false());
    BOOST_CHECK((e>v3).is_true());

    BOOST_CHECK((v0>e).is_false());
    BOOST_CHECK((v1>e).is_uncertain());
    BOOST_CHECK((v2>e).is_true());
    BOOST_CHECK((v3>e).is_false());
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_approximation_greater_or_equal, T, test_types)
{
    Approximation<T> a(0,1), b(0,0), c(1,1), d(-1,2), e(0, 1);
    T v0{0}, v1{1}, v2{2}, v3{-1};

    BOOST_CHECK((a>=a).is_true());
    BOOST_CHECK((b>=a).is_uncertain());
    BOOST_CHECK((c>=a).is_true());
    BOOST_CHECK((d>=a).is_uncertain());
    BOOST_CHECK((e>=a).is_uncertain());

    BOOST_CHECK((a>=v0).is_true());
    BOOST_CHECK((a>=v1).is_uncertain());
    BOOST_CHECK((a>=v2).is_false());
    BOOST_CHECK((a>=v3).is_true());

    BOOST_CHECK((v0>=a).is_uncertain());
    BOOST_CHECK((v1>=a).is_true());
    BOOST_CHECK((v2>=a).is_true());
    BOOST_CHECK((v3>=a).is_false());

    BOOST_CHECK((a>=b).is_true());
    BOOST_CHECK((b>=b).is_true());
    BOOST_CHECK((c>=b).is_true());
    BOOST_CHECK((d>=b).is_uncertain());
    BOOST_CHECK((e>=b).is_true());

    BOOST_CHECK((b>=v0).is_true());
    BOOST_CHECK((b>=v1).is_false());
    BOOST_CHECK((b>=v2).is_false());
    BOOST_CHECK((b>=v3).is_true());

    BOOST_CHECK((v0>=b).is_true());
    BOOST_CHECK((v1>=b).is_true());
    BOOST_CHECK((v2>=b).is_true());
    BOOST_CHECK((v3>=b).is_false());

    BOOST_CHECK((c>=a).is_true());
    BOOST_CHECK((c>=b).is_true());
    BOOST_CHECK((c>=c).is_true());
    BOOST_CHECK((c>=d).is_uncertain());
    BOOST_CHECK((c>=e).is_true());

    BOOST_CHECK((c>=v0).is_true());
    BOOST_CHECK((c>=v1).is_true());
    BOOST_CHECK((c>=v2).is_false());
    BOOST_CHECK((c>=v3).is_true());

    BOOST_CHECK((v0>=c).is_false());
    BOOST_CHECK((v1>=c).is_true());
    BOOST_CHECK((v2>=c).is_true());
    BOOST_CHECK((v3>=c).is_false());

    BOOST_CHECK((d>=a).is_uncertain());
    BOOST_CHECK((d>=b).is_uncertain());
    BOOST_CHECK((d>=c).is_uncertain());
    BOOST_CHECK((d>=d).is_true());
    BOOST_CHECK((d>=e).is_uncertain());

    BOOST_CHECK((d>=v0).is_uncertain());
    BOOST_CHECK((d>=v1).is_uncertain());
    BOOST_CHECK((d>=v2).is_uncertain());
    BOOST_CHECK((d>=v3).is_true());

    BOOST_CHECK((v0>=d).is_uncertain());
    BOOST_CHECK((v1>=d).is_uncertain());
    BOOST_CHECK((v2>=d).is_true());
    BOOST_CHECK((v3>=d).is_uncertain());

    BOOST_CHECK((a>=e).is_uncertain());
    BOOST_CHECK((b>=e).is_uncertain());
    BOOST_CHECK((c>=e).is_true());
    BOOST_CHECK((d>=e).is_uncertain());
    BOOST_CHECK((e>=e).is_true());

    BOOST_CHECK((e>=v0).is_true());
    BOOST_CHECK((e>=v1).is_uncertain());
    BOOST_CHECK((e>=v2).is_false());
    BOOST_CHECK((e>=v3).is_true());

    BOOST_CHECK((v0>=e).is_uncertain());
    BOOST_CHECK((v1>=e).is_true());
    BOOST_CHECK((v2>=e).is_true());
    BOOST_CHECK((v3>=e).is_false());
}

template<typename T, typename R=Approximation<T>>
R sum(const T& init_value, const T& value, const size_t num_of_terms)
{
    R a{init_value};
    for (size_t i=0; i<num_of_terms; ++i) {
        a += value;
    }

    return a;
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_approximation_sum, T, test_types)
{
    using mantissa_type = typename IEEE754Rounding<T>::mantissa_type;

    T large_value{T(mantissa_type(1) << (IEEE754Rounding<T>::mantissa_size+1))};

    const std::vector<T> larges = {large_value, 1, T(1)/3000, 1, large_value, 1};
    const std::vector<T> smalls = {1, T(1)/large_value, 1, T(1)/3000, T(1)/3000, 1};
    constexpr size_t num_of_terms{1<<7};

    for (size_t i=0; i<larges.size(); ++i) {
        for (const auto& l_sign: {1, -1}) {
            for (const auto& s_sign: {1, -1}) {
                const auto large = l_sign*larges[i];
                const auto small = s_sign*smalls[i];
                const auto approx = sum<T>(large, small, num_of_terms);
                auto fp_sum = large+T(num_of_terms)*small;
                approx.contains(fp_sum);
                BOOST_CHECK_MESSAGE(approx.contains(fp_sum), approx << " does not contain " << fp_sum);
                BOOST_CHECK_MESSAGE(approx.error()>0, "Null error: "<< large << "+" << num_of_terms <<"*" 
                                    << small << "= "<< large << " + Sum_{i=1}^{"<<num_of_terms<< "} " 
                                    << small );
            }
        }
    }
}

template<typename T, typename R=Approximation<T>>
R difference(const T& init_value, const T& value, const size_t num_of_terms)
{
    R a{init_value};
    for (size_t i=0; i<num_of_terms; ++i) {
        a -= value;
    }

    return a;
}


BOOST_AUTO_TEST_CASE_TEMPLATE(test_approximation_subtraction, T, test_types)
{
    using mantissa_type = typename IEEE754Rounding<T>::mantissa_type;

    T large_value{T(mantissa_type(1) << (IEEE754Rounding<T>::mantissa_size+1))};

    const std::vector<T> larges = {large_value, 1, T(1)/3000, 1, large_value, 1};
    const std::vector<T> smalls = {1, T(1)/large_value, 1, T(1)/3000, T(1)/3000, 1};
    constexpr size_t num_of_terms{1<<7};

    for (size_t i=0; i<larges.size(); ++i) {
        for (const auto& l_sign: {1, -1}) {
            for (const auto& s_sign: {1, -1}) {
                const auto large = l_sign*larges[i];
                const auto small = s_sign*smalls[i];
                const auto approx = difference<T>(large, small, num_of_terms);
                const auto fp_difference = large-T(num_of_terms)*small;
                approx.contains(fp_difference);
                BOOST_CHECK_MESSAGE(approx.contains(fp_difference), 
                                    approx << " does not contain " << fp_difference);
                
                BOOST_CHECK_MESSAGE(approx.error()>0, "Null error: "<< large << "-" << num_of_terms <<"*" 
                                    << small << "= "<< large << " - Sum_{i=1}^{"<<num_of_terms<< "} " 
                                    << small );
            }
        }
    }
}

template<typename T, typename R=Approximation<T>>
R product(const T& init_value, const T& value, const size_t num_of_terms)
{
    R a{init_value};
    for (size_t i=0; i<num_of_terms; ++i) {
        a *= value;
    }

    return a;
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_approximation_product, T, test_types)
{
    using mantissa_type = typename IEEE754Rounding<T>::mantissa_type;

    T large_value{T(mantissa_type(1) << (IEEE754Rounding<T>::mantissa_size+1))};

    const std::vector<T> larges = {large_value, 1, 5, 1023, 0, 1, large_value};
    const std::vector<T> smalls = {1, T(1)/large_value, T(1)/large_value, T(1)/3000, 
                                   T(1)/large_value, T(1)/3000, large_value};
    constexpr size_t num_of_terms{1<<7};

    for (size_t i=0; i<larges.size(); ++i) {
        for (const auto& l_sign: {1, -1}) {
            for (const auto& s_sign: {1, -1}) {
                const auto large = l_sign*larges[i];
                const auto small = s_sign*(smalls[i]);
                const auto approx = product(large, small, num_of_terms);
                const auto fp_product = large*pow(small, num_of_terms);

                BOOST_CHECK_MESSAGE(approx.contains(fp_product), 
                                    approx <<  " does not contain " << fp_product);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(test_strict_aliasing_issue)
{
    Approximation<double> a{2.51581e-08}, b{2.51581e-08};

    auto c = a-b;
    BOOST_CHECK(is_true(c==0));
}
