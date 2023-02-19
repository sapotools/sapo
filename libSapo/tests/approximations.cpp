#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE approximations

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include "Approximation.h"

#ifdef HAVE_GMP
#include <gmpxx.h>

//typedef boost::mpl::list<double, mpq_class> test_types;

typedef boost::mpl::list<double> test_types;
#else
typedef boost::mpl::list<double> test_types;
#endif

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


BOOST_AUTO_TEST_CASE_TEMPLATE(test_approximation_strictly_contains, T, test_types)
{
    Approximation<T> a(0,1), b(0,0), c(1,1), d(-1,2), e(0, 1);

    BOOST_CHECK(!a.strictly_contains(a));
    BOOST_CHECK(!a.strictly_contains(b));
    BOOST_CHECK(!a.strictly_contains(c));
    BOOST_CHECK(!a.strictly_contains(d));
    BOOST_CHECK(!a.strictly_contains(e));

    BOOST_CHECK(!b.strictly_contains(a));
    BOOST_CHECK(!b.strictly_contains(b));
    BOOST_CHECK(!b.strictly_contains(c));
    BOOST_CHECK(!b.strictly_contains(d));
    BOOST_CHECK(!b.strictly_contains(e));

    BOOST_CHECK(!c.strictly_contains(a));
    BOOST_CHECK(!c.strictly_contains(b));
    BOOST_CHECK(!c.strictly_contains(c));
    BOOST_CHECK(!c.strictly_contains(d));
    BOOST_CHECK(!c.strictly_contains(e));

    BOOST_CHECK(d.strictly_contains(a));
    BOOST_CHECK(d.strictly_contains(b));
    BOOST_CHECK(d.strictly_contains(c));
    BOOST_CHECK(!d.strictly_contains(d));
    BOOST_CHECK(d.strictly_contains(e));

    BOOST_CHECK(!e.strictly_contains(a));
    BOOST_CHECK(!e.strictly_contains(b));
    BOOST_CHECK(!e.strictly_contains(c));
    BOOST_CHECK(!e.strictly_contains(d));
    BOOST_CHECK(!e.strictly_contains(e));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_approximation_equality_inequality, T, test_types)
{
    Approximation<T> a(0,1), b(0,0), c(1,1), d(-1,2), e(0, 1);

    BOOST_CHECK(a==a);
    BOOST_CHECK(a!=b);
    BOOST_CHECK(a!=c);
    BOOST_CHECK(a!=d);
    BOOST_CHECK(a==e);

    BOOST_CHECK(b!=a);
    BOOST_CHECK(b==b);
    BOOST_CHECK(b!=c);
    BOOST_CHECK(b!=d);
    BOOST_CHECK(b!=e);

    BOOST_CHECK(c!=a);
    BOOST_CHECK(c!=b);
    BOOST_CHECK(c==c);
    BOOST_CHECK(c!=d);
    BOOST_CHECK(c!=e);

    BOOST_CHECK(d!=a);
    BOOST_CHECK(d!=b);
    BOOST_CHECK(d!=c);
    BOOST_CHECK(d==d);
    BOOST_CHECK(d!=e);

    BOOST_CHECK(e==a);
    BOOST_CHECK(e!=b);
    BOOST_CHECK(e!=c);
    BOOST_CHECK(e!=d);
    BOOST_CHECK(e==e);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_approximation_lesser, T, test_types)
{
    Approximation<T> a(0,1), b(0,0), c(1,1), d(-1,2), e(0, 1);
    T v0{0}, v1{1}, v2{2}, v3{-1};

    BOOST_CHECK(!(a<a));
    BOOST_CHECK(!(a<b));
    BOOST_CHECK(!(a<c));
    BOOST_CHECK(!(a<d));
    BOOST_CHECK(!(a<e));

    BOOST_CHECK(!(a<v0));
    BOOST_CHECK(!(a<v1));
    BOOST_CHECK(a<v2);
    BOOST_CHECK(!(a<v3));

    BOOST_CHECK(!(v0<a));
    BOOST_CHECK(!(v1<a));
    BOOST_CHECK(!(v2<a));
    BOOST_CHECK(v3<a);

    BOOST_CHECK(!(b<a));
    BOOST_CHECK(!(b<b));
    BOOST_CHECK(b<c);
    BOOST_CHECK(!(b<d));
    BOOST_CHECK(!(b<e));

    BOOST_CHECK(!(b<v0));
    BOOST_CHECK(b<v1);
    BOOST_CHECK(b<v2);
    BOOST_CHECK(!(b<v3));

    BOOST_CHECK(!(v0<b));
    BOOST_CHECK(!(v1<b));
    BOOST_CHECK(!(v2<b));
    BOOST_CHECK(v3<b);

    BOOST_CHECK(!(c<a));
    BOOST_CHECK(!(c<b));
    BOOST_CHECK(!(c<c));
    BOOST_CHECK(!(c<d));
    BOOST_CHECK(!(c<e));

    BOOST_CHECK(!(c<v0));
    BOOST_CHECK(!(c<v1));
    BOOST_CHECK(c<v2);
    BOOST_CHECK(!(c<v3));

    BOOST_CHECK(v0<c);
    BOOST_CHECK(!(v1<c));
    BOOST_CHECK(!(v2<c));
    BOOST_CHECK(v3<c);

    BOOST_CHECK(!(d<a));
    BOOST_CHECK(!(d<b));
    BOOST_CHECK(!(d<c));
    BOOST_CHECK(!(d<d));
    BOOST_CHECK(!(d<e));

    BOOST_CHECK(!(d<v0));
    BOOST_CHECK(!(d<v1));
    BOOST_CHECK(!(d<v2));
    BOOST_CHECK(!(d<v3));

    BOOST_CHECK(!(v0<d));
    BOOST_CHECK(!(v1<d));
    BOOST_CHECK(!(v2<d));
    BOOST_CHECK(!(v3<d));

    BOOST_CHECK(!(e<a));
    BOOST_CHECK(!(e<b));
    BOOST_CHECK(!(e<c));
    BOOST_CHECK(!(e<d));
    BOOST_CHECK(!(e<e));

    BOOST_CHECK(!(e<v0));
    BOOST_CHECK(!(e<v1));
    BOOST_CHECK(e<v2);
    BOOST_CHECK(!(e<v3));

    BOOST_CHECK(!(v0<e));
    BOOST_CHECK(!(v1<e));
    BOOST_CHECK(!(v2<e));
    BOOST_CHECK(v3<e);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_approximation_lesser_or_equal, T, test_types)
{
    Approximation<T> a(0,1), b(0,0), c(1,1), d(-1,2), e(0, 1);
    T v0{0}, v1{1}, v2{2}, v3{-1};

    BOOST_CHECK(!(a<=a));
    BOOST_CHECK(!(a<=b));
    BOOST_CHECK(a<=c);
    BOOST_CHECK(!(a<=d));
    BOOST_CHECK(!(a<=e));

    BOOST_CHECK(!(a<=v0));
    BOOST_CHECK(a<=v1);
    BOOST_CHECK(a<=v2);
    BOOST_CHECK(!(a<=v3));

    BOOST_CHECK(v0<=a);
    BOOST_CHECK(!(v1<=a));
    BOOST_CHECK(!(v2<=a));
    BOOST_CHECK(v3<=a);

    BOOST_CHECK(b<=a);
    BOOST_CHECK(b<=b);
    BOOST_CHECK(b<=c);
    BOOST_CHECK(!(b<d));
    BOOST_CHECK(b<=e);

    BOOST_CHECK(b<=v0);
    BOOST_CHECK(b<=v1);
    BOOST_CHECK(b<=v2);
    BOOST_CHECK(!(b<=v3));

    BOOST_CHECK(v0<=b);
    BOOST_CHECK(!(v1<=b));
    BOOST_CHECK(!(v2<=b));
    BOOST_CHECK(v3<=b);

    BOOST_CHECK(!(c<=a));
    BOOST_CHECK(!(c<=b));
    BOOST_CHECK(c<=c);
    BOOST_CHECK(!(c<=d));
    BOOST_CHECK(!(c<=e));

    BOOST_CHECK(!(c<=v0));
    BOOST_CHECK(c<=v1);
    BOOST_CHECK(c<=v2);
    BOOST_CHECK(!(c<=v3));

    BOOST_CHECK(v0<=c);
    BOOST_CHECK(v1<=c);
    BOOST_CHECK(!(v2<=c));
    BOOST_CHECK(v3<=c);

    BOOST_CHECK(!(d<=a));
    BOOST_CHECK(!(d<=b));
    BOOST_CHECK(!(d<=c));
    BOOST_CHECK(!(d<=d));
    BOOST_CHECK(!(d<=e));

    BOOST_CHECK(!(d<=v0));
    BOOST_CHECK(!(d<=v1));
    BOOST_CHECK(d<=v2);
    BOOST_CHECK(!(d<=v3));

    BOOST_CHECK(!(v0<=d));
    BOOST_CHECK(!(v1<=d));
    BOOST_CHECK(!(v2<=d));
    BOOST_CHECK(v3<=d);

    BOOST_CHECK(!(e<=a));
    BOOST_CHECK(!(e<=b));
    BOOST_CHECK(e<=c);
    BOOST_CHECK(!(e<=d));
    BOOST_CHECK(!(e<=e));

    BOOST_CHECK(!(e<=v0));
    BOOST_CHECK(e<=v1);
    BOOST_CHECK(e<=v2);
    BOOST_CHECK(!(e<=v3));

    BOOST_CHECK(v0<=e);
    BOOST_CHECK(!(v1<=e));
    BOOST_CHECK(!(v2<=e));
    BOOST_CHECK(v3<=a);
}


BOOST_AUTO_TEST_CASE_TEMPLATE(test_approximation_greater, T, test_types)
{
    Approximation<T> a(0,1), b(0,0), c(1,1), d(-1,2), e(0, 1);
    T v0{0}, v1{1}, v2{2}, v3{-1};

    BOOST_CHECK(!(a>a));
    BOOST_CHECK(!(b>a));
    BOOST_CHECK(!(c>a));
    BOOST_CHECK(!(d>a));
    BOOST_CHECK(!(e>a));

    BOOST_CHECK(!(a<v0));
    BOOST_CHECK(!(a<v1));
    BOOST_CHECK(a<v2);
    BOOST_CHECK(!(a<v3));

    BOOST_CHECK(!(a>b));
    BOOST_CHECK(!(b>b));
    BOOST_CHECK(c>b);
    BOOST_CHECK(!(d>b));
    BOOST_CHECK(!(e>b));

    BOOST_CHECK(!(b<v0));
    BOOST_CHECK(b<v1);
    BOOST_CHECK(b<v2);
    BOOST_CHECK(!(b<v3));

    BOOST_CHECK(!(c<a));
    BOOST_CHECK(!(c<b));
    BOOST_CHECK(!(c<c));
    BOOST_CHECK(!(c<d));
    BOOST_CHECK(!(c<e));

    BOOST_CHECK(!(c<v0));
    BOOST_CHECK(!(c<v1));
    BOOST_CHECK(c<v2);
    BOOST_CHECK(!(c<v3));

    BOOST_CHECK(!(d<a));
    BOOST_CHECK(!(d<b));
    BOOST_CHECK(!(d<c));
    BOOST_CHECK(!(d<d));
    BOOST_CHECK(!(d<e));

    BOOST_CHECK(!(d<v0));
    BOOST_CHECK(!(d<v1));
    BOOST_CHECK(!(d<v2));
    BOOST_CHECK(!(d<v3));

    BOOST_CHECK(!(e<a));
    BOOST_CHECK(!(e<b));
    BOOST_CHECK(!(e<c));
    BOOST_CHECK(!(e<d));
    BOOST_CHECK(!(e<e));

    BOOST_CHECK(!(e<v0));
    BOOST_CHECK(!(e<v1));
    BOOST_CHECK(e<v2);
    BOOST_CHECK(!(e<v3));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_approximation_greater_or_equal, T, test_types)
{
    Approximation<T> a(0,1), b(0,0), c(1,1), d(-1,2), e(0, 1);
    T v0{0}, v1{1}, v2{2}, v3{-1};

    BOOST_CHECK(!(a<=a));
    BOOST_CHECK(!(a<=b));
    BOOST_CHECK(a<=c);
    BOOST_CHECK(!(a<=d));
    BOOST_CHECK(!(a<=e));

    BOOST_CHECK(!(a<=v0));
    BOOST_CHECK(a<=v1);
    BOOST_CHECK(a<=v2);
    BOOST_CHECK(!(a<=v3));

    BOOST_CHECK(b<=a);
    BOOST_CHECK(b<=b);
    BOOST_CHECK(b<=c);
    BOOST_CHECK(!(b<d));
    BOOST_CHECK(b<=e);

    BOOST_CHECK(b<=v0);
    BOOST_CHECK(b<=v1);
    BOOST_CHECK(b<=v2);
    BOOST_CHECK(!(b<=v3));

    BOOST_CHECK(!(c<=a));
    BOOST_CHECK(!(c<=b));
    BOOST_CHECK(c<=c);
    BOOST_CHECK(!(c<=d));
    BOOST_CHECK(!(c<=e));

    BOOST_CHECK(!(c<=v0));
    BOOST_CHECK(c<=v1);
    BOOST_CHECK(c<=v2);
    BOOST_CHECK(!(c<=v3));

    BOOST_CHECK(!(d<=a));
    BOOST_CHECK(!(d<=b));
    BOOST_CHECK(!(d<=c));
    BOOST_CHECK(!(d<=d));
    BOOST_CHECK(!(d<=e));

    BOOST_CHECK(!(d<=v0));
    BOOST_CHECK(!(d<=v1));
    BOOST_CHECK(d<=v2);
    BOOST_CHECK(!(d<=v3));

    BOOST_CHECK(!(e<=a));
    BOOST_CHECK(!(e<=b));
    BOOST_CHECK(e<=c);
    BOOST_CHECK(!(e<=d));
    BOOST_CHECK(!(e<=e));

    BOOST_CHECK(!(e<=v0));
    BOOST_CHECK(e<=v1);
    BOOST_CHECK(e<=v2);
    BOOST_CHECK(!(e<=v3));
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

BOOST_AUTO_TEST_CASE_TEMPLATE(test_approximation_sum_subtraction, T, test_types)
{
    const std::vector<T> larges = {1e16, 1};
    const std::vector<T> smalls = {1, 1e-16};
    const std::vector<size_t> num_of_terms = {(size_t)1e5, (size_t)1e5};

    for (size_t i=0; i<larges.size(); ++i) {
        for (const T l_sign: std::vector<T>{-1,1}) {
            for (const T s_sign: std::vector<T>{-1,1}) {
                T large = l_sign*larges[i];
                T small = s_sign*smalls[i];
                auto approx = sum<T>(large, small, num_of_terms[i]);
                auto fp_sum = large+num_of_terms[i]*small;
                approx.contains(fp_sum);
                BOOST_CHECK(approx.contains(fp_sum));
                
                // verify that the standard floating point types produces the wrong answer
                if constexpr (std::is_floating_point_v<T>) {
                    Approximation<T> wrong_approx = sum<T, T>(large, small, num_of_terms[i]);

                    wrong_approx.contains(fp_sum);
                    BOOST_CHECK(!(wrong_approx.contains(fp_sum)));
                }
            }
        }
    }
}
