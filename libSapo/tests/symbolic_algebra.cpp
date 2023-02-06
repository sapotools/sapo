#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE symbolic_algebra

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <gmpxx.h>

#include <sstream>

#include "SymbolicAlgebra.h"

typedef boost::mpl::list<double, mpq_class> test_types;

using namespace SymbolicAlgebra;

BOOST_AUTO_TEST_CASE_TEMPLATE(test_symbols, T, test_types)
{
    Symbol<T> x("x"), y("y");

    BOOST_REQUIRE(x.get_id()!=y.get_id());
    BOOST_REQUIRE(x.get_id()==Symbol<T>("x").get_id());
    BOOST_REQUIRE(y.get_id()==Symbol<T>("y").get_id());
    BOOST_REQUIRE(y.get_name()=="y");
    BOOST_REQUIRE(x.get_name()=="x");
    BOOST_REQUIRE(x.get_name()==Symbol<T>::get_symbol_name(x.get_id()));
    BOOST_REQUIRE_THROW(Symbol<T>(""), std::domain_error);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_constants, T, test_types)
{
    std::vector<T> tests{3, 3.4, 0, -1, static_cast<T>(7)/3};

    for (auto t_it = std::begin(tests); t_it != std::end(tests); ++t_it) {
        BOOST_REQUIRE(Expression<T>(*t_it).evaluate()==*t_it);
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_sums, T, test_types)
{
    Symbol<T> x("x"), y("y");

    std::vector<std::pair<Expression<T>, std::string>> tests{
        {x+1, "1 + x"},
        {3-4, "-1"},
        {x+3-4+x+0+y, "-1 + x + x + y"},
        {x+3-4-x+y, "-1 + x + -x + y"},
    };

    for (auto t_it = std::begin(tests); t_it != std::end(tests); ++t_it) {
        std::ostringstream ss;
        ss << t_it->first;
        BOOST_REQUIRE_MESSAGE(ss.str() == t_it->second, t_it->first << " != " << t_it->second);
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_prods, T, test_types)
{
    Symbol<T> x("x"), y("y");

    std::vector<std::pair<Expression<T>, std::string>> tests{
        {x*1, "x"},
        {1*x, "x"},
        {3*-4, "-12"},
        {3*-4*x, "-12*x"},
        {x*3*-4*x*0*y, "0"},
        {x*3*-4*-x*y, "12*x*x*y"},
        {x/y, "x/y"},
        {-x/y+2*y/x, "2*y/x + -x/y"},
    };

    for (auto t_it = std::begin(tests); t_it != std::end(tests); ++t_it) {
        std::ostringstream ss;
        ss << t_it->first;
        BOOST_REQUIRE_MESSAGE(ss.str() == t_it->second, t_it->first << " != " << t_it->second);
    }
}

template<typename T>
bool operator==(const std::set<Symbol<T>>& set_a,
                const std::set<Symbol<T>>& set_b)
{
    if (set_a.size() != set_b.size()) {
        return false;
    }

    auto a_it = std::begin(set_a);
    auto b_it = std::begin(set_b);

    for (; a_it != std::end(set_a); ++a_it, ++b_it) {
        if (a_it->get_id() != b_it->get_id()) {
            return false;
        }
    }

    return true;
}

template<typename T>
bool are_equivalent(const std::map<int, Expression<T>>& map_a,
                const std::map<int, Expression<T>>& map_b)
{
    if (map_a.size() != map_b.size()) {
        return false;
    }

    auto a_it = std::begin(map_a);
    auto b_it = std::begin(map_b);

    for (; a_it != std::end(map_a); ++a_it, ++b_it) {
        if (!((*a_it).first == (*b_it).first &&
                are_equivalent((*a_it).second,(*b_it).second))) {
            return false;
        }
    }

    return true;
}

template<typename T1, typename T2>
std::ostream& operator<<(std::ostream& out, const std::pair<T1, T2> &P)
{
    out << "(" << P.first << "," << P.second << ")";

    return out;
}

template<typename T1, typename T2>
std::ostream& operator<<(std::ostream& out, const std::map<T1, T2> &M)
{
    out << "{";

    for (auto m_it = std::begin(M); m_it != std::end(M); ++m_it) {
        if (m_it != std::begin(M)) {
             out << ",";
        }

         out << *m_it;
    }

    out << "}";
    return out;
}

template<typename T>
std::ostream& operator<<(std::ostream& out, const std::set<T> &S)
{
    out << "{";

    for (auto s_it = std::begin(S); s_it != std::end(S); ++s_it) {
        if (s_it != std::begin(S)) {
             out << ",";
        }

         out << *s_it;
    }

    out << "}";

    return out;
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_symbols_in_expr, T, test_types)
{
    Symbol<T> x("x"), y("y"), z("z");

    std::vector<std::pair<Expression<T>, std::set<Symbol<T>>>> tests{
        {3+4*y-x*y+z*(z*y)*z, {x, y, z}},
        {(3+4*y-x+z)*y, {x, y, z}},
        {y*(3+4*y), {y}},
        {3-5*x, {x}},
        {-x, {x}},
        {4+5, {}},
    };

    for (auto t_it = std::begin(tests); t_it != std::end(tests); ++t_it) {
        auto exp_symbols = t_it->first.get_symbols();
        bool b_eval = (exp_symbols == t_it->second);
        std::ostringstream ss;
        ss << t_it->first << " != " << t_it->second;
        BOOST_REQUIRE_MESSAGE(b_eval, ss.str());
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_expr_coeffs, T, test_types)
{
    Symbol<T> x("x"), y("y"), z("z");

    std::vector<std::pair<Expression<T>,
                std::vector<std::pair<Symbol<T>, std::map<int, Expression<T>>>>>> tests
    {
        {3+4*y-x*y+z*(z*y)*z, {{y, {{0,3},{1,4 - x + z*z*z}}},
                               {x, {{0,3 + 4*y + y*z*z*z},{1,-y}}},
                               {z, {{0,3 + 4*y - x*y},{3,y}}}
                               }},
        {2*(y + z*y), {{y, {{1,2*(1 + z)}}},
                       {z, {{0,2*y},{1,2*y}}},
                       {x, {{0,2*(y*z + y)}}}} },
        {x*(y + z*y), {{y, {{1,(1 + z)*x}}},
                       {z, {{0,x*y},{1,x*y}}},
                       {x, {{1,y + z*y}}} }},
        {3 + 4*y + z - x*(y + z*(z*y)*z),
                               {{y, {{0,3 + z},{1,4 - x*(1+z*z*z)}}},
                                {x, {{0,3 + 4*y + z},{1,-y-y*z*z*z}}},
                                {z, {{0,3 + 4*y-x*y},{1,1},{3,-x*y}}}
                               }},
    };

    for (auto t_it = std::begin(tests); t_it != std::end(tests); ++t_it) {
        for (auto s_it = std::begin(t_it->second); s_it != std::end(t_it->second); ++s_it) {
            auto coeffs = t_it->first.get_coeffs((*s_it).first);
            bool b_eval = are_equivalent(coeffs,(*s_it).second);
            std::ostringstream ss;
            ss << "("<< t_it->first << ").get_coeffs("<< (*s_it).first << ") = " << coeffs
               << " != " << (*s_it).second;
            BOOST_REQUIRE_MESSAGE(b_eval, ss.str());
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_equivalence, T, test_types)
{
    Symbol<T> x("x"), y("y"), z("z");

    std::vector<std::pair<std::pair<Expression<T>,Expression<T>>, bool>> tests
    {
        {{3+4*y-x*y+z*(z*y)*z, 3+4*y-x*y+z*(z*y)*z}, true},
        {{-y*(x-z*(z*z)-4)+3, -y*(x-z*(z*z)-4)+3}, true},
        {{3+4*y-x*y+z*(z*y)*z, -y*(x-z*(z*z)-4)+3}, true},
        {{2*(y + z*y), y*(2+2*z)}, true},
        {{x*(y + z*y), (z+1)*x*y}, true},
        {{3 + 4*y + z - x*(y + z*(z*y)*z), (3+z)+(4-(z*z*z+1)*x)*y}, true},
        {{3+4*y-x*y+z*(z*y)*z, (3+z)+(4-(z*z*z+1)*x)*y}, false},
        {{2*(y + z*y), (z+1)*x*y}, false},
        {{x*(y + z*y), y*(2+2*z)}, false},
        {{3 + 4*y + z - x*(y + z*(z*y)*z), -y*(x-z*(z*z)-4)+3}, false},
    };

    for (auto t_it = std::begin(tests); t_it != std::end(tests); ++t_it) {
        bool b_eval = are_equivalent(t_it->first.second,t_it->first.first)==t_it->second;
        std::ostringstream ss;
        ss << t_it->first.first << (t_it->second ? " != ": " == ") << t_it->first.second;
        BOOST_REQUIRE_MESSAGE(b_eval, ss.str());
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_rational_form, T, test_types)
{
    Symbol<T> x("x"), y("y"), z("z");

    std::vector<std::pair<Expression<T>, std::pair<Expression<T>,Expression<T>>>> tests
    {
        {3*x+4/y, {4+3*x*y, y}},
        {3*x/y+4/(x/y), {3*x*x + 4*y*y, y*x}},
        {4, {4,1}},
        {4+x/x-y, {4*x + x + -y*x, x}},
    };

    for (auto t_it = std::begin(tests); t_it != std::end(tests); ++t_it) {
        auto r_poly = t_it->first.get_rational_form();
        {
            std::ostringstream ss;
            ss << "The numerator of " << r_poly << " should be " << t_it->second.first;
            BOOST_REQUIRE_MESSAGE(are_equivalent(r_poly.get_numerator(), 
                                                 t_it->second.first), ss.str());
        }
        {
            std::ostringstream ss;
            ss << "The denominator of " << r_poly << " should be " << t_it->second.second;
            BOOST_REQUIRE_MESSAGE(are_equivalent(r_poly.get_denominator(), 
                                                 t_it->second.second), ss.str());
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_apply, T, test_types)
{
    Symbol<T> x("x"), y("y"), z("z");

    std::vector<std::pair<std::pair<Expression<T>, std::map<Symbol<T>,T>>, T>> tests
    {
        {{3+x, {{x, 1},{y, 0},{z, 2}}}, 4},
        {{3+4*x, {{x, 1},{y, 0},{z, 2}}}, 7},
        {{3+4*x-x*y, {{x, 1},{y, 0},{z, 2}}}, 7},
        {{z*z, {{x, 1},{y, 0},{z, 2}}}, 4},
        {{3+4*y-x*y+z*(z*y)*z, {{x, 1},{y, 0},{z, 2}}}, 3},
        {{3+4*y-x*y+z*(z*y)*z, {{x, 1},{y, 1},{z, 2}}}, 14},
        {{-y*(x-z*(z*z)-4)+3, {{x, 1},{y, 0},{z, 2}}}, 3},
        {{-y*(x-z*(z*z)-4)+3, {{x, 1},{y, 1},{z, 2}}}, 14},
        {{2*(y + z*y), {{y, 1},{z, 2}}}, 6},
        {{2*(y + z*y), {{y, 1},{z, 0}}}, 2},
        {{x*(y + z*y), {{x, 1},{y, 0},{z, 2}}}, 0},
        {{x*(y + z*y), {{x, 1},{y, 1},{z, 2}}}, 3},
        {{3 + 4*y + z - x*(y + z*(z*y)*z), {{x, 1},{y, 0},{z, 2}}}, 5},
        {{3 + 4*y + z - x*(y + z*(z*y)*z), {{x, 1},{y, 1},{z, 2}}}, 0},
        {{(3+z)+(4-(z*z*z+1)*x)*y, {{x, 1},{y, 0},{z, 2}}}, 5},
        {{(3+z)+(4-(z*z*z+1)*x)*y, {{x, 1},{y, 1},{z, 2}}}, 0},
    };

    for (auto t_it = std::begin(tests); t_it != std::end(tests); ++t_it) {
        auto eval = t_it->first.first.apply(t_it->first.second);
        std::ostringstream ss;
        ss << "("<< t_it->first.first << ").apply(" << t_it->first.second<< ") == "<< eval << "!=" << t_it->second;
        BOOST_REQUIRE_MESSAGE((eval == t_it->second), ss.str());
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_derivative, T, test_types)
{
    Symbol<T> x("x"), y("y"), z("z");

    std::list<std::pair<Expression<T>, std::list<std::pair<Symbol<T>, Expression<T>>>>> tests
    {
        {3+x, {{x, 1}, {y, 0}, {z, 0}}},
        {3+4*x, {{x, 4}, {y, 0}, {z, 0}}},
        {-x*y, {{x, -y}, {y, -x}, {z, 0}}},
        {3+4*x-x*y, {{x, 4-y}, {y, -x}, {z, 0}}},
        {z*z, {{x, 0}, {y, 0}, {z, 2*z}}},
        {3+4*y-x*y+z*(z*y)*z, {{x, -y}, {y, 4-x+z*z*z}, {z, 3*z*z*y}}},
        {-y*(x-z*(z*z)-4)+3, {{x, -y}, {y, 4-x+z*z*z}, {z, 3*z*z*y}}},
        {2*(y + z*y), {{x, 0}, {y, 2+2*z}, {z, 2*y}}},
        /* libSapo also supports quotient derivatives, but `are_equivalent`
           does not support rational expressions yet
        {x/y, {{x, 1/y}, {y, -x/(2*y*y)}, {z, 0}}},
        {2*(x*y + x*z*y)/(z+x), {{x, 0}, {y, 2+2*z}, {z, 2*y}}},
        */
    };

    for (auto t_it = std::begin(tests); t_it != std::end(tests); ++t_it) {
        for (auto d_it = std::begin(t_it->second); d_it != std::end(t_it->second); ++d_it) {
            auto eval = t_it->first.get_derivative_wrt(d_it->first);
            std::ostringstream ss;
            ss << "("<< t_it->first << ").get_derivative_wrt(" << d_it->first<< ") == "<< eval << "!=" << d_it->second;
            BOOST_REQUIRE_MESSAGE(are_equivalent(eval,d_it->second), ss.str());
        }
    }
}