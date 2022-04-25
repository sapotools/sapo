#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE symbolic_algebra

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <sstream>

#include "SymbolicAlgebra.h"

typedef boost::mpl::list<double> test_types;

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
    std::vector<std::pair<std::vector<Expression<T>>, std::string>> tests{
        {{Symbol<T>("x"), 1}, "1 + x"},
        {{3, -4}, "-1"},
        {{Symbol<T>("x"), 3, -4, Symbol<T>("x"), 0, Symbol<T>("y")}, "-1 + x + x + y"},
        {{Symbol<T>("x"), 3, -4, Symbol<T>("x"), Symbol<T>("y")}, "-1 + x + x + y"},
    };

    for (auto t_it = std::begin(tests); t_it != std::end(tests); ++t_it) {
        auto e_it = std::begin(t_it->first);
        Expression<T> e = *e_it;
        for (++e_it; e_it != std::end(t_it->first); ++e_it) {
            e = e + *e_it;
        }

        std::ostringstream ss;
        ss << e;
        BOOST_REQUIRE_MESSAGE(ss.str() == t_it->second, e << " != " << t_it->second);
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_prods, T, test_types)
{
    std::vector<std::pair<std::vector<Expression<T>>, std::string>> tests{
        {{Symbol<T>("x"), 1}, "x"},
        {{1, Symbol<T>("x")}, "x"},
        {{3, -4}, "-12"},
        {{Symbol<T>("x"), 3, -4, Symbol<T>("x"), 0, Symbol<T>("y")}, "0"},
        {{Symbol<T>("x"), 3, -4, Symbol<T>("x"), Symbol<T>("y")}, "-12*x*x*y"},
    };

    for (auto t_it = std::begin(tests); t_it != std::end(tests); ++t_it) {
        auto e_it = std::begin(t_it->first);
        Expression<T> e = *e_it;
        for (++e_it; e_it != std::end(t_it->first); ++e_it) {
            e = e * *e_it;
        }

        std::ostringstream ss;
        ss << e;
        BOOST_REQUIRE_MESSAGE(ss.str() == t_it->second, e << " != " << t_it->second);
    }
}
