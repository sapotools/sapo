#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE tribool

#include <boost/test/unit_test.hpp>

#include "TriBool.h"

BOOST_AUTO_TEST_CASE(test_tribool_init)
{
    BOOST_CHECK(!(TriBool().is_true())); 
    BOOST_CHECK(TriBool().is_false());
    BOOST_CHECK(!(TriBool().is_uncertain()));
 
    BOOST_CHECK(TriBool(TriBool::TRUE).is_true());
    BOOST_CHECK(!(TriBool(TriBool::TRUE).is_false()));
    BOOST_CHECK(!(TriBool(TriBool::TRUE).is_uncertain()));

    BOOST_CHECK(!(TriBool(TriBool::FALSE).is_true()));
    BOOST_CHECK(TriBool(TriBool::FALSE).is_false());
    BOOST_CHECK(!(TriBool(TriBool::FALSE).is_uncertain()));

    BOOST_CHECK(!(TriBool(TriBool::UNCERTAIN).is_true())); 
    BOOST_CHECK(!(TriBool(TriBool::UNCERTAIN).is_false())); 
    BOOST_CHECK(TriBool(TriBool::UNCERTAIN).is_uncertain()); 

    BOOST_CHECK(TriBool(true).is_true());
    BOOST_CHECK(!(TriBool(true).is_false()));
    BOOST_CHECK(!(TriBool(true).is_uncertain()));

    BOOST_CHECK(!(TriBool(false).is_true()));
    BOOST_CHECK(TriBool(false).is_false());
    BOOST_CHECK(!(TriBool(false).is_uncertain()));

    TriBool a_true(TriBool::TRUE), a_false(TriBool::FALSE),
            a_uncertain(TriBool::UNCERTAIN);

    BOOST_CHECK(TriBool(a_true).is_true());
    BOOST_CHECK(!(TriBool(a_true).is_false()));
    BOOST_CHECK(!(TriBool(a_true).is_uncertain()));

    BOOST_CHECK(!(TriBool(a_false).is_true()));
    BOOST_CHECK(TriBool(a_false).is_false());
    BOOST_CHECK(!(TriBool(a_false).is_uncertain()));
    
    BOOST_CHECK(!(TriBool(a_uncertain).is_true())); 
    BOOST_CHECK(!(TriBool(a_uncertain).is_false())); 
    BOOST_CHECK(TriBool(a_uncertain).is_uncertain()); 
}

BOOST_AUTO_TEST_CASE(test_tribool_eq)
{
    TriBool a_true(TriBool::TRUE), a_false(TriBool::FALSE),
            a_uncertain(TriBool::UNCERTAIN);

    {
        TriBool a_true2(a_true);

        BOOST_CHECK((a_true == a_true2).is_true());
        BOOST_CHECK((a_true == a_false).is_false());
        BOOST_CHECK((a_true == a_uncertain).is_uncertain());
    }

    {
        TriBool a_false2(a_false);

        BOOST_CHECK((a_false == a_true).is_false());
        BOOST_CHECK((a_false == a_false2).is_true());
        BOOST_CHECK((a_false == a_uncertain).is_uncertain());
    }

    {
        TriBool a_uncertain2(a_uncertain);

        BOOST_CHECK((a_uncertain == a_true).is_uncertain());
        BOOST_CHECK((a_uncertain == a_false).is_uncertain());
        BOOST_CHECK((a_uncertain == a_uncertain2).is_uncertain());
    }

    BOOST_CHECK((true == a_true).is_true());
    BOOST_CHECK((true == a_false).is_false());
    BOOST_CHECK((true == a_uncertain).is_uncertain());

    BOOST_CHECK((a_true == true).is_true());
    BOOST_CHECK((a_false == true).is_false());
    BOOST_CHECK((a_uncertain == true).is_uncertain());

    BOOST_CHECK((false == a_true).is_false());
    BOOST_CHECK((false == a_false).is_true());
    BOOST_CHECK((false == a_uncertain).is_uncertain());

    BOOST_CHECK((a_true == false).is_false());
    BOOST_CHECK((a_false == false).is_true());
    BOOST_CHECK((a_uncertain == false).is_uncertain());


    BOOST_CHECK((TriBool::TRUE == a_true).is_true());
    BOOST_CHECK((TriBool::TRUE == a_false).is_false());
    BOOST_CHECK((TriBool::TRUE == a_uncertain).is_uncertain());

    BOOST_CHECK((a_true == TriBool::TRUE).is_true());
    BOOST_CHECK((a_false == TriBool::TRUE).is_false());
    BOOST_CHECK((a_uncertain == TriBool::TRUE).is_uncertain());

    BOOST_CHECK((TriBool::FALSE == a_true).is_false());
    BOOST_CHECK((TriBool::FALSE == a_false).is_true());
    BOOST_CHECK((TriBool::FALSE == a_uncertain).is_uncertain());

    BOOST_CHECK((a_true == TriBool::FALSE).is_false());
    BOOST_CHECK((a_false == TriBool::FALSE).is_true());
    BOOST_CHECK((a_uncertain == TriBool::FALSE).is_uncertain());
}

BOOST_AUTO_TEST_CASE(test_tribool_ineq)
{
    TriBool a_true(TriBool::TRUE), a_false(TriBool::FALSE),
            a_uncertain(TriBool::UNCERTAIN);

    {
        TriBool a_true2(a_true);

        BOOST_CHECK((a_true != a_true2).is_false());
        BOOST_CHECK((a_true != a_false).is_true());
        BOOST_CHECK((a_true != a_uncertain).is_uncertain());
    }

    {
        TriBool a_false2(a_false);

        BOOST_CHECK((a_false != a_true).is_true());
        BOOST_CHECK((a_false != a_false2).is_false());
        BOOST_CHECK((a_false != a_uncertain).is_uncertain());
    }

    {
        TriBool a_uncertain2(a_uncertain);

        BOOST_CHECK((a_uncertain != a_true).is_uncertain());
        BOOST_CHECK((a_uncertain != a_false).is_uncertain());
        BOOST_CHECK((a_uncertain != a_uncertain2).is_uncertain());
    }

    BOOST_CHECK((true != a_true).is_false());
    BOOST_CHECK((true != a_false).is_true());
    BOOST_CHECK((true != a_uncertain).is_uncertain());

    BOOST_CHECK((a_true != true).is_false());
    BOOST_CHECK((a_false != true).is_true());
    BOOST_CHECK((a_uncertain != true).is_uncertain());

    BOOST_CHECK((false != a_true).is_true());
    BOOST_CHECK((false != a_false).is_false());
    BOOST_CHECK((false != a_uncertain).is_uncertain());

    BOOST_CHECK((a_true != false).is_true());
    BOOST_CHECK((a_false != false).is_false());
    BOOST_CHECK((a_uncertain != false).is_uncertain());

    BOOST_CHECK((TriBool::TRUE != a_true).is_false());
    BOOST_CHECK((TriBool::TRUE != a_false).is_true());
    BOOST_CHECK((TriBool::TRUE != a_uncertain).is_uncertain());

    BOOST_CHECK((a_true != TriBool::TRUE).is_false());
    BOOST_CHECK((a_false != TriBool::TRUE).is_true());
    BOOST_CHECK((a_uncertain != TriBool::TRUE).is_uncertain());

    BOOST_CHECK((TriBool::FALSE != a_true).is_true());
    BOOST_CHECK((TriBool::FALSE != a_false).is_false());
    BOOST_CHECK((TriBool::FALSE != a_uncertain).is_uncertain());

    BOOST_CHECK((a_true != TriBool::FALSE).is_true());
    BOOST_CHECK((a_false != TriBool::FALSE).is_false());
    BOOST_CHECK((a_uncertain != TriBool::FALSE).is_uncertain());
}

BOOST_AUTO_TEST_CASE(test_tribool_neg)
{
    TriBool a_true(TriBool::TRUE), a_false(TriBool::FALSE),
            a_uncertain(TriBool::UNCERTAIN);

    
    BOOST_CHECK((!a_true).is_false());
    BOOST_CHECK((!a_false).is_true());
    BOOST_CHECK((!a_uncertain).is_uncertain());
}

BOOST_AUTO_TEST_CASE(test_tribool_and)
{
    TriBool a_true(TriBool::TRUE), a_false(TriBool::FALSE),
            a_uncertain(TriBool::UNCERTAIN);
    TriBool a_true2(TriBool::TRUE), a_false2(TriBool::FALSE),
            a_uncertain2(TriBool::UNCERTAIN);

    BOOST_CHECK((a_true&&a_true).is_true());
    BOOST_CHECK((a_true&&a_true2).is_true());
    BOOST_CHECK((a_true2&&a_true).is_true());
    BOOST_CHECK((a_true&&a_false).is_false());
    BOOST_CHECK((a_false&&a_true).is_false());
    BOOST_CHECK((a_true&&a_uncertain).is_uncertain());
    BOOST_CHECK((a_uncertain&&a_true).is_uncertain());

    BOOST_CHECK((a_false&&a_false).is_false());
    BOOST_CHECK((a_false&&a_false2).is_false());
    BOOST_CHECK((a_false2&&a_false).is_false());
    BOOST_CHECK((a_false&&a_uncertain).is_false());
    BOOST_CHECK((a_uncertain&&a_false).is_false());

    BOOST_CHECK((a_uncertain&&a_uncertain).is_uncertain());
    BOOST_CHECK((a_uncertain&&a_uncertain2).is_uncertain());
    BOOST_CHECK((a_uncertain2&&a_uncertain).is_uncertain());
}

BOOST_AUTO_TEST_CASE(test_tribool_or)
{
    TriBool a_true(TriBool::TRUE), a_false(TriBool::FALSE),
            a_uncertain(TriBool::UNCERTAIN);
    TriBool a_true2(TriBool::TRUE), a_false2(TriBool::FALSE),
            a_uncertain2(TriBool::UNCERTAIN);

    BOOST_CHECK((a_true||a_true).is_true());
    BOOST_CHECK((a_true||a_true2).is_true());
    BOOST_CHECK((a_true2||a_true).is_true());
    BOOST_CHECK((a_true||a_false).is_true());
    BOOST_CHECK((a_false||a_true).is_true());
    BOOST_CHECK((a_true||a_uncertain).is_true());
    BOOST_CHECK((a_uncertain||a_true).is_true());

    BOOST_CHECK((a_false||a_false).is_false());
    BOOST_CHECK((a_false||a_false2).is_false());
    BOOST_CHECK((a_false2||a_false).is_false());
    BOOST_CHECK((a_false||a_uncertain).is_uncertain());
    BOOST_CHECK((a_uncertain||a_false).is_uncertain());

    BOOST_CHECK((a_uncertain||a_uncertain).is_uncertain());
    BOOST_CHECK((a_uncertain||a_uncertain2).is_uncertain());
    BOOST_CHECK((a_uncertain2||a_uncertain).is_uncertain());
}
