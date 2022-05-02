/**
 * @file STL.cpp
 * STL formula
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#include "STL.h"

namespace STL
{

STL::STL(formula_type type): type(type) {}

STL::~STL() {}

TimeInterval STL::time_bounds() const
{
   return TimeInterval(1);
}

std::ostream &operator<<(std::ostream &os, const STL &formula)
{
    return formula.print(os);
}

}
