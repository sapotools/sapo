/**
 * @file STL.cpp
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @author Alberto Casagrande <acasagrande@units.it>
 * @brief STL formula
 * @version 0.2
 * @date 2022-05-04
 * 
 * @copyright Copyright (c) 2015-2022
 */

#include "STL.h"

namespace STL
{

STL::STL(formula_type type): _type(type) {}

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
