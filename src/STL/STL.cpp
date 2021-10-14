/**
 * @file STL.cpp
 * STL formula
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#include "STL.h"

STL::STL(const formula_type type, const int options): type(type),  options(options) {}

const STL& STL::set_options(const int options) 
{
    this->options = options;

    return *this;
}

STL::~STL() {}