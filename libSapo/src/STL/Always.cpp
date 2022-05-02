/**
 * @file Always.cpp
 * Always STL formula
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#include <iostream>

#include "Always.h"

namespace STL
{

/**
 * Constructor that instantiates an Always formula (G_[a,b]f)
 *
 * @param[in] begin is beginning of temporal interval
 * @param[in] end is the end of temporal interval
 * @param[in] f subformula
 */
Always::Always(const int begin, const int end, const std::shared_ptr<STL> f):
    STL(ALWAYS), f(f), t_itvl(begin, end)
{
}

/**
 * Print the formula
 */
std::ostream &Always::print(std::ostream &os) const
{
  os << "G " << this->t_itvl << " (" << *(this->f) << ")";

  return os;
}

Always::~Always() {}

}