/**
 * @file Until.cpp
 * Until STL formula
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#include "Until.h"

/**
 * Constructor that instantiates an Until formula (f1 U_[a,b] f2)
 *
 * @param[in] f1 left subformula
 * @param[in] begin is beginning of temporal interval
 * @param[in] end is the end of temporal interval
 * @param[in] f2 right subformula
 */
Until::Until(const std::shared_ptr<STL> f1, const int begin, const int end,
             const std::shared_ptr<STL> f2):
    STL(UNTIL),
    f1(f1), f2(f2), t_itvl(begin, end)
{
}

/**
 * Print the formula
 */
std::ostream &Until::print(std::ostream &os) const
{
  return os << "(" << *f1 << ") U " << t_itvl << " (" << *f2 << ")";
}

Until::~Until() {}