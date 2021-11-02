/**
 * @file Always.cpp
 * Always STL formula
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#include "Always.h"

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
void Always::print() const
{
  using namespace std;

  cout << "always_" << this->t_itvl << " (";
  this->f->print();
  cout << ")";
}

Always::~Always() {}
