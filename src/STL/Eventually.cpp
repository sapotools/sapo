/**
 * @file Eventually.cpp
 * Eventually STL formula
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#include "Eventually.h"

/**
 * Constructor that instantiates an Eventually formula (F_[a,b]f)
 *
 * @param[in] begin is beginning of temporal interval
 * @param[in] end is the end of temporal interval
 * @param[in] f subformula
 */
Eventually::Eventually(const int begin, const int end,
                       const std::shared_ptr<STL> f):
    STL(EVENTUALLY),
    f(f), t_itvl(begin, end)
{
}

/**
 * Print the formula
 */
void Eventually::print() const
{
  using namespace std;

  cout << "Eventually_" << this->t_itvl << " (";
  this->f->print();
  cout << ")";
}

Eventually::~Eventually() {}
