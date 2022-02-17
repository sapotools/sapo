/**
 * @file Conjunction.cpp
 * Conjunction STL formula
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#include "Conjunction.h"

#include <iostream>

/**
 * Constructor that instantiates a Conjunction STL formula (f1 /\ f2)
 *
 * @param[in] f1 left conjunct
 * @param[in] f2 right conjunct
 */
Conjunction::Conjunction(const std::shared_ptr<STL> f1,
                         const std::shared_ptr<STL> f2):
    STL(CONJUNCTION),
    f1(f1), f2(f2)
{
}

/**
 * Print the formula
 */
void Conjunction::print() const
{
  using namespace std;

  cout << "(";
  this->f1->print();
  cout << ") and (";
  this->f2->print();
  cout << ")";
}

TimeInterval Conjunction::time_bounds() const
{
  const TimeInterval ti1 = f1->time_bounds();
  const TimeInterval ti2 = f2->time_bounds();

  return TimeInterval(std::min(ti1.begin(),ti2.begin()),
                      std::max(ti1.end(),ti2.end()));
}

Conjunction::~Conjunction() {}
