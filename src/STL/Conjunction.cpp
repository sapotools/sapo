/**
 * @file Conjunction.cpp
 * Conjunction STL formula
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#include "Conjunction.h"

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
std::ostream &Conjunction::print(std::ostream &os) const
{
  return os << "(" << f1 << ") && (" << f2 << ")";
}

Conjunction::~Conjunction() {}
