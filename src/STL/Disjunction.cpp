/**
 * @file Disjunction.cpp
 * Disjunction STL formula
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#include "Disjunction.h"

#include <iostream>

/**
 * Constructor that instantiates a Disjunction STL formula (f1 /\ f2)
 *
 * @param[in] f1 left disjunct
 * @param[in] f2 right disjunct
 */
Disjunction::Disjunction(const std::shared_ptr<STL> f1,
                         const std::shared_ptr<STL> f2):
    STL(DISJUNCTION),
    f1(f1), f2(f2)
{
}

/**
 * Print the formula
 */
std::ostream &Disjunction::print(std::ostream &os) const
{
  return os << "(" << f1 << ") || (" << f2 << ")";
}

Disjunction::~Disjunction() {}
