/**
 * @file Atom.cpp
 * Atomic STL formula
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#include "Atom.h"

int Atom::num_of_atoms = 0;

/**
 * Constructor that instantiates a atomic predicate
 *
 * @param[in] predicate a symbolic expression
 * @param[in] id an identifier for the atomic formula
 */
Atom::Atom(const GiNaC::ex &predicate):
    STL(ATOM), predicate(predicate), id(this->num_of_atoms++)
{
}

Atom::~Atom()
{
  // TODO Auto-generated destructor stub
}
