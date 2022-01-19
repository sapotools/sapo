/**
 * @file Atom.h
 * Atomic STL formula
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#ifndef ATOM_H_
#define ATOM_H_

#include "SymbolicAlgebra.h"

#include "STL.h"

class Atom : public STL
{
private:
  static int num_of_atoms;

  SymbolicAlgebra::Expression<> predicate; //!< predicate
  int id;                                  //!< atom identifier

public:
  Atom(const SymbolicAlgebra::Expression<> &predicate);

  const SymbolicAlgebra::Expression<> &getPredicate() const
  {
    return predicate;
  };

  void print() const
  {
    std::cout << this->predicate << " <= 0";
  }

  int getID() const
  {
    return this->id;
  }

  ~Atom();
};

#endif /* ATOM_H_ */
