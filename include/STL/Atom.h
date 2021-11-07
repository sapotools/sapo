/**
 * @file Atom.h
 * Atomic STL formula
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#ifndef ATOM_H_
#define ATOM_H_

#include <ginac/ginac.h>

#include "STL.h"

class Atom : public STL
{
private:
  static int num_of_atoms;

  GiNaC::ex predicate; // predicate
  std::vector<GiNaC::lst>
      predicateControlPts; // control points associated to this atom
  int id;                  // identifier

public:
  Atom(const GiNaC::ex &predicate);

  const GiNaC::ex &getPredicate() const
  {
    return predicate;
  };

  /**
   * Returns the control points associated with this atom
   *
   * @ returns vector of control points
   */
  const std::vector<GiNaC::lst> &getPredicateControlPts() const
  {
    return this->predicateControlPts;
  }

  /**
   * Associate a vector of control points to this atom
   *
   * @param[in] predicateControlPts vector of control points
   */
  void
  setPredicateControlPts(const std::vector<GiNaC::lst> &predicateControlPts)
  {
    this->predicateControlPts = predicateControlPts;
  }

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
