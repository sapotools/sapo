/**
 * @file Disjunction.h
 * Disjunction STL formula
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#ifndef DISJUNCTION_H_
#define DISJUNCTION_H_

#include "STL.h"

class Disjunction : public STL
{

private:
  const std::shared_ptr<STL> f1, f2; // subformulas

public:
  Disjunction(const std::shared_ptr<STL> f1, const std::shared_ptr<STL> f2);

  const std::shared_ptr<STL> getLeftSubFormula() const
  {
    return f1;
  }
  const std::shared_ptr<STL> getRightSubFormula() const
  {
    return f2;
  }

  void print() const;

  ~Disjunction();
};

#endif /* CONJUNCTION_H_ */
