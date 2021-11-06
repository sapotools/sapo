/**
 * @file Conjunction.h
 * Conjunction STL formula
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#ifndef CONJUNCTION_H_
#define CONJUNCTION_H_

#include <memory>

#include "STL.h"

class Conjunction : public STL
{

private:
  const std::shared_ptr<STL> f1, f2; // subformulas

public:
  // constructor
  Conjunction(const std::shared_ptr<STL> f1, const std::shared_ptr<STL> f2);

  const std::shared_ptr<STL> getLeftSubFormula() const
  {
    return f1;
  }
  const std::shared_ptr<STL> getRightSubFormula() const
  {
    return f2;
  }

  void print() const;

  ~Conjunction();
};

#endif /* CONJUNCTION_H_ */
