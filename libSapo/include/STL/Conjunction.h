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


namespace STL
{

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

  const std::shared_ptr<STL> simplify() const
  {
    return std::make_shared<Conjunction>(f1->simplify(), f2->simplify());
  }

  std::ostream &print(std::ostream &os) const;
  TimeInterval time_bounds() const;

  ~Conjunction();
};

}

#endif /* CONJUNCTION_H_ */
