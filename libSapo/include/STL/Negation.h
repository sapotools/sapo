/**
 * @file Negation.h
 * Negation STL formula
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#ifndef __NEGATION_H__
#define __NEGATION_H__

#include <iostream>
#include <memory>

#include "STL.h"
#include "Atom.h"
#include "Conjunction.h"
#include "Disjunction.h"
#include "Negation.h"
#include "Always.h"
#include "Eventually.h"
#include "Until.h"
#include "TimeInterval.h"

namespace STL 
{

class Negation : public STL
{

private:
  const std::shared_ptr<STL> f; // subformula

public:
  Negation(const std::shared_ptr<STL> f);

  const std::shared_ptr<STL> getSubFormula() const
  {
    return f;
  }

  const std::shared_ptr<STL> simplify() const;

  std::ostream &print(std::ostream &os) const;

  TimeInterval time_bounds() const
  {
    return f->time_bounds();
  }

  ~Negation();
};

}
#endif /* Eventually_H */
