/**
 * @file STL.h
 * STL formula
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#ifndef STL_H_
#define STL_H_

#include <iostream>
#include <string>
#include <memory>

#include "TimeInterval.h"

namespace STL {

enum formula_type {
  ATOM,
  CONJUNCTION,
  DISJUNCTION,
  UNTIL,
  ALWAYS,
  EVENTUALLY,
  NEGATION
};

class STL
{
protected:
  formula_type type;

  STL(formula_type type);

  virtual std::ostream &print(std::ostream &os) const = 0;
public:

  inline const formula_type &get_type() const
  {
    return type;
  }

  const virtual std::shared_ptr<STL> simplify() const = 0;

  virtual TimeInterval time_bounds() const;

  virtual ~STL();

  friend std::ostream &operator<<(std::ostream &os, const STL &formula);
};

}
#endif /* STL_H_ */
