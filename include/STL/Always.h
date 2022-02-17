/**
 * @file Always.h
 * Always STL formula
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#ifndef ALWAYS_H_
#define ALWAYS_H_

#include <memory>

#include "STL.h"
#include "TimeInterval.h"

class Always : public STL
{

private:
  const std::shared_ptr<STL> f; // subformula
  TimeInterval t_itvl;          // temporal formula bounds

public:
  Always(const int begin, const int end, const std::shared_ptr<STL> f);

  const std::shared_ptr<STL> getSubFormula() const
  {
    return this->f;
  }

  TimeInterval time_bounds() const
  {
    return t_itvl;
  }

  void print() const;

  virtual ~Always();
};

#endif /* ALWAYS_H */
