/**
 * @file Until.h
 * Until STL formula
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#ifndef UNTIL_H_
#define UNTIL_H_

#include <iostream>
#include <memory>

#include "STL.h"
#include "TimeInterval.h"

class Until : public STL
{

private:
  const std::shared_ptr<STL> f1, f2; // subformulas
  TimeInterval t_itvl;               // temporal formula bounds

public:
  Until(const std::shared_ptr<STL> f1, const int begin, const int end,
        const std::shared_ptr<STL> f2);

  const std::shared_ptr<STL> getLeftSubFormula() const
  {
    return f1;
  }
  const std::shared_ptr<STL> getRightSubFormula() const
  {
    return f2;
  }

  const TimeInterval &time_bounds() const
  {
    return t_itvl;
  }

  const std::shared_ptr<STL> simplify() const
	{
		return std::make_shared<Until>(f1->simplify(), t_itvl.begin(), t_itvl.end(), f2->simplify());
	}

  std::ostream &print(std::ostream &os) const;

  ~Until();
};

#endif /* UNTIL_H */
