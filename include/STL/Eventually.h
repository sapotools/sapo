/**
 * @file Eventually.h
 * Eventually STL formula
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#ifndef EVENTUALLY_H_
#define EVENTUALLY_H_

#include <memory>

#include "STL.h"
#include "TimeInterval.h"

class Eventually : public STL
{
private:
  std::shared_ptr<STL> f; // subformula
  TimeInterval t_itvl;          // temporal formula bounds

public:
  Eventually(const int begin, const int end, const std::shared_ptr<STL> f);

  const std::shared_ptr<STL> getSubFormula() const
  {
    return f;
  }

  TimeInterval time_bounds() const
  {
    return t_itvl;
  }

  const std::shared_ptr<STL> simplify() const
	{
		return std::make_shared<Eventually>(t_itvl.begin(), t_itvl.end(), f->simplify());
	}

  std::ostream &print(std::ostream &os) const;

  ~Eventually();
};

#endif /* Eventually_H */
