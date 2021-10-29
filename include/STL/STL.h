/**
 * @file STL.h
 * STL formula
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#ifndef STL_H_
#define STL_H_

#include "Common.h"
#include <string.h>

class STL
{
  const formula_type type;

protected:
  int options;

  STL(const formula_type type);

public:
  const STL &set_options(const int options);

  const formula_type &getType() const
  {
    return type;
  }

  virtual void print() const = 0;
  virtual ~STL();
};

#endif /* STL_H_ */
