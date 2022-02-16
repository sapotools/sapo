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
	friend std::ostream &operator<<(std::ostream &os, const STL &f)
	{
		return f.print(os);
	}
	
  formula_type type;

protected:
  int options;

  STL(formula_type type);

public:
  const STL &set_options(const int options);

  const formula_type &getType() const
  {
    return type;
  }
  
  void setType(formula_type t)
	{
		type = t;
	}

  const virtual std::shared_ptr<STL> simplify() const = 0;
	
  virtual std::ostream &print(std::ostream &os) const = 0;
  virtual ~STL();
};

#endif /* STL_H_ */
