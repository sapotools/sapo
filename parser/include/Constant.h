#ifndef __CONSTANT_H__
#define __CONSTANT_H__

#include "SymbolicAlgebra.h"

namespace AbsSyn
{

class Constant
{
  friend std::ostream &operator<<(std::ostream &os, Constant &c)
  {
    return os << c.s;
  }

public:
  Constant(SymbolicAlgebra::Symbol<> n, double v): s(n), val(v) {}

  ~Constant() {}

  std::string getName() const
  {
    return s.get_symbol_name(s.get_id());
  }

  SymbolicAlgebra::Symbol<> getSymbol() const
  {
    return s;
  }

  double getValue() const
  {
    return val;
  }

protected:
  SymbolicAlgebra::Symbol<> s;
  double val;
};

}

#endif
