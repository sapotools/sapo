#ifndef __PARAMETER_H__
#define __PARAMETER_H__

#include "SymbolicAlgebra.h"

namespace AbsSyn
{

class Parameter
{
  friend std::ostream &operator<<(std::ostream &os, const Parameter &p)
  {
    return os << p.s;
  }

public:
  Parameter(const SymbolicAlgebra::Symbol<> &n): s(n) {}

  ~Parameter() {}

  const SymbolicAlgebra::Symbol<> &getSymbol() const
  {
    return s;
  }

  const std::string &getName() const
  {
    return s.get_symbol_name(s.get_id());
  }

  bool isCovered()
  {
    return covered;
  }

  void setCovered()
  {
    this->covered = true;
  }

protected:
  SymbolicAlgebra::Symbol<> s; // name of the parameter
  bool covered;
};

}

#endif
