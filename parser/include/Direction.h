#ifndef __DIRECTION_H__
#define __DIRECTION_H__

#include <cmath>

#include "SymbolicAlgebra.h"
#include "Expr.h"
#include "AbsSynIO.h"

namespace AbsSyn
{

class Direction
{
  friend std::ostream &operator<<(std::ostream &os, const Direction &d);

public:
  enum Type {
    LT, // <
    LE, // <=
    GT, // >
    GE, // >=
    EQ, // =
    INT // lhs in [a,b]
  };

  Direction(SymbolicAlgebra::Expression<> e1, SymbolicAlgebra::Expression<> e2,
            Type t, double lb = -std::numeric_limits<double>::infinity(),
            double ub = std::numeric_limits<double>::infinity(),
            SymbolicAlgebra::Symbol<> *sym = NULL):
      lhs(e1),
      rhs(e2), type(t), LB(lb), UB(ub), s(sym)
  {
  }

  ~Direction() {}

  std::vector<double>
  getDirectionVector(std::vector<SymbolicAlgebra::Symbol<>> symbols)
      const; // returns the vector representing the direction

  double
  getOffset() const; // return the offset of the direction, if type is not INT

  std::string getName() const
  {
    if (s != NULL) {
      return SymbolicAlgebra::Symbol<>::get_symbol_name(s->get_id());
    } else {
      throw std::runtime_error("Direction has no name");
    }
  }

  SymbolicAlgebra::Symbol<> *getSymbol() const
  {
    return s;
  }

  void setSymbol(SymbolicAlgebra::Symbol<> *sym)
  {
    s = sym;
  }

  double getLB() const;
  double getUB() const;

  bool hasLB() const
  {
    return type == Type::INT || type == Type::EQ
           || LB != -std::numeric_limits<double>::infinity();
  }
  bool hasUB() const
  {
    return type == Type::INT || type == Type::EQ
           || UB != std::numeric_limits<double>::infinity();
  }

  void setLB(double val);
  void setUB(double val);

  SymbolicAlgebra::Expression<> getLHS()
  {
    return lhs;
  }
  SymbolicAlgebra::Expression<> getRHS()
  {
    return rhs;
  }

  Type getType() const
  {
    return type;
  }

  bool contains(std::vector<SymbolicAlgebra::Symbol<>> symbols) const
  {
    return AbsSyn::contains(lhs, symbols)
           || (type != Type::INT && AbsSyn::contains(rhs, symbols));
  }

  Direction *copy() const; // deep copy of direction

  Direction *getComplementary() const; // returns the negated direction

  bool compare(Direction *d) const; // comparison between directions

  bool covers(const Symbol<> &s)
      const; // checks if the symbol named "name" is present in the direction

protected:
  SymbolicAlgebra::Expression<> lhs, rhs;
  Type type;
  double LB, UB; // used only if type is INT
  SymbolicAlgebra::Symbol<> *s;
};

}

#endif
