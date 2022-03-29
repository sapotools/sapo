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
    IN  // lhs in [a,b]
  };

  Direction(SymbolicAlgebra::Expression<> e1, double lb, double ub,
            SymbolicAlgebra::Symbol<> *sym = NULL):
      lhs(e1),
      rhs(-e1.get_constant_term()), type(AbsSyn::Direction::Type::IN), LB(lb),
      UB(ub), s(sym)
  {
    lhs += rhs;

    const double const_term = rhs.evaluate<double>();
    LB += const_term;
    UB += const_term;
  }

  Direction(SymbolicAlgebra::Expression<> e1, SymbolicAlgebra::Expression<> e2,
            Type t, double lb, double ub, SymbolicAlgebra::Symbol<> *sym):
      lhs(e1),
      rhs(e2), type(t), LB(lb), UB(ub), s(sym)
  {
  }

  Direction(SymbolicAlgebra::Expression<> e1, SymbolicAlgebra::Expression<> e2,
            Type t):
      lhs(e1 - e2),
      rhs(), type(t), LB(-std::numeric_limits<double>::infinity()),
      UB(std::numeric_limits<double>::infinity())
  {
    rhs = -lhs.get_constant_term();
    lhs += rhs;

    const double rhs_value = rhs.evaluate<double>();

    switch (type) {
    case Type::LE:
    case Type::LT:
      UB = rhs_value;
      break;
    case Type::GE:
    case Type::GT:
      LB = rhs_value;
      break;
    case Type::EQ:
      UB = rhs_value;
      LB = UB;
      break;
    case Type::IN:
      throw std::domain_error("Unhandled direction type (IN)");
      break;
    default:
      throw std::domain_error("Unknown direction type");
    }
  }

  ~Direction() {}

  std::vector<double>
  getConstraintVector(std::vector<SymbolicAlgebra::Symbol<>> symbols)
      const; // returns the vector representing the constraint

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
    /* Why must be an interval or an equation? */
    return LB != -std::numeric_limits<double>::infinity();
  }

  bool hasUB() const
  {
    return UB != std::numeric_limits<double>::infinity();
  }

  void setLB(double val);
  void setUB(double val);

  const SymbolicAlgebra::Expression<> &getLHS() const
  {
    return lhs;
  }

  const SymbolicAlgebra::Expression<> &getRHS() const
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
           || (type != Type::IN && AbsSyn::contains(rhs, symbols));
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

std::ostream &operator<<(std::ostream &out,
                         const AbsSyn::Direction &direction);

#endif
