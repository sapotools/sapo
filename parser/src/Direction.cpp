#include "Direction.h"

namespace AbsSyn
{

std::ostream &operator<<(std::ostream &os, const Direction &d)
{
  if (d.s != NULL) {
    os << *d.s << ": ";
  }

  os << d.lhs;

  switch (d.type) {
  case Direction::Type::LT:
    os << " < ";
    break;
  case Direction::Type::LE:
    os << " <= ";
    break;
  case Direction::Type::GT:
    os << " > ";
    break;
  case Direction::Type::GE:
    os << " >= ";
    break;
  case Direction::Type::EQ:
    os << " = ";
    break;
  case Direction::Type::INT:
    return os << " in [" << d.LB << ", " << d.UB << "]";
  default:
    throw std::logic_error("unsupported direction type");
    break;
  }
  os << d.rhs;
  return os;
}

double Direction::getLB() const
{
  if (type == Type::INT) {
    return LB;
  } else if (type == Type::GT || type == Type::GE || type == Type::EQ) {
    return -this->getOffset();
  } else {
    return -std::numeric_limits<double>::infinity();
  }
}
double Direction::getUB() const
{
  if (type == Type::INT) {
    return UB;
  } else if (type == Type::LT || type == Type::LE || type == Type::EQ) {
    return this->getOffset();
  } else {
    return std::numeric_limits<double>::infinity();
  }
}

void Direction::setLB(double val)
{
  UB = this->getUB();

  if (this->hasUB()) {
    type = Type::INT;
  }

  LB = (val == 0 ? 0 : val);
}
void Direction::setUB(double val)
{
  LB = this->getLB();

  if (this->hasLB()) {
    type = Type::INT;
  }

  UB = (val == 0 ? 0 : val);
}

Direction *Direction::copy() const
{
  SymbolicAlgebra::Expression<> new_lhs(lhs), new_rhs(rhs);
  return new Direction(new_lhs, new_rhs, type, LB, UB, s);
}

Direction *Direction::getComplementary() const
{
  Direction::Type newType;

  switch (type) {
  case Direction::Type::LT:
    newType = Direction::Type::GT;
    break;
  case Direction::Type::LE:
    newType = Direction::Type::GE;
    break;
  case Direction::Type::GT:
    newType = Direction::Type::LT;
    break;
  case Direction::Type::GE:
    newType = Direction::Type::LE;
    break;
  case Direction::Type::EQ:
    newType = Direction::Type::EQ;
    break;
  case Direction::Type::INT:
    newType = Direction::Type::INT;
    break;
  default:
    throw std::logic_error("undefined direction type");
    break;
  }

  SymbolicAlgebra::Expression<> new_lhs(lhs), new_rhs(rhs);
  return new Direction(-new_lhs, -new_rhs, newType, -UB, -LB, s);
}

bool Direction::compare(Direction *d) const
{
  std::vector<SymbolicAlgebra::Symbol<>> symbols{};
  std::set<SymbolicAlgebra::Symbol<>> ids
      = (lhs + rhs + d->getLHS() + d->getRHS()).get_symbols();

  for (auto it = ids.begin(); it != ids.end(); it++) {
    symbols.push_back(*it);
  }

  std::vector<double> d1{}, d2{};

  for (unsigned i = 0; i < symbols.size(); i++) {
    d1.push_back(getCoefficient(lhs - rhs, symbols[i]));
    d2.push_back(getCoefficient(d->getLHS() - d->getRHS(), symbols[i]));
  }

  double tol = 1E-8;

  // compute length of vectors
  double l1 = 0, l2 = 0;
  for (unsigned i = 0; i < d1.size(); i++) {
    l1 += d1[i] * d1[i];
    l2 += d2[i] * d2[i];
  }
  l1 = std::sqrt(l1);
  l2 = std::sqrt(l2);

  // check normalized difference
  for (unsigned i = 0; i < d1.size(); i++) {
    if (abs(d1[i] / l1 - d2[i] / l2) > tol) {
      return false;
    }
  }
  return true;
}

std::vector<double> Direction::getDirectionVector(
    std::vector<SymbolicAlgebra::Symbol<>> symbols) const
{
  std::vector<double> res{};

  // if inequality is > or >=, flip sign
  int coeff;
  if (type == Type::LE || type == Type::LT) {
    coeff = 1;
  } else if (type == Type::GE || type == Type::GT) {
    coeff = -1;
  } else if (type == Type::INT) {
    coeff = 1;
  } else if (type == Type::EQ) {
    coeff = 1;
  } else {
    throw std::logic_error("Unsupported inequality type");
  }

  for (unsigned i = 0; i < symbols.size(); i++) {
    res.push_back(coeff * getCoefficient(lhs - rhs, symbols[i]));
  }

  return res;
}

double Direction::getOffset() const
{
  if (type == Type::LE || type == Type::LT) {
    SymbolicAlgebra::Expression<> e = rhs - lhs;
    return AbsSyn::getOffset(e);
  } else if (type == Type::GE || type == Type::GT) {
    SymbolicAlgebra::Expression<> e = lhs - rhs;
    return AbsSyn::getOffset(e);
  } else if (type == Type::EQ) {
    SymbolicAlgebra::Expression<> e = rhs - lhs;
    return AbsSyn::getOffset(e);
  } else {
    throw std::logic_error("unsupported inequality type");
  }
}

bool Direction::covers(const SymbolicAlgebra::Symbol<> &s) const
{
  SymbolicAlgebra::Expression<> e = lhs - rhs;
  return getCoefficient(e, s) != 0;
}

}
