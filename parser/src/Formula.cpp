#include "Formula.h"

namespace AbsSyn
{

ostream &operator<<(std::ostream &os, const Formula &f)
{
  switch (f.type) {
  case Formula::formulaType::ATOM:
    return os << f.ex << " <= 0";
    break;
  case Formula::formulaType::CONJ:
    return os << *(f.f1) << " && " << *(f.f2);
    break;
  case Formula::formulaType::DISJ:
    return os << *(f.f1) << " || " << *(f.f2);
    break;
  case Formula::formulaType::F_NEG:
    return os << "!" << *(f.f1);
    break;
  case Formula::formulaType::ALW:
    return os << "G" << f.i << *(f.f1);
    break;
  case Formula::formulaType::EVENT:
    return os << "F" << f.i << *(f.f1);
    break;
  case Formula::formulaType::UNTIL:
    return os << *(f.f1) << "U" << f.i << *(f.f2);
    break;
  }

  return os;
}

Formula *Formula::conj(Formula *f)
{
  Formula *res = new Formula();
  res->f1 = this;
  res->f2 = f;
  res->type = formulaType::CONJ;
  return res;
}

Formula *Formula::disj(Formula *f)
{
  Formula *res = new Formula();
  res->f1 = this;
  res->f2 = f;
  res->type = formulaType::DISJ;
  return res;
}

Formula *Formula::neg()
{
  Formula *res = new Formula();
  res->f1 = this;
  res->type = formulaType::F_NEG;
  return res;
}

Formula *Formula::always(pair<int, int> in)
{
  Formula *res = new Formula();
  res->f1 = this;
  res->i = in;
  res->type = formulaType::ALW;
  return res;
}

Formula *Formula::eventually(pair<int, int> in)
{
  Formula *res = new Formula();
  res->f1 = this;
  res->i = in;
  res->type = formulaType::EVENT;
  return res;
}

Formula *Formula::until(pair<int, int> in, Formula *f)
{
  Formula *res = new Formula();
  res->f1 = this;
  res->f2 = f;
  res->i = in;
  res->type = formulaType::UNTIL;
  return res;
}

bool Formula::isLinear(const Context &ctx) const
{
  switch (type) {
  case formulaType::ATOM:
    return getVarDegree(ex, ctx) <= 1;

  // boolean combination
  case formulaType::CONJ:
  case formulaType::DISJ:
    return f1->isLinear(ctx) && f2->isLinear(ctx);

  // temporal operators, not boolean combinations
  case formulaType::ALW:
  case formulaType::EVENT:
  case formulaType::UNTIL:
    return false;

  default:
    std::logic_error("Unsupported formula");
  }
  // should not get here
  return false;
}

bool Formula::simplify()
{
  int v = this->simplifyRec();
  while (v != 0) {
    if (v == 2)
      return false;

    v = this->simplifyRec();
  }

  return true;
}

int Formula::simplifyRec()
{
  if (type == formulaType::F_NEG) {
    Formula *child = f1;
    switch (f1->type) {
    case formulaType::ATOM:
      type = formulaType::ATOM;
      ex = -f1->ex;
      free(f1);
      f1 = NULL;
      break;
    case formulaType::CONJ:
      type = DISJ;
      f1 = f1->neg();
      f2 = f2->neg();
      break;
    case formulaType::DISJ:
      type = CONJ;
      f1 = f1->neg();
      f2 = f2->neg();
      break;
    case formulaType::F_NEG:
      type = child->f1->type;
      f1 = child->f1->f1;
      f2 = child->f1->f1;
      ex = child->f1->ex;
      i = child->f1->i;
      free(child->f1);
      free(child);
      break;
    case formulaType::ALW:
      type = formulaType::EVENT;
      i = f1->i;
      f1->type = formulaType::F_NEG;
      break;
    case formulaType::EVENT:
      type = formulaType::ALW;
      i = f1->i;
      f1->type = formulaType::F_NEG;
      break;
    case formulaType::UNTIL:
      return 2;
    }
    int v1 = f1 == NULL ? 0 : f1->simplifyRec();
    int v2 = f2 == NULL ? 0 : f2->simplifyRec();
    int max = 1;
    if (v1 > max)
      max = v1;
    if (v2 > max)
      max = v2;
    return max;
  } else if (type == formulaType::ATOM)
    return 0;
  else {
    int v1 = f1 == NULL ? 0 : f1->simplifyRec();
    int v2 = f2 == NULL ? 0 : f2->simplifyRec();
    return v1 > v2 ? v1 : v2;
  }
}

// no negations, were eliminated in simplify
std::shared_ptr<STL>
Formula::toSTL(const Context &ctx,
               const std::vector<SymbolicAlgebra::Symbol<>> &vars,
               const std::vector<SymbolicAlgebra::Symbol<>> &params) const
{
  using namespace SymbolicAlgebra;
  switch (type) {
  case formulaType::ATOM: {
    return std::make_shared<Atom>(ex);
  }
  case formulaType::CONJ:
    return std::make_shared<Conjunction>(f1->toSTL(ctx, vars, params),
                                         f2->toSTL(ctx, vars, params));
  case formulaType::DISJ:
    return std::make_shared<Disjunction>(f1->toSTL(ctx, vars, params),
                                         f2->toSTL(ctx, vars, params));
  case formulaType::ALW:
    return std::make_shared<Always>(i.first, i.second,
                                    f1->toSTL(ctx, vars, params));
  case formulaType::EVENT:
    return std::make_shared<Eventually>(i.first, i.second,
                                        f1->toSTL(ctx, vars, params));
  case formulaType::UNTIL:
    return std::make_shared<Until>(f1->toSTL(ctx, vars, params), i.first,
                                   i.second, f2->toSTL(ctx, vars, params));
  default:
    throw std::logic_error("Unsupported formula type");
  }
}

}
