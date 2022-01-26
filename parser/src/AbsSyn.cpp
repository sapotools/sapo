#include "../include/AbsSyn.h"

#include <limits>

#include "../include/AbsSynIO.h"

#include "../../include/SymbolicAlgebra.h"

using namespace std;

namespace std
{
ostream &operator<<(ostream &os, const AbsSyn::problemType t)
{
  if (t == AbsSyn::problemType::P_UNDEF)
    return os << "UNDEF";
  else if (t == AbsSyn::problemType::REACH)
    return os << "reachability";
  else
    return os << "synthesis";
}

ostream &operator<<(ostream &os, const AbsSyn::modeType t)
{
  if (t == AbsSyn::modeType::M_UNDEF)
    return os << "UNDEF";
  else if (t == AbsSyn::modeType::BOX)
    return os << "boxes";
  else if (t == AbsSyn::modeType::PARAL)
    return os << "parallelotopes";
  else
    return os << "polytopes";
}
}
namespace AbsSyn
{

std::ostream &Expr::prettyPrint(std::ostream &os, const int level) const
{
  if (level < type)
    os << "(";

  switch (type) {
  case exprType::NUM_ATOM:
    os << val;
    break;
  case exprType::ID_ATOM:
    os << name;
    break;
  case exprType::SUM:
    left->prettyPrint(os, type);
    os << " + ";
    right->prettyPrint(os, type);
    break;
  case exprType::SUB:
    left->prettyPrint(os, type);
    os << " - ";
    right->prettyPrint(os, type);
    break;
  case exprType::MUL:
    left->prettyPrint(os, type);
    os << " * ";
    right->prettyPrint(os, type);
    break;
  case exprType::DIV:
    left->prettyPrint(os, type);
    os << " / ";
    right->prettyPrint(os, type);
    break;
  case exprType::NEG:
    os << " - ";
    left->prettyPrint(os, type);
    break;
  default:
    std::logic_error("Unsupported expression");
    break;
  }

  if (level < type)
    os << ")";

  return os;
}

Expr *Expr::mul(Expr *e)
{
  Expr *res = new Expr();
  res->left = this;
  res->right = e;
  res->type = exprType::MUL;
  return res;
}

Expr *Expr::div(Expr *e)
{
  Expr *res = new Expr();
  res->left = this;
  res->right = e;
  res->type = exprType::DIV;
  return res;
}

Expr *Expr::sum(Expr *e)
{
  Expr *res = new Expr();
  res->left = this;
  res->right = e;
  res->type = exprType::SUM;
  return res;
}

Expr *Expr::sub(Expr *e)
{
  Expr *res = new Expr();
  res->left = this;
  res->right = e;
  res->type = exprType::SUB;
  return res;
}

Expr *Expr::neg()
{
  Expr *res = new Expr();
  res->left = this;
  res->type = exprType::NEG;
  return res;
}

int Expr::getDegree(const InputData &id, bool variable) const
{
	switch (type) {
		case Expr::NUM_ATOM:
			return 0;
		case Expr::ID_ATOM:
			if ((variable && id.isVarDefined(name)) || (!variable && id.isParamDefined(name))) {
				return 1;
			} else if (id.isDefDefined(name)) {
				return id.getDef(name)->getValue()->getDegree(id, variable);
			} else {
				return 0;
			}
		case Expr::NEG:
		case Expr::DIV:			// division is allowed only if denominator is numeric
			return left->getDegree(id, variable);
		case Expr::SUM:
		case Expr::SUB:
			return std::max(left->getDegree(id, variable), right->getDegree(id, variable));
		case Expr::MUL:
			return left->getDegree(id, variable) + right->getDegree(id, variable);
		default:
			std::logic_error("Unsupported expression");
	}
	// should not get here
	return -1;
}

double Expr::getCoefficient(const InputData &id, std::string symbolName) const
{
	switch (type) {
		case Expr::NUM_ATOM:
			return 0;
		case Expr::ID_ATOM:
			if (id.isSymbolDefined(name) && name == symbolName) {
				return 1;
			} else if (id.isDefDefined(name)) {
				return id.getDef(name)->getValue()->getCoefficient(id, symbolName);
			} else {
				return 0;
			}
		case Expr::NEG:
			return -left->getCoefficient(id, symbolName);
		case Expr::DIV:			// division is allowed only if denominator is numeric
			return left->getCoefficient(id, symbolName) / right->evaluate(id);
		case Expr::SUM:
			return left->getCoefficient(id, symbolName) + right->getCoefficient(id, symbolName);
		case Expr::SUB:
			return left->getCoefficient(id, symbolName) - right->getCoefficient(id, symbolName);
		case Expr::MUL: {
			// one of the two has no variables in it
			double l = left->getCoefficient(id, symbolName);
			double r = right->getCoefficient(id, symbolName);
			
			if (l == 0 && r == 0) {
				return 0;
			} else if (l == 0) {
				return r * left->evaluate(id);
			} else {
				return l * right->evaluate(id);
			}
		}
		default:
			std::logic_error("Unsupported expression");
	}
	// should not get here
	return -1;
}

double Expr::getOffset(const InputData &id) const
{
  switch (type) {
  case Expr::NUM_ATOM:
    return val;
  case Expr::ID_ATOM:
    return 0;
  case Expr::NEG:
    return -left->getOffset(id);
  case Expr::DIV: // division is allowed only if denominator is numeric
    return left->getOffset(id) / right->evaluate(id);
  case Expr::SUM:
    return left->getOffset(id) + right->getOffset(id);
  case Expr::SUB:
    return left->getOffset(id) - right->getOffset(id);
  case Expr::MUL: {
    return left->getOffset(id) * right->getOffset(id);
  }
  default:
    std::logic_error("Unsupported expression");
  }
  // should not get here
  return -1;
}

Expr *Expr::copy() const
{
  Expr *res = NULL;

  switch (type) {
  case Expr::NUM_ATOM:
    res = new Expr(val);
    break;
  case Expr::ID_ATOM:
    res = new Expr(name);
    break;
  case Expr::NEG:
    res = left->copy()->neg();
    break;
  case Expr::MUL:
    res = left->copy()->mul(right->copy());
    break;
  case Expr::DIV:
    res = left->copy()->div(right->copy());
    break;
  case Expr::SUM:
    res = left->copy()->sum(right->copy());
    break;
  case Expr::SUB:
    res = left->copy()->sub(right->copy());
    break;
  default:
    std::logic_error("Unsupported expression");
  }

  return res;
}

bool Expr::isNumeric(const InputData &im) const
{
  if (type == exprType::ID_ATOM) {
    if (im.isDefDefined(name) && im.getDef(name)->getValue()->isNumeric(im))
      return true;
    else
      return false;
  }
  if (type == exprType::NUM_ATOM)
    return true;

  return left->isNumeric(im) && (right != NULL ? right->isNumeric(im) : true);
}

bool Expr::hasVars(const InputData &id) const
{
  if (type == exprType::ID_ATOM) {
    if (id.isVarDefined(name)) {
      return true;
		} else {
      return false;
		}
  }
  if (type == exprType::NUM_ATOM) {
    return false;
	}

  return left->hasVars(id) && (right != NULL ? right->hasVars(id) : true);
}

bool Expr::hasParams(const InputData &id) const
{
  if (type == exprType::ID_ATOM) {
    if (id.isParamDefined(name)) {
      return true;
		} else {
      return false;
		}
  }
  if (type == exprType::NUM_ATOM) {
    return false;
	}

  return left->hasParams(id) && (right != NULL ? right->hasParams(id) : true);
}

bool Expr::contains(const std::string symbolName) const
{
	if (type == exprType::ID_ATOM) {
		if (name == symbolName) {
			return true;
		} else {
			return false;
		}
	}
	if (type == exprType::NUM_ATOM) {
		return false;
	}
	
	return left->contains(symbolName) && (right != nullptr ? right->contains(symbolName) : false);
}

double Expr::evaluate(const InputData &im) const
{
  switch (type) {
  case exprType::NUM_ATOM:
    return val;
  case exprType::ID_ATOM:
    return im.getDef(name)->getValue()->evaluate(im);
  case exprType::NEG:
    return -left->evaluate(im);
  case exprType::MUL:
    return left->evaluate(im) * right->evaluate(im);
  case exprType::DIV:
    return left->evaluate(im) / right->evaluate(im);
  case exprType::SUM:
    return left->evaluate(im) + right->evaluate(im);
  case exprType::SUB:
    return left->evaluate(im) - right->evaluate(im);
  default:
    std::cerr << "Unknown expression type" << std::endl;
    exit(EXIT_FAILURE);
  }
}

SymbolicAlgebra::Expression<>
Expr::toEx(const InputData &m,
           const std::vector<SymbolicAlgebra::Symbol<>> &vars,
           const std::vector<SymbolicAlgebra::Symbol<>> &params) const
{
  using namespace SymbolicAlgebra;

  switch (type) {
  case exprType::NUM_ATOM:
    return Expression<>(val);
  case exprType::ID_ATOM:
    if (m.isVarDefined(name))
      return vars[m.getVarPos(name)];
    if (m.isParamDefined(name))
      return params[m.getParamPos(name)];
    if (m.isDefDefined(name))
      return m.getDef(name)->getValue()->toEx(m, vars, params);
    throw std::logic_error("Unknown ID");
  case exprType::SUM:
    return left->toEx(m, vars, params) + right->toEx(m, vars, params);
  case exprType::MUL:
    return left->toEx(m, vars, params) * right->toEx(m, vars, params);
  case exprType::DIV:
    return left->toEx(m, vars, params) / right->toEx(m, vars, params);
  case exprType::SUB:
    return left->toEx(m, vars, params) - right->toEx(m, vars, params);
  case exprType::NEG:
    return -left->toEx(m, vars, params);
  }

  throw std::logic_error("Unknown expression type");
}

ostream &operator<<(ostream &os, const Formula &f)
{
  switch (f.type) {
  case Formula::formulaType::ATOM:
    return os << *(f.ex) << " <= 0";
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

bool Formula::isLinear(const InputData &id) const
{
  switch (type) {
  case formulaType::ATOM:
    return ex->getDegree(id) <= 1;

  // boolean combination
  case formulaType::CONJ:
  case formulaType::DISJ:
    return f1->isLinear(id) && f2->isLinear(id);

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
      ex = f1->ex->neg();
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
Formula::toSTL(const InputData &m,
               const std::vector<SymbolicAlgebra::Symbol<>> &vars,
               const std::vector<SymbolicAlgebra::Symbol<>> &params) const
{
  using namespace SymbolicAlgebra;
  switch (type) {
  case formulaType::ATOM: {
    return std::make_shared<Atom>(ex->toEx(m, vars, params));
  }
  case formulaType::CONJ:
    return std::make_shared<Conjunction>(f1->toSTL(m, vars, params),
                                         f2->toSTL(m, vars, params));
  case formulaType::DISJ:
    return std::make_shared<Disjunction>(f1->toSTL(m, vars, params),
                                         f2->toSTL(m, vars, params));
  case formulaType::ALW:
    return std::make_shared<Always>(i.first, i.second,
                                    f1->toSTL(m, vars, params));
  case formulaType::EVENT:
    return std::make_shared<Eventually>(i.first, i.second,
                                        f1->toSTL(m, vars, params));
  case formulaType::UNTIL:
    return std::make_shared<Until>(f1->toSTL(m, vars, params), i.first,
                                   i.second, f2->toSTL(m, vars, params));
  default:
    throw std::logic_error("Unsupported formula type");
  }
}


/*
 ***************************
 *        DIRECTION        *
 ***************************
 */

std::ostream &operator<<(std::ostream &os, const Direction &d)
{
	if (d.name != "") {
		os << d.name << ": ";
	}
	
	os << *(d.lhs);
	
	switch(d.type) {
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
			throw logic_error("unsupported direction type");
			break;
	}
	os << *(d.rhs);
	return os;
}

double Direction::getLB(const InputData &id) const
{
	if (type == Type::INT) {
		return LB;
	} else if (type == Type::GT || type == Type::GE || type == Type::EQ) {
		return -this->getOffset(id);
	} else {
		return -std::numeric_limits<double>::infinity();
	}
}
double Direction::getUB(const InputData &id) const
{
	if (type == Type::INT) {
		return UB;
	} else if (type == Type::LT || type == Type::LE || type == Type::EQ) {
		return this->getOffset(id);
	} else {
		return std::numeric_limits<double>::infinity();
	}
}

void Direction::setLB(const InputData &id, double val)
{
	UB = this->getUB(id);
	
	if (this->hasUB()) {
		type = Type::INT;
	}
	
	LB = (val == 0 ? 0 : val);
}
void Direction::setUB(const InputData &id, double val)
{
	LB = this->getLB(id);
	
	if (this->hasLB()) {
		type = Type::INT;
	}
	
	UB = (val == 0 ? 0 : val);
}


Direction *Direction::copy() const
{
	return new Direction(lhs->copy(), rhs->copy(), type, LB, UB, name);
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
			throw logic_error("undefined direction type");
			break;
	}
	
	return new Direction(lhs->copy()->neg(), rhs->copy()->neg(), newType, -UB, -LB, name);
}

bool Direction::compare(Direction *d, const InputData &id, bool variable) const
{
	std::vector<SymbolicAlgebra::Symbol<>> vars{}, params{};
	for (unsigned i = 0; i < id.getVarNum(); i++) {
		SymbolicAlgebra::Symbol<> s(id.getVar(i)->getName());
		vars.push_back(s);
	}
	for (unsigned i = 0; i < id.getParamNum(); i++) {
		SymbolicAlgebra::Symbol<> s(id.getParam(i)->getName());
		params.push_back(s);
	}
	
	std::vector<double> d1{}, d2{};
	if (variable) {
		Expr *e1 = lhs->sub(rhs);
		Expr *e2 = d->getLHS()->sub(d->getRHS());
		for (unsigned i = 0; i < id.getVarNum(); i++) {
			d1.push_back(e1->getCoefficient(id, id.getVar(i)->getName()));
			d2.push_back(e2->getCoefficient(id, id.getVar(i)->getName()));
		}
	} else {
		Expr *e1 = lhs->sub(rhs);
		Expr *e2 = d->getLHS()->sub(d->getRHS());
		for (unsigned i = 0; i < id.getParamNum(); i++) {
			d1.push_back(e1->getCoefficient(id, id.getParam(i)->getName()));
			d2.push_back(e2->getCoefficient(id, id.getParam(i)->getName()));
		}
	}
	
	double tol = 1E-8;
	
	// compute length of vectors
	double l1 = 0, l2 = 0;
	for (unsigned i = 0; i < d1.size(); i++) {
		l1 += d1[i]*d1[i];
		l2 += d2[i]*d2[i];
	}
	l1 = sqrt(l1);
	l2 = sqrt(l2);

	// check normalized difference
	for (unsigned i = 0; i < d1.size(); i++) {
		if (abs(d1[i]/l1 - d2[i]/l2) > tol) {
			return false;
		}
	}
	return true;
}



std::vector<double> Direction::getDirectionVector(const InputData &id, bool variables) const
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
		throw logic_error("Unsupported inequality type");
	}
	
	if (variables) {				// Inequality has only variables
		for (unsigned i = 0; i < id.getVarNum(); i++) {
			res.push_back(coeff * lhs->sub(rhs)->getCoefficient(id, id.getVar(i)->getName()));
		}
	} else {								// Inequality has only parameters
		for (unsigned i = 0; i < id.getParamNum(); i++) {
			res.push_back(coeff * lhs->sub(rhs)->getCoefficient(id, id.getParam(i)->getName()));
		}
	}
	
	return res;
}

double Direction::getOffset(const InputData &id) const
{
	if (type == Type::LE || type == Type::LT) {
		return rhs->sub(lhs)->getOffset(id);
	} else if (type == Type::GE || type == Type::GT) {
		return lhs->sub(rhs)->getOffset(id);
	} else if (type == Type::EQ) {
		return rhs->sub(lhs)->getOffset(id);
	} else {
		throw logic_error("unsupported inequality type");
	}
}

bool Direction::covers(const InputData &id, const std::string name) const
{
	return lhs->sub(rhs)->getCoefficient(id, name) != 0;
}

/*
 ***********************
 *        MODEL        *
 ***********************
 */

InputData::InputData():
    problem(problemType::P_UNDEF), varMode(modeType::M_UNDEF),
    paramMode(modeType::M_UNDEF), iterations(0), iter_set(false),
    max_param_splits(0), presplits(false),
    max_bundle_magnitude(std::numeric_limits<double>::max()), vars(), params(),
    consts(), defs(), assumptions(), spec(NULL), directions(),
    templateMatrix(), paramDirections(),
    trans(transType::T_UNDEF), decomp(false), decomp_defined(false), alpha(-1)
{
}

InputData::~InputData()
{
  for (auto it = std::begin(vars); it != std::end(vars); ++it)
    delete *it;

  for (auto it = std::begin(params); it != std::end(params); ++it)
    delete *it;

  for (auto it = std::begin(consts); it != std::end(consts); ++it)
    delete *it;

  for (auto it = std::begin(defs); it != std::end(defs); ++it)
    delete *it;

  for (auto it = std::begin(assumptions); it != std::end(assumptions); ++it)
    delete *it;
	
  for (auto it = std::begin(directions); it != std::end(directions); ++it)
    delete *it;
	
  for (auto it = std::begin(paramDirections); it != std::end(paramDirections); ++it)
    delete *it;

  delete spec;
}
ostream &operator<<(ostream &os, const InputData &m)
{
  os << "Problem: " << m.problem << endl;
  os << "Iterations: " << m.iterations << endl;

  os << endl;
  os << "Variables: " << endl;
  for (unsigned i = 0; i < m.vars.size(); i++) {
    os << "\t" << m.vars[i]->getName() << ": " << *(m.vars[i]->getDynamic())
       << endl;
  }
  os << endl;

  os << "Parameters: " << endl;
  for (unsigned i = 0; i < m.params.size(); i++) {
    os << "\t" << m.params[i]->getName() << endl;
  }
  os << endl;

  os << "Constants: " << endl;
  for (unsigned i = 0; i < m.consts.size(); i++) {
    os << "\t" << m.consts[i]->getName() << " = " << m.consts[i]->getValue()
       << endl;
  }
  os << endl;

  os << "Defines: " << endl;
  for (unsigned i = 0; i < m.defs.size(); i++) {
    os << "\t" << m.defs[i]->getName() << " = " << *(m.defs[i]->getValue())
       << endl;
  }
  os << endl;

  os << "spec: ";
  if (m.spec == NULL)
    os << "NULL" << endl;
  else
    os << *(m.spec) << endl;
  os << endl;

  os << "assumptions: ";
  for (unsigned i = 0; i < m.assumptions.size(); i++)
    os << *(m.assumptions[i]) << endl;
  os << endl;

  os << endl;
  os << "Directions:" << endl << "{" << endl;
  for (unsigned i = 0; i < m.directions.size(); i++){
		os << "\t" << *(m.directions[i]) << endl;
	}
//    os << "\t<" << m.directions[i] << "> in [" << m.LBoffsets[i] << ", "
//       << m.UBoffsets[i] << "]" << endl;
  os << "}" << endl;

  os << endl;
  os << "Template:" << endl << "{" << endl;
  for (unsigned i = 0; i < m.templateMatrix.size(); i++)
    os << "\t{" << m.templateMatrix[i] << "}" << endl;
  os << "}" << endl;

  os << endl;
  os << "Parameter directions:" << endl << "{" << endl;
  for (unsigned i = 0; i < m.paramDirections.size(); i++) {
		os << *(m.paramDirections[i]) << (i == m.paramDirections.size() - 1 ? "" : ",") << endl;
	}
  os << "}" << endl;

  return os;
}

bool InputData::isVarDefined(const string &name) const
{
  for (unsigned i = 0; i < vars.size(); i++)
    if (vars[i]->getName() == name)
      return true;

  return false;
}

bool InputData::isParamDefined(const string &name) const
{
  for (unsigned i = 0; i < params.size(); i++)
    if (params[i]->getName() == name)
      return true;

  return false;
}

bool InputData::isConstDefined(const string &name) const
{
  for (unsigned i = 0; i < consts.size(); i++)
    if (consts[i]->getName() == name)
      return true;

  return false;
}

bool InputData::isDefDefined(const string &name) const
{
  for (unsigned i = 0; i < defs.size(); i++)
    if (defs[i]->getName() == name)
      return true;

  return false;
}

bool InputData::isDirectionDefined(const std::string &name) const
{
	for (unsigned i = 0; i < directions.size(); i++) {
		if (directions[i]->getName() == name) {
			return true;
		}
	}
	return false;
}

bool InputData::isSymbolDefined(const string &name) const
{
  return isVarDefined(name) || isParamDefined(name) || isConstDefined(name)
         || isDefDefined(name) || isDirectionDefined(name);
}

const Variable *InputData::getVar(const string &name) const
{
  for (unsigned i = 0; i < vars.size(); i++)
    if (vars[i]->getName() == name)
      return vars[i];

  return NULL;
}

Variable *InputData::getVar(const string &name)
{
  for (unsigned i = 0; i < vars.size(); i++)
    if (vars[i]->getName() == name)
      return vars[i];

  return NULL;
}

int InputData::getVarPos(const string &name) const
{
  for (unsigned i = 0; i < vars.size(); i++)
    if (vars[i]->getName() == name)
      return i;

  return -1;
}

const Parameter *InputData::getParam(const string &name) const
{
  for (unsigned i = 0; i < params.size(); i++)
    if (params[i]->getName() == name)
      return params[i];

  return NULL;
}

int InputData::getParamPos(const string &name) const
{
  for (unsigned i = 0; i < params.size(); i++)
    if (params[i]->getName() == name)
      return i;

  return -1;
}

const Constant *InputData::getConst(const string &name) const
{
  for (unsigned i = 0; i < consts.size(); i++)
    if (consts[i]->getName() == name)
      return consts[i];

  return NULL;
}

const Definition *InputData::getDef(const string &name) const
{
  // TODO: replaced this and analoguous code by using maps
  for (unsigned i = 0; i < defs.size(); i++)
    if (defs[i]->getName() == name)
      return defs[i];

  return NULL;
}

int InputData::getDefPos(const string &name) const
{
  for (unsigned i = 0; i < defs.size(); i++)
    if (defs[i]->getName() == name)
      return i;

  return -1;
}
/*
void InputData::addDirection(vector<double> d, double LB, double UB)
{
  directions.push_back(d);
  LBoffsets.push_back(LB);
  UBoffsets.push_back(UB);
}*/

int find(std::vector<Direction *> M, Direction *v, const InputData &id, bool variable = true)
{
	unsigned pos = 0;
	while (pos < M.size()) {
		if (M[pos]->compare(v, id, variable)) {
			return pos;
		} else {
			pos++;
		}
	}
	// not found
	return -1;
}

void InputData::addDirectionConstraint(Direction *d)
{
//	std::cout << "Direction " << *d << std::endl;
	
	if (d->getType() == Direction::Type::INT || d->getType() == Direction::Type::EQ) {		// constraint is an interval specification
		this->addDirectionConstraint(new Direction(d->getLHS()->copy(), new Expr(d->getUB(*this)), Direction::Type::LE, 0, 0, d->getName()));
		this->addDirectionConstraint(new Direction(d->getLHS()->copy(), new Expr(d->getLB(*this)), Direction::Type::GE, 0, 0, d->getName()));
		delete(d);
		return;
	}
	
	Direction *new_dir = d;
	Direction *negated_dir = new_dir->getComplementary();
	
//	std::cout << "direction: " << *new_dir << ", negated: " << *negated_dir << std::endl;
	
	// check if new direction is already present
	int pos = find(directions, new_dir, *this);
	int negated_pos = find(directions, negated_dir, *this);
	
	if (pos != -1) {
		
//		std::cout << "already present, pos = " << pos << std::endl;
		// direction already present
		if (!directions[pos]->hasUB() || new_dir->getUB(*this) < directions[pos]->getUB(*this)) {
			directions[pos]->setUB(*this, new_dir->getUB(*this));
			directions[pos]->setName(new_dir->getName());
		}
		if (!directions[pos]->hasLB() || new_dir->getLB(*this) > directions[pos]->getLB(*this)) {
			directions[pos]->setLB(*this, new_dir->getLB(*this));
			directions[pos]->setName(new_dir->getName());
		}
		delete(new_dir);
		delete(negated_dir);
		
	} else if (negated_pos != -1) {
		
//		std::cout << "negated present, pos = " << negated_pos << ", offset = " << -i->getOffset(*this) << std::endl;
		// negation of direction already present
		if (!directions[negated_pos]->hasUB() || negated_dir->getUB(*this) < directions[negated_pos]->getUB(*this)) {
			directions[negated_pos]->setUB(*this, negated_dir->getUB(*this));
			directions[negated_pos]->setName(negated_dir->getName());
		}
		if (!directions[negated_pos]->hasLB() || negated_dir->getLB(*this) > directions[negated_pos]->getLB(*this)) {
			directions[negated_pos]->setLB(*this, negated_dir->getLB(*this));
			directions[negated_pos]->setName(negated_dir->getName());
		}
		delete(new_dir);
		delete(negated_dir);
		
	} else {
		
//		std::cout << "not present" << std::endl;
		// new direction
		directions.push_back(new_dir);
		
		// cover variables
		for (unsigned i = 0; i < vars.size(); i++) {
			if (new_dir->covers(*this, vars[i]->getName())) {
				vars[i]->setCovered();
			}
		}
		
		delete(negated_dir);
	}
/*	std::cout << "directions: {";
	for (unsigned i = 0; i < directions.size(); i++) {
		std::cout << *directions[i] << std::endl;
	}
	std::cout << "}" << std::endl;*/
}

/*void InputData::defaultDirections()
{
  directions.resize(vars.size(), vector<double>(vars.size(), 0));
  for (unsigned i = 0; i < vars.size(); i++)
    directions[i][i] = 1;
}*/

int InputData::findDirectionPos(const std::string &name) const
{
	for (unsigned i = 0; i < directions.size(); i++) {
		if (directions[i]->getName() == name) {
			return i;
		}
	}
	return -1;
}


void InputData::defaultTemplate()
{
  templateMatrix.resize(1, vector<int>(min(vars.size(), directions.size())));
  iota(templateMatrix[0].begin(), templateMatrix[0].end(), 0);
  //	templateMatrix.resize(1, vector<int>(vars.size(), 1));
}

void InputData::addParamDirectionConstraint(Direction *d)
{
//	std::cout << "Direction " << *d << std::endl;
	
	if (d->getType() == Direction::Type::INT) {		// constraint is an interval specification
		this->addParamDirectionConstraint(new Direction(d->getLHS()->copy(), new Expr(d->getUB(*this)), Direction::Type::LE, 0, 0, d->getName()));
		this->addParamDirectionConstraint(new Direction(d->getLHS()->copy(), new Expr(d->getLB(*this)), Direction::Type::GE, 0, 0, d->getName()));
		delete(d);
		return;
	}
	
	Direction *new_dir = d;
	Direction *negated_dir = new_dir->getComplementary();
	
	// check if new direction is already present
	int pos = find(paramDirections, new_dir, *this, false);
	int negated_pos = find(paramDirections, negated_dir, *this, false);
	
	if (pos != -1) {
		// direction already present
		if (!paramDirections[pos]->hasUB() || new_dir->getUB(*this) < paramDirections[pos]->getUB(*this)) {
			paramDirections[pos]->setUB(*this, new_dir->getUB(*this));
			paramDirections[pos]->setName(new_dir->getName());
		}
		if (!paramDirections[pos]->hasLB() || new_dir->getLB(*this) > paramDirections[pos]->getLB(*this)) {
			paramDirections[pos]->setLB(*this, new_dir->getLB(*this));
			paramDirections[pos]->setName(new_dir->getName());
		}
		delete(new_dir);
		delete(negated_dir);
		
	} else if (negated_pos != -1) {
		// negation of direction already present		
		if (!paramDirections[negated_pos]->hasLB() || negated_dir->getLB(*this) > paramDirections[negated_pos]->getLB(*this)) {
			paramDirections[negated_pos]->setLB(*this, negated_dir->getLB(*this));
			paramDirections[negated_pos]->setName(negated_dir->getName());
		}
		if (!paramDirections[negated_pos]->hasUB() || negated_dir->getUB(*this) < paramDirections[pos]->getUB(*this)) {
			paramDirections[negated_pos]->setUB(*this, negated_dir->getUB(*this));
			paramDirections[negated_pos]->setName(negated_dir->getName());
		}
		delete(new_dir);
		delete(negated_dir);
		
	} else {
		// new direction
		paramDirections.push_back(new_dir);
		
		// cover parameters
		for (unsigned i = 0; i < params.size(); i++) {
			if (new_dir->covers(*this, params[i]->getName())) {
				params[i]->setCovered();
			}
		}
	}
}

/*void InputData::addParamDirection(vector<double> d, double LB, double UB)
{
  paramDirections.push_back(d);
  paramLBoffsets.push_back(LB);
  paramUBoffsets.push_back(UB);
}

void InputData::defaultParamDirections()
{
  paramDirections.resize(params.size(), vector<double>(params.size(), 0));
  for (unsigned i = 0; i < params.size(); i++)
    paramDirections[i][i] = 1;
}*/

bool InputData::check()
{
  bool res = true;

  if (!isIterationSet()) {
    cerr << "Number of iterations is mandatory" << endl;
    res = false;
  }

  // each variable must have its dynamic
  for (unsigned i = 0; i < vars.size(); i++) {
    if (!vars[i]->isDynamicDefined()) {
      cerr << "Variable " << vars[i]->getName() << " has not a dynamic"
           << endl;
      res = false;
    }
	}
	
	// each var must be covered
	for (unsigned i = 0; i < vars.size(); i++) {
    if (!vars[i]->isCovered()) {
      cerr << "Variable " << vars[i]->getName() << " is not covered by any direction"
           << endl;
      res = false;
    }
	}
	
	// for each variable, check if it appears  positively (or negatively) in any constraint.
	// If so, we conclude that it is upper (or lower) bounded
	{	// put a block to avoid namespace pollution
		std::vector<bool> hasLB(vars.size(), false), hasUB(vars.size(), false);
		for (unsigned d = 0; d < directions.size(); d++) {
			std::vector<double> dirVector = directions[d]->getDirectionVector(*this, true);
			for (unsigned v = 0; v < vars.size(); v++) {
				if ((dirVector[v] > 0 && directions[d]->hasUB()) || (dirVector[v] < 0 && directions[d]->hasLB())) {
					hasUB[v] = true;
				}
				if ((dirVector[v] < 0 && directions[d]->hasUB()) || (dirVector[v] > 0 && directions[d]->hasLB())) {
					hasLB[v] = true;
				}
			}
		}
		
		for (unsigned v = 0; v < vars.size(); v++) {
			if (!hasLB[v]) {
				cerr << "Variable " << vars[v]->getName() << " has no finite lower bound" << endl;
				res = false;
			} else if (!hasUB[v]) {
				cerr << "Variable " << vars[v]->getName() << " has no finite upper bound" << endl;
				res = false;
			}
		}
	}
	
	
	// set directions LB where needed
	if (res) {
		vector<vector<double>> A{};
		vector<double> b{};
		for (unsigned i = 0; i < directions.size(); i++) {
			// add only bounded directions
			if (directions[i]->hasUB()) {
				A.push_back(directions[i]->getDirectionVector(*this, true));
				b.push_back(directions[i]->getUB(*this));
			}
		}
		for (unsigned i = 0; i < directions.size(); i++) {
			if (directions[i]->hasLB()) {
				A.push_back(-directions[i]->getDirectionVector(*this, true));
				b.push_back(-directions[i]->getLB(*this));
			}
		}
		LinearSystem LS(A, b);
		
		std::vector<SymbolicAlgebra::Symbol<>> symbols{};
		for (unsigned i = 0; i < vars.size(); i++) {
			SymbolicAlgebra::Symbol<> s(vars[i]->getName());
			symbols.push_back(s);
		}
		
		for (unsigned i = 0; i < directions.size(); i++) {
			if (!directions[i]->hasLB()) {
				std::vector<double> dirVector = directions[i]->getDirectionVector(*this, true);
				SymbolicAlgebra::Expression<> obj_function = 0;
				for (unsigned j = 0; j < dirVector.size(); j++) {
					obj_function += dirVector[j] * symbols[j];
				}
				
				double minVal = LS.minimize(symbols, obj_function);
				directions[i]->setLB(*this, minVal);
			}
			
			if (!directions[i]->hasUB()) {
				std::vector<double> dirVector = directions[i]->getDirectionVector(*this, true);
				SymbolicAlgebra::Expression<> obj_function = 0;
				for (unsigned j = 0; j < dirVector.size(); j++) {
					obj_function += dirVector[j] * symbols[j];
				}
				
				double maxVal = LS.maximize(symbols, obj_function);
				directions[i]->setUB(*this, maxVal);
			}
		}
	}
	
	// each param must be covered
	for (unsigned i = 0; i < params.size(); i++) {
    if (!params[i]->isCovered()) {
      cerr << "Parameter " << params[i]->getName() << " is not covered by any direction"
           << endl;
      res = false;
    }
	}
	
	// for each parameter, check if it appears  positively (or negatively) in any constraint.
	// If so, we conclude that it is upper (or lower) bounded
	{	// put a block to avoid namespace pollution
//		double infinity = std::numeric_limits<double>::max();
//		double negInfinity = -std::numeric_limits<double>::infinity();
		
		std::vector<bool> hasLB(params.size(), false), hasUB(params.size(), false);
		for (unsigned d = 0; d < paramDirections.size(); d++) {
			std::vector<double> dirVector = paramDirections[d]->getDirectionVector(*this, false);
			for (unsigned p = 0; p < params.size(); p++) {
				if ((dirVector[p] < 0 && paramDirections[d]->hasLB()) || (dirVector[p] > 0 && paramDirections[d]->hasUB())) {
					hasUB[p] = true;
				}
				if ((dirVector[p] > 0 && paramDirections[d]->hasLB()) || (dirVector[p] < 0 && paramDirections[d]->hasUB())) {
					hasLB[p] = true;
				}
			}
		}
	
		for (unsigned p = 0; p < params.size(); p++) {
			if (!hasLB[p]) {
				cerr << "Parameter " << params[p]->getName() << " has no finite lower bound" << endl;
				res = false;
			} else if (!hasUB[p]) {
				cerr << "Parameter " << params[p]->getName() << " has no finite upper bound" << endl;
				res = false;
			}
		}
	}
	
	
	
	// set param directions bounds where needed
	if (res) {
		// prepare linear system
		vector<vector<double>> A{};
		vector<double> b{};
		for (unsigned i = 0; i < paramDirections.size(); i++) {
			if (paramDirections[i]->hasUB()) {
				A.push_back(paramDirections[i]->getDirectionVector(*this, false));
				b.push_back(paramDirections[i]->getUB(*this));
			}
			if (paramDirections[i]->hasLB()) {
				A.push_back(paramDirections[i]->getComplementary()->getDirectionVector(*this, false));
				b.push_back(paramDirections[i]->getLB(*this));
			}
		}
		LinearSystem LS(A, b);
		
		std::vector<SymbolicAlgebra::Symbol<>> symbols{};
		for (unsigned i = 0; i < params.size(); i++) {
			SymbolicAlgebra::Symbol<> s(params[i]->getName());
			symbols.push_back(s);
		}
		
		for (unsigned i = 0; i < paramDirections.size(); i++) {
			if (!paramDirections[i]->hasLB()) {
				SymbolicAlgebra::Expression<> obj_function = 0;
				std::vector<double> dirVector = paramDirections[i]->getDirectionVector(*this, false);
				for (unsigned j = 0; j < dirVector.size(); j++) {
					obj_function += dirVector[j] * symbols[j];
				}
				
				double minVal = LS.minimize(symbols, obj_function);
				paramDirections[i]->setLB(*this, minVal);
			}
			
			if (!paramDirections[i]->hasUB()) {
				SymbolicAlgebra::Expression<> obj_function = 0;
				std::vector<double> dirVector = paramDirections[i]->getDirectionVector(*this, false);
				for (unsigned j = 0; j < dirVector.size(); j++) {
					obj_function += dirVector[j] * symbols[j];
				}
				
				double maxVal = LS.maximize(symbols, obj_function);
				paramDirections[i]->setUB(*this, maxVal);
			}
		}
	}
	
	/* 
	 * check that the number of parameter directions equals the number
	 * of parameters. If they are less, paramter set is not bounded (TODO: check),
	 * if they are more, we would need a polytope, which is not supported
	 */
	if (paramDirections.size() < params.size()) {
		cerr << "Too few directions for parameters, set is unbounded" << endl;
		res = false;
	} else if (paramDirections.size() > params.size()) {
		cerr << "Too much directions for parameters, polytopes are not supported" << endl;
		cerr << "parameters: {" << endl;
		for (unsigned i = 0; i < params.size(); i++) {
			cerr << "\t" << *(params[i]) << endl;
		}
		cerr << "}" << endl;
		cerr << "parameter directions: {" << endl;
		for (unsigned i = 0; i < paramDirections.size(); i++) {
			cerr << "\t" << *(paramDirections[i]) << endl;
		}
		cerr << "}" << endl;
		res = false;
	}

  // directions
  /*if (varMode != modeType::BOX && directions.size() == 0) {
    if (directions.size() == 0) {
      cerr << "Directions must be provided if variable modality is " << varMode
           << endl;
      res = false;
    } else if (directions.size() < vars.size()) {
      cerr << "Number of directions must be at least equal to the number of "
              "variables"
           << endl;
      res = false;
    }
  }*/

  // template
  /*if (varMode == modeType::POLY && templateMatrix.size() == 0) {
    cerr << "Template matrix must be provided if variable modality is "
         << varMode << endl;
    res = false;
  }*/

  // param directions
  if (paramMode == modeType::PARAL && paramDirections.size() == 0) {
    if (paramDirections.size() == 0) {
      cerr << "Parameter directions must be provided if parameter modality is "
           << paramMode << endl;
      res = false;
    } else if (paramDirections.size() < params.size()) {
      cerr << "Number of parameter directions must be at least equal to the "
              "number of parameters"
           << endl;
      res = false;
    }
  }

  // specs
  if (problem == problemType::SYNTH && spec == NULL) {
    cerr << "If problem is synthesis, a formula must be provided as spec"
         << endl;
    res = false;
  }

  return res;
}
}
