#include "../include/AbsSyn.h"

#include "../include/AbsSynIO.h"

using namespace std;
using namespace GiNaC;


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

ostream &operator<<(ostream &os, const Expr &e)
{
  switch (e.type) {
  case Expr::exprType::NUM_ATOM:
    return os << e.val;
    break;
  case Expr::exprType::ID_ATOM:
    return os << e.name;
    break;
  case Expr::exprType::SUM:
    return os << *(e.left) << " + " << *(e.right);
    break;
  case Expr::exprType::SUB:
    return os << *(e.left) << " - " << *(e.right);
    break;
  case Expr::exprType::MUL:
    return os << *(e.left) << " * " << *(e.right);
    break;
  case Expr::exprType::DIV:
    return os << *(e.left) << " / " << *(e.right);
    break;
  case Expr::exprType::NEG:
    return os << "-(" << *(e.left) << ")";
    break;
  }
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

ex Expr::toEx(const InputData &m, const lst &vars, const lst &params) const
{
  switch (type) {
  case exprType::NUM_ATOM:
    return numeric(std::to_string(val).c_str());
  case exprType::ID_ATOM:
    if (m.isVarDefined(name))
      return vars[m.getVarPos(name)];
    else if (m.isParamDefined(name))
      return params[m.getParamPos(name)];
    else if (m.isDefDefined(name))
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
std::shared_ptr<STL> Formula::toSTL(const InputData &m, const lst &vars,
                                    const lst &params) const
{
  switch (type) {
  case formulaType::ATOM:
    return std::make_shared<Atom>(ex->toEx(m, vars, params));
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
 ***********************
 *        MODEL        *
 ***********************
 */

InputData::InputData()
{
  problem = problemType::P_UNDEF;
  varMode = modeType::M_UNDEF;
  paramMode = modeType::M_UNDEF;

  iterations = 0;
  iter_set = false;

  max_param_splits = 0;

  spec = NULL;

  vars.resize(0);
  params.resize(0);
  consts.resize(0);
  defs.resize(0);

  directions.resize(0);
  LBoffsets.resize(0);
  UBoffsets.resize(0);

  templateMatrix.resize(0);

  paramDirections.resize(0);
  paramLBoffsets.resize(0);
  paramUBoffsets.resize(0);

  trans = transType::T_UNDEF;
  decomp = false;
  decomp_defined = false;
  alpha = -1;
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

  delete spec;
}
ostream &operator<<(ostream &os, const InputData &m)
{
  // TODO: implement
  os << "Problem: " << m.problem << endl;
  os << "VarMode: " << m.varMode << endl;
  os << "ParamMode: " << m.paramMode << endl;
  os << "Iterations: " << m.iterations << endl;

  os << endl;
  os << "Variables: " << endl;
  for (unsigned i = 0; i < m.vars.size(); i++)
    os << "\t" << m.vars[i]->getName() << ": " << *(m.vars[i]->getDynamic())
       << endl;
  os << endl;

  os << "Parameters: " << endl;
  for (unsigned i = 0; i < m.params.size(); i++)
    os << "\t" << m.params[i]->getName() << endl;
  os << endl;

  os << "Constants: " << endl;
  for (unsigned i = 0; i < m.consts.size(); i++)
    os << "\t" << m.consts[i]->getName() << " = " << m.consts[i]->getValue()
       << endl;
  os << endl;

  os << endl;
  os << "spec: ";
  if (m.spec == NULL)
    os << "NULL" << endl;
  else
    os << *(m.spec) << endl;

  os << endl;
  os << "Directions:" << endl << "{" << endl;
  for (unsigned i = 0; i < m.directions.size(); i++)
    os << "\t<" << m.directions[i] << "> in [" << m.LBoffsets[i] << ", "
       << m.UBoffsets << "]" << endl;
  os << "}" << endl;

  os << endl;
  os << "Template:" << endl << "{" << endl;
  for (unsigned i = 0; i < m.templateMatrix.size(); i++)
    os << "\t{" << m.templateMatrix[i] << "}" << endl;
  os << "}" << endl;

  os << endl;
  os << "Parameter directions:" << endl << "{" << endl;
  for (unsigned i = 0; i < m.paramDirections.size(); i++)
    os << "\t<" << m.paramDirections[i] << "> in [" << m.paramLBoffsets[i]
       << ", " << m.paramUBoffsets << "]" << endl;
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

bool InputData::isSymbolDefined(const string &name) const
{
  return isVarDefined(name) || isParamDefined(name) || isConstDefined(name)
         || isDefDefined(name);
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

void InputData::addDirection(vector<double> d, double LB, double UB)
{
  directions.push_back(d);
  LBoffsets.push_back(LB);
  UBoffsets.push_back(UB);
}

void InputData::defaultDirections()
{
  directions.resize(vars.size(), vector<double>(vars.size(), 0));
  for (unsigned i = 0; i < vars.size(); i++)
    directions[i][i] = 1;
}

void InputData::defaultTemplate()
{
  templateMatrix.resize(1, vector<int>(vars.size()));
  iota(templateMatrix[0].begin(), templateMatrix[0].end(), 0);
  //	templateMatrix.resize(1, vector<int>(vars.size(), 1));
}

void InputData::addParamDirection(vector<double> d, double LB, double UB)
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
}

bool InputData::check()
{
  bool res = true;

  if (!isIterationSet()) {
    cerr << "Number of iterations is mandatory" << endl;
    res = false;
  }

  // each variable must have its dynamic
  for (unsigned i = 0; i < vars.size(); i++)
    if (!vars[i]->isDynamicDefined()) {
      cerr << "Variable " << vars[i]->getName() << " has not a dynamic"
           << endl;
      res = false;
    }

  // directions
  if (varMode != modeType::BOX && directions.size() == 0) {
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
  }

  // template
  if (varMode == modeType::POLY && templateMatrix.size() == 0) {
    cerr << "Template matrix must be provided if variable modality is "
         << varMode << endl;
    res = false;
  }

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

