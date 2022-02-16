#include "InputData.h"

using namespace std;

namespace AbsSyn
{

InputData::InputData():
    problem(problemType::P_UNDEF), varMode(modeType::M_UNDEF),
    paramMode(modeType::M_UNDEF), iterations(0), iter_set(false),
    max_param_splits(0), presplits(false),
    max_bundle_magnitude(std::numeric_limits<double>::max()), vars(), params(),
    consts(), defs(), assumptions(), spec(NULL), directions(),
    templateMatrix(), paramDirections(),
    trans(transType::T_UNDEF), decomp(false), decomp_defined(false),
    alpha(0.5), alphaDefined(false)
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

  spec.reset();
}
ostream &operator<<(ostream &os, const InputData &m)
{
  os << "Problem: " << m.problem << endl;
  os << "Iterations: " << m.iterations << endl;

  os << endl;
  os << "Variables: " << endl;
  for (unsigned i = 0; i < m.vars.size(); i++) {
    os << "\t" << m.vars[i]->getName() << ": " << m.vars[i]->getDynamic()
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
    os << "\t" << m.defs[i]->getName() << " = " << m.defs[i]->getValue()
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

int find(std::vector<Direction *> M, Direction *v, const Context &ctx, bool variable = true)
{
	unsigned pos = 0;
	while (pos < M.size()) {
		if (M[pos]->compare(v, ctx, variable)) {
			return pos;
		} else {
			pos++;
		}
	}
	// not found
	return -1;
}

void InputData::addDirectionConstraint(Direction *d, Context ctx)
{
	if (d->getType() == Direction::Type::INT) {		// constraint is an interval specification
		this->addDirectionConstraint(new Direction(d->getLHS(), d->getUB(ctx), Direction::Type::LE, 0, 0, d->getName()), ctx);
		this->addDirectionConstraint(new Direction(d->getLHS(), d->getLB(ctx), Direction::Type::GE, 0, 0, d->getName()), ctx);
		delete(d);
		return;
	} else if (d->getType() == Direction::Type::EQ) {
		this->addDirectionConstraint(new Direction(d->getLHS(), d->getUB(ctx), Direction::Type::LE, 0, 0, d->getName()), ctx);
		this->addDirectionConstraint(new Direction(d->getLHS(), -d->getLB(ctx), Direction::Type::GE, 0, 0, d->getName()), ctx);
		delete(d);
		return;
	}
	
	Direction *new_dir = d;
	Direction *negated_dir = new_dir->getComplementary();
	
	// check if new direction is already present
	int pos = find(directions, new_dir, ctx);
	int negated_pos = find(directions, negated_dir, ctx);
	
	if (pos != -1) {
		
		// direction already present
		if (!directions[pos]->hasUB() || new_dir->getUB(ctx) < directions[pos]->getUB(ctx)) {
			directions[pos]->setUB(ctx, new_dir->getUB(ctx));
			directions[pos]->setName(new_dir->getName());
		}
		if (!directions[pos]->hasLB() || new_dir->getLB(ctx) > directions[pos]->getLB(ctx)) {
			directions[pos]->setLB(ctx, new_dir->getLB(ctx));
			directions[pos]->setName(new_dir->getName());
		}
		delete(new_dir);
		delete(negated_dir);
		
	} else if (negated_pos != -1) {
		
		// negation of direction already present
		if (!directions[negated_pos]->hasUB() || negated_dir->getUB(ctx) < directions[negated_pos]->getUB(ctx)) {
			directions[negated_pos]->setUB(ctx, negated_dir->getUB(ctx));
			directions[negated_pos]->setName(negated_dir->getName());
		}
		if (!directions[negated_pos]->hasLB() || negated_dir->getLB(ctx) > directions[negated_pos]->getLB(ctx)) {
			directions[negated_pos]->setLB(ctx, negated_dir->getLB(ctx));
			directions[negated_pos]->setName(negated_dir->getName());
		}
		delete(new_dir);
		delete(negated_dir);
		
	} else {
		
		// new direction
		directions.push_back(new_dir);
		
		// cover variables
		for (unsigned i = 0; i < vars.size(); i++) {
			if (new_dir->covers(ctx, vars[i]->getName())) {
				vars[i]->setCovered();
			}
		}
		
		delete(negated_dir);
	}
}


int InputData::findDirectionPos(const std::string &name) const
{
	for (unsigned i = 0; i < directions.size(); i++) {
		if (directions[i]->getName() == name) {
			return i;
		}
	}
	return -1;
}


void InputData::addParamDirectionConstraint(Direction *d, Context ctx)
{
//	std::cout << "Direction " << *d << std::endl;
	
	if (d->getType() == Direction::Type::INT) {		// constraint is an interval specification
		this->addParamDirectionConstraint(new Direction(d->getLHS(), d->getUB(ctx), Direction::Type::LE, 0, 0, d->getName()), ctx);
		this->addParamDirectionConstraint(new Direction(d->getLHS(), d->getLB(ctx), Direction::Type::GE, 0, 0, d->getName()), ctx);
		delete(d);
		return;
	}
	
	Direction *new_dir = d;
	Direction *negated_dir = new_dir->getComplementary();
	
	// check if new direction is already present
	int pos = find(paramDirections, new_dir, ctx, false);
	int negated_pos = find(paramDirections, negated_dir, ctx, false);
	
	if (pos != -1) {
		// direction already present
		if (!paramDirections[pos]->hasUB() || new_dir->getUB(ctx) < paramDirections[pos]->getUB(ctx)) {
			paramDirections[pos]->setUB(ctx, new_dir->getUB(ctx));
			paramDirections[pos]->setName(new_dir->getName());
		}
		if (!paramDirections[pos]->hasLB() || new_dir->getLB(ctx) > paramDirections[pos]->getLB(ctx)) {
			paramDirections[pos]->setLB(ctx, new_dir->getLB(ctx));
			paramDirections[pos]->setName(new_dir->getName());
		}
		delete(new_dir);
		delete(negated_dir);
		
	} else if (negated_pos != -1) {
		// negation of direction already present		
		if (!paramDirections[negated_pos]->hasLB() || negated_dir->getLB(ctx) > paramDirections[negated_pos]->getLB(ctx)) {
			paramDirections[negated_pos]->setLB(ctx, negated_dir->getLB(ctx));
			paramDirections[negated_pos]->setName(negated_dir->getName());
		}
		if (!paramDirections[negated_pos]->hasUB() || negated_dir->getUB(ctx) < paramDirections[pos]->getUB(ctx)) {
			paramDirections[negated_pos]->setUB(ctx, negated_dir->getUB(ctx));
			paramDirections[negated_pos]->setName(negated_dir->getName());
		}
		delete(new_dir);
		delete(negated_dir);
		
	} else {
		// new direction
		paramDirections.push_back(new_dir);
		
		// cover parameters
		for (unsigned i = 0; i < params.size(); i++) {
			if (new_dir->covers(ctx, params[i]->getName())) {
				params[i]->setCovered();
			}
		}
	}
}


bool InputData::check(Context ctx)
{
  bool res = true;

	if (!isProblemDefined()) {
		cerr << "Problem type must be defined" << endl;
		res = false;
	}
	
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
	/*for (unsigned i = 0; i < vars.size(); i++) {
    if (!vars[i]->isCovered()) {
      cerr << "Variable " << vars[i]->getName() << " is not covered by any direction"
           << endl;
      res = false;
    }
	}*/
	
	// for each variable, check if it appears  positively (or negatively) in any constraint.
	// If so, we conclude that it is upper (or lower) bounded
	{	// put a block to avoid namespace pollution
		std::vector<bool> hasLB(vars.size(), false), hasUB(vars.size(), false);
		for (unsigned d = 0; d < directions.size(); d++) {
			std::vector<double> dirVector = directions[d]->getDirectionVector(ctx, true);
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
			}
			if (!hasUB[v]) {
				cerr << "Variable " << vars[v]->getName() << " has no finite upper bound" << endl;
				res = false;
			}
		}
	}
	
	
	// set variable directions LB where needed
	if (res) {
		vector<vector<double>> A{};
		vector<double> b{};
		for (unsigned i = 0; i < directions.size(); i++) {
			// add only bounded directions
			if (directions[i]->hasUB()) {
				A.push_back(directions[i]->getDirectionVector(ctx, true));
				b.push_back(directions[i]->getUB(ctx));
			}
		}
		for (unsigned i = 0; i < directions.size(); i++) {
			if (directions[i]->hasLB()) {
				A.push_back(-directions[i]->getDirectionVector(ctx, true));
				b.push_back(-directions[i]->getLB(ctx));
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
				std::vector<double> dirVector = directions[i]->getDirectionVector(ctx, true);
				SymbolicAlgebra::Expression<> obj_function = 0;
				for (unsigned j = 0; j < dirVector.size(); j++) {
					obj_function += dirVector[j] * symbols[j];
				}
				
				double minVal = LS.minimize(symbols, obj_function);
				directions[i]->setLB(ctx, minVal);
			}
			
			if (!directions[i]->hasUB()) {
				std::vector<double> dirVector = directions[i]->getDirectionVector(ctx, true);
				SymbolicAlgebra::Expression<> obj_function = 0;
				for (unsigned j = 0; j < dirVector.size(); j++) {
					obj_function += dirVector[j] * symbols[j];
				}
				
				double maxVal = LS.maximize(symbols, obj_function);
				directions[i]->setUB(ctx, maxVal);
			}
		}
	}
	
	// each template row must be bounded
	{
		for (unsigned i = 0; i < templateMatrix.size(); i++) {
			vector<vector<double>> M{};
			for (unsigned j = 0; j < templateMatrix[i].size(); j++) {
				M.push_back(directions[templateMatrix[i][j]]->getDirectionVector(ctx, true));
			}
			
			vector<double> zeroes(templateMatrix[i].size(), 0);
			
			DenseLinearAlgebra::PLU_Factorization<double> PLU(M);
			try {
				vector<double> res = PLU.solve(zeroes);
			} catch (domain_error &e) {
				// directions are dependent, parallelotope is not bounded
				cerr << "Template row " << templateMatrix[i] << " defines an unbounded parallelotope" << endl;
				res = false;
			}
		}
	}
	
	// each param must be covered
	/*for (unsigned i = 0; i < params.size(); i++) {
    if (!params[i]->isCovered()) {
      cerr << "Parameter " << params[i]->getName() << " is not covered by any direction"
           << endl;
      res = false;
    }
	}*/
	
	// for each parameter, check if it appears  positively (or negatively) in any constraint.
	// If so, we conclude that it is upper (or lower) bounded
	{	// put a block to avoid namespace pollution
//		double infinity = std::numeric_limits<double>::max();
//		double negInfinity = -std::numeric_limits<double>::infinity();
		
		std::vector<bool> hasLB(params.size(), false), hasUB(params.size(), false);
		for (unsigned d = 0; d < paramDirections.size(); d++) {
			std::vector<double> dirVector = paramDirections[d]->getDirectionVector(ctx, false);
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
			}
			if (!hasUB[p]) {
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
				A.push_back(paramDirections[i]->getDirectionVector(ctx, false));
				b.push_back(paramDirections[i]->getUB(ctx));
			}
			if (paramDirections[i]->hasLB()) {
				A.push_back(paramDirections[i]->getComplementary()->getDirectionVector(ctx, false));
				b.push_back(paramDirections[i]->getLB(ctx));
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
				std::vector<double> dirVector = paramDirections[i]->getDirectionVector(ctx, false);
				for (unsigned j = 0; j < dirVector.size(); j++) {
					obj_function += dirVector[j] * symbols[j];
				}
				
				double minVal = LS.minimize(symbols, obj_function);
				paramDirections[i]->setLB(ctx, minVal);
			}
			
			if (!paramDirections[i]->hasUB()) {
				SymbolicAlgebra::Expression<> obj_function = 0;
				std::vector<double> dirVector = paramDirections[i]->getDirectionVector(ctx, false);
				for (unsigned j = 0; j < dirVector.size(); j++) {
					obj_function += dirVector[j] * symbols[j];
				}
				
				double maxVal = LS.maximize(symbols, obj_function);
				paramDirections[i]->setUB(ctx, maxVal);
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
  
  // trans type (AFO or OFO), if undefined set to AFO
  if (trans == transType::T_UNDEF) {
		trans = transType::AFO;
	}

  return res;
}

} // end namespace AbsSyn
