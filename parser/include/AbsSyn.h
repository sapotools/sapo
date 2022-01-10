#ifndef __ABSSYN_H__
#define __ABSSYN_H__

#include <iostream>
#include <numeric> // iota
#include <string>
#include <utility> // pair
#include <vector>

// STL formulas
#include "Always.h"
#include "Atom.h"
#include "Conjunction.h"
#include "Disjunction.h"
#include "Eventually.h"
#include "STL.h"
#include "Until.h"

#include "SymbolicAlgebra.h"

namespace AbsSyn
{

// forward declaration needed in Expr
class InputData;

// defines the types of problems solved by SAPO
enum problemType {
  P_UNDEF, // undefined
  REACH,   // reachability
  SYNTH    // parameter synthesis
};

// defines modalities for representing variables and parameters
enum modeType {
  M_UNDEF, // undefined
  BOX,     // boxes
  PARAL,   // parallelotopes
  POLY     // polytopes
};

// defines types of transformation
enum transType {
  T_UNDEF, // not yet defined
  AFO,     // All-For-One
  OFO      // One-For-One
};

/*
 ************************
 *      EXPRESSION      *
 ************************
 */

class Expr
{
  friend std::ostream &operator<<(std::ostream &os, const Expr &e)
	{
		return e.prettyPrint(os, 10);
	}

public:
	// order of types follows order of priority (needed by Expr::prettyPrint)
  enum exprType {
    NUM_ATOM = 0, // atomic expression (single number)
    ID_ATOM,  // atomic expression (single identifier)
    NEG,       // unary minus (-x)
    DIV,      // division (x / y)
    MUL,      // multiplication (x * y)
    SUB,      // subtraction (x - y)
    SUM      // sum (x + y)
  };

  Expr(std::string n)
  {
    type = exprType::ID_ATOM;
    name = n;
    left = NULL;
    right = NULL;
  } // atomic identifier
  Expr(double v)
  {
    type = exprType::NUM_ATOM;
    val = v;
    left = NULL;
    right = NULL;
  } // atomic number
    //	Expr(Expr *l, Expr *r, exprType t) { type = t; left = l; right = r; }
  //// SUM, SUB or MUL 	Expr(Expr *e) { type = exprType::NEG; left = e; right =
  // NULL; }
  // // NEG

  Expr *mul(Expr *e);
  Expr *div(Expr *e);
  Expr *sum(Expr *e);
  Expr *sub(Expr *e);
  Expr *neg();

  ~Expr()
  {
    delete (left);
    delete (right);
  }

  exprType getType()
  {
    return type;
  }

  double getVal()
  {
    return val;
  }

  std::string getName()
  {
    return name;
  }

  Expr *getLeftOp()
  {
    return left;
  }

  Expr *getRightOp()
  {
    return right;
  }
  
  int getDegree(const InputData &id)
			const; // return the degree of the polynomial expression considering only vars

  Expr *copy() const; // deep copy of expression

  bool isNumeric(const InputData &im)
      const; // checks if the expression contains only numbers
  double evaluate(const InputData &im)
      const; // evaluates the value of a numeric expression
	
	bool hasParams(const InputData &im)
			const; // checks if the expression contains parameter names
      
	double getCoefficient(const InputData &id, const std::string name)
			const;	 /*
								* returns the coefficient of the variable "name"
								* in the simplified expression (in which each variable
								* apperas only one time)
								* Applies only to linear expressions without parameters
								*/
	
	double getOffset(const InputData &id)
			const;	 /*
								* returns the numerical term
								* in the simplified expression (in which each variable
								* apperas only one time)
								* Applies only to linear expressions without parameters
								*/

  SymbolicAlgebra::Expression<>
  toEx(const InputData &m, const std::vector<SymbolicAlgebra::Symbol<>> &vars,
       const std::vector<SymbolicAlgebra::Symbol<>> &params) const;

protected:
  Expr()
  {
    left = NULL;
    right = NULL;
  };

  exprType type;    // type of expression
  std::string name; // name of ident (if ATOM)
  double val;       // value of number (if ATOM)
  Expr *left;       // left operand (if not ATOM)
  Expr *right;      // right operand (if not ATOM)
  
  // utility for printing expressions without too many parenthesis
  std::ostream &prettyPrint(std::ostream &os, const int level) const;
};

/*
 ***********************
 *       FORMULA       *
 ***********************
 */

class Formula
{
  friend std::ostream &operator<<(std::ostream &os, const Formula &f);

public:
  enum formulaType {
    ATOM,  // atomic formula (ex <= 0)
    CONJ,  // conjunction (f1 && f2)
    DISJ,  // disjunction (f1 || f2)
    F_NEG, // negation (! f1)
    ALW,   // always (G[i] f1)
    EVENT, // eventually (F[i] f1)
    UNTIL  // until (f1 U[i] f2)
  };

  Formula(Expr *e)
  {
    type = formulaType::ATOM;
    ex = e;
    f1 = NULL;
    f2 = NULL;
  } // atomic formula
    //	Formula(Formula *l, Formula *r, formulaType t) { type = t; ex = NULL;
    // f1 = l; f2 = r; }		// state formula (CONJ, DISJ, F_NEG
    //	Formula(Formula *l, Formula *r, formulaType t, pair<int,int> in) { type
    //= t; ex = NULL; f1 = l; f2 = r; i = in; }		// path formula (other
    // types)

  Formula *conj(Formula *f);
  Formula *disj(Formula *f);
  Formula *neg();
  Formula *always(std::pair<int, int> in);
  Formula *eventually(std::pair<int, int> in);
  Formula *until(std::pair<int, int> in, Formula *f);

  ~Formula()
  {
    delete (ex);
    delete (f1);
    delete (f2);
  }

  Expr *getEx()
  {
    return ex;
  }

  Formula *getLeft()
  {
    return f1;
  }

  Formula *getRight()
  {
    return f2;
  }

  std::pair<int, int> getInterval()
  {
    return i;
  }
  
  formulaType getType() const
	{
		return type;
	}
  
  bool isLinear(const InputData &id)
			const; // checks if the formula is a boolean combination
						 // of linear inequalities

  /*
   * simplification, removes negations
   */
  bool simplify();

  /**
   * @brief Transforms a Formula into a SAPO STL formula.
   *
   * @param m is the input data object.
   * @param vars is the vector of variable symbols.
   * @param params is the vector of parameter symbols.
   * @return a pointer to the STL formula.
   */
  std::shared_ptr<STL>
  toSTL(const InputData &m, const std::vector<SymbolicAlgebra::Symbol<>> &vars,
        const std::vector<SymbolicAlgebra::Symbol<>> &params) const;

protected:
  Formula()
  {
    ex = NULL;
    f1 = NULL;
    f2 = NULL;
  }

  formulaType type;      // type of formula
  Expr *ex;              // expression of atomic formula (ex >= 0)
  Formula *f1;           // left formula
  Formula *f2;           // right formula
  std::pair<int, int> i; // interval of application

  /*
   * utility for semplification
   * returns 0 if no further calls are needed
   * returns 1 if more are required
   * returns 2 if a negation of an until is found
   */
  int simplifyRec();
};

/*
 ************************
 *       VARIABLE       *
 ************************
 */

class Variable
{
  friend std::ostream &operator<<(std::ostream &os, Variable &v)
  {
    return os << v.name;
  }

public:
  Variable(const std::string &n): name(n), dynamic(NULL) {}

  ~Variable()
  {
    delete (dynamic);
  }

  const std::string &getName() const
  {
    return name;
  }

  const Expr *getDynamic() const
  {
    return dynamic;
  }

  void setDynamic(Expr *e)
  {
    dynamic = e;
  }

  // checks if dynamic has already been set
  bool isDynamicDefined()
  {
    return dynamic != NULL;
  }

protected:
  std::string name; // name of the variable
  Expr *dynamic;    // dynamic associated with variable
};

/*
 *************************
 *       PARAMETER       *
 *************************
 */

class Parameter
{
  friend std::ostream &operator<<(std::ostream &os, const Parameter &p)
  {
    return os << p.name;
  }

public:
  Parameter(const std::string &n): name(n) {}

  ~Parameter() {}

  const std::string &getName() const
  {
    return name;
  }

protected:
  std::string name; // name of the parameter
};

/*
 ************************
 *       CONSTANT       *
 ************************
 */

class Constant
{
  friend std::ostream &operator<<(std::ostream &os, Constant &c)
  {
    return os << c.name;
  }

public:
  Constant(std::string n, double v): name(n), val(v) {}

  ~Constant() {}

  std::string getName() const
  {
    return name;
  }
  double getValue() const
  {
    return val;
  }

protected:
  std::string name; // name of the constant
  double val;       // value
};


/*
 ************************
 *      DEFINITION      *
 ************************
 */

class Definition
{
  friend std::ostream &operator<<(std::ostream &os, const Definition &d)
  {
    return os << d.name;
  }

public:
  Definition(std::string id, Expr *e): name(id), value(e) {}

  ~Definition()
  {
    delete value;
  }

  std::string getName()
  {
    return name;
  }
  Expr *getValue()
  {
    return value;
  }
  const Expr *getValue() const
  {
    return value;
  }

private:
  std::string name;
  Expr *value;
};

/*
 *************************
 *       ASSERTION       *
 *************************
 */

class Assumption
{
  friend std::ostream &operator<<(std::ostream &os, const Assumption &a)
	{
		return os << *(a.ex);
	}

public:
	
	Assumption(Expr *e): ex(e) {}
	
	~Assumption()
	{
		delete ex;
	}
	
	Expr *getExpr() const
	{
		return ex;
	}
	
	std::vector<double> getDirection(const InputData &id)
				const;		// return the direction corresponding to the linear constraint
	
	double getOffset(const InputData &id) const		// return the offset of the constraint*/
	{
		return -ex->getOffset(id);
	}
	
protected:
	Expr *ex;		// constraint of the form ex <= 0
};

/*
 ***********************
 *        MODEL        *
 ***********************
 */

class InputData
{
  friend std::ostream &operator<<(std::ostream &os, const InputData &m);

public:
  InputData();

  ~InputData();

  bool isProblemDefined() const
  {
    return problem != problemType::P_UNDEF;
  }
  const problemType &getProblem() const
  {
    return problem;
  }
  void setProblem(problemType t)
  {
    problem = t;
  }

  bool isVarModeDefined() const
  {
    return varMode != modeType::M_UNDEF;
  }
  const modeType &getVarMode() const
  {
    return varMode;
  }
  void setVarMode(modeType t)
  {
    varMode = t;
  }

  bool isParamModeDefined() const
  {
    return paramMode != modeType::M_UNDEF;
  }
  const modeType &getParamMode() const
  {
    return paramMode;
  }
  void setParamMode(modeType t)
  {
    paramMode = t;
  }

  bool isIterationSet() const
  {
    return iter_set;
  }

  const unsigned int &getIterations() const
  {
    return iterations;
  }

  void setIterations(unsigned int n)
  {
    iter_set = true;
    iterations = n;
  }

  const unsigned int &getMaxParameterSplits() const
  {
    return max_param_splits;
  }

  void setMaxParameterSplits(unsigned int n)
  {
    max_param_splits = n;
  }

  const bool &isPreSplitsSet() const
  {
    return presplits;
  }

  void setPreSplits(const bool presplits)
  {
    this->presplits = presplits;
  }

  unsigned getVarNum() const
  {
    return vars.size();
  }
  unsigned getParamNum() const
  {
    return params.size();
  }
  unsigned getConstNum() const
  {
    return consts.size();
  }
  unsigned getDefNumber() const
  {
    return defs.size();
  }
  unsigned getAssumptionsNumber() const
  {
		return assumptions.size();
	}

  bool isVarDefined(const std::string &name)
      const; // checks if a variable named 'name' already exists
  bool isParamDefined(const std::string &name)
      const; // checks if a parameter named 'name' already exists
  bool isConstDefined(const std::string &name)
      const; // checks if a constant named 'name' already exists
  bool isDefDefined(const std::string &name)
      const; // checks if a definition named 'name' already exists
  bool isSymbolDefined(
      const std::string &name) const; // checks if a symbol (var, param, const
                                      // or def) named 'name' already exists

  const Variable *getVar(int i) const
  {
    return vars[i];
  }
  const Variable *getVar(const std::string &name)
      const; // return the variable named 'name', which must exist
  Variable *getVar(int i)
  {
    return vars[i];
  }
  Variable *getVar(const std::string &name);
  int getVarPos(const std::string &name)
      const; // return an index such that vars[i] has name 'name'
  const Parameter *getParam(int i) const
  {
    return params[i];
  }
  const Parameter *getParam(const std::string &name)
      const; // return the parameter named 'name', which must exist
  int getParamPos(const std::string &name)
      const; // return an index i such that params[i] has name 'name'
  const Constant *getConst(int i) const
  {
    return consts[i];
  }
  const Constant *getConst(const std::string &name)
      const; // return the constant named 'name', which must exist
  const Definition *getDef(int i) const
  {
    return defs[i];
  }
  const Definition *getDef(const std::string &name)
      const; // return the definition named 'name', which must exist
  int getDefPos(const std::string &name)
      const; // return an index i such that defs[i] has name 'name'
	
	const Assumption *getAssumption(int i) const // return the assumption in position i
	{
		return assumptions[i];
	}

  void addVariable(Variable *v)
  {
    vars.push_back(v);
  } // adds a new variable, which name is not already used
  void addParameter(Parameter *p)
  {
    params.push_back(p);
  } // adds a new parameter, which name is not already used
  void addConstant(Constant *c)
  {
    consts.push_back(c);
  } // adds a new constant, which name is not already used
  void addDefinition(Definition *d)
  {
    defs.push_back(d);
  } // adds a new definition, which name is not already used
  void addAssumption(Assumption *a)
	{
		assumptions.push_back(a);
	} // adds a new assumption

  bool isSpecDefined() const
  {
    return spec != NULL;
  }
  void addSpec(Formula *f)
  {
    if (spec == NULL)
      spec = f;
    else
      spec = spec->conj(f);
  }
  const Formula *getSpec() const
  {
    return spec;
  }

  unsigned int directionsNum() const
  {
    return directions.size();
  }
  const std::vector<std::vector<double>> &getDirections() const
  {
    return directions;
  }
  void addDirection(std::vector<double> d, double LB, double UB);
  void addDirection(std::vector<double> d)
  {
    directions.push_back(d);
  }
  void addBounds(double LB, double UB)
  {
    LBoffsets.push_back(LB);
    UBoffsets.push_back(UB);
  }
  void defaultDirections();
  const std::vector<double> &getLB() const
  {
    return LBoffsets;
  }
  const std::vector<double> &getUB() const
  {
    return UBoffsets;
  }

  unsigned templateRows() const
  {
    return templateMatrix.size();
  }
  unsigned templateCols() const
  {
    return templateMatrix[0].size();
  }
  void setTemplate(std::vector<std::vector<int>> m)
  {
    templateMatrix = m;
  }
  void defaultTemplate();
  const std::vector<std::vector<int>> &getTemplate() const
  {
    return templateMatrix;
  }

  unsigned paramDirectionsNum() const
  {
    return paramDirections.size();
  }
  void addParamDirection(std::vector<double> d, double LB, double UB);
  void addParamDirection(std::vector<double> d)
  {
    paramDirections.push_back(d);
  }
  void addParamBounds(double LB, double UB)
  {
    paramLBoffsets.push_back(LB);
    paramUBoffsets.push_back(UB);
  }

  void defaultParamDirections();

  const std::vector<std::vector<double>> &getParameterDirections() const
  {
    return paramDirections;
  }

  const std::vector<double> &getParamDirection(int i) const
  {
    return paramDirections[i];
  }

  const std::vector<double> &getParamLB() const
  {
    return paramLBoffsets;
  }

  const std::vector<double> &getParamUB() const
  {
    return paramUBoffsets;
  }

  bool isTransModeDefined() const
  {
    return trans != transType::T_UNDEF;
  }

  void setTransMode(transType t)
  {
    trans = t;
  }

  const transType &getTransMode() const
  {
    return trans;
  }

  unsigned char
  getTransValue() const // returns int value used by sapo (AFO = 1, OFO = 0)
  {
    if (trans == transType::AFO)
      return 1;
    else
      return 0;
  }

  const bool &isDecompositionDefined() const
  {
    return decomp_defined;
  }
  void setDecomposition()
  {
    decomp = true;
  }
  const bool &getDecomposition() const
  {
    return decomp;
  }

  bool isAlphaDefined() const
  {
    return alpha >= 0;
  }
  void setAlpha(double a)
  {
    alpha = a;
  }
  const double &getAlpha() const
  {
    return alpha;
  }

  bool check(); // checks for errors in model

protected:
  problemType problem;

  modeType varMode;
  modeType paramMode;

  unsigned int iterations;
  bool iter_set;

  unsigned int max_param_splits;
  bool presplits;

  std::vector<Variable *> vars;
  std::vector<Parameter *> params;
  std::vector<Constant *> consts;
  std::vector<Definition *> defs;
	
	std::vector<Assumption *> assumptions;

  Formula *spec;

  std::vector<std::vector<double>> directions;
  std::vector<double> LBoffsets;
  std::vector<double> UBoffsets;

  std::vector<std::vector<int>> templateMatrix;

  std::vector<std::vector<double>> paramDirections;
  std::vector<double> paramLBoffsets;
  std::vector<double> paramUBoffsets;

  // SAPO options
  transType trans;
  bool decomp, decomp_defined;
  double alpha;
};

}
#endif
