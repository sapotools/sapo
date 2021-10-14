#ifndef __ABSSYN_H__
#define __ABSSYN_H__

#include<iostream>
#include <vector>
#include <string>
#include <utility>	// pair
#include <numeric>	// iota
#include <ginac/ginac.h>

// STL formulas
#include "Always.h"
#include "Eventually.h"
#include "Until.h"
#include "Conjunction.h"
#include "Disjunction.h"
#include "STL.h"
#include "Atom.h"

using namespace std;
using namespace GiNaC;

namespace AbsSyn
{

// forward declaration needed in Expr
class InputData;

// defines the types of problems solved by SAPO
enum problemType
{
	P_UNDEF,		// undefined
	REACH,			// reachability
	SYNTH				// parameter synthesis
};

// defines modalities for representing variables and parameters
enum modeType
{
	M_UNDEF,		// undefined
	BOX,				// boxes
	PARAL,			// parallelotopes
	POLY				// polytopes
};


// defines types of transformation
enum transType
{
	T_UNDEF,		// not yet defined
	AFO,				// All-For-One
	OFO					// One-For-One
};

/*
 ************************
 *      EXPRESSION      *
 ************************
 */

class Expr
{
	friend ostream& operator<<(ostream& os, Expr& e);

public:
	enum exprType
	{
		NUM_ATOM,			// atomic expression (single number)
		ID_ATOM,			// atomic expression (single identifier)
		SUM,					// sum (x + y)
		SUB,					// subtraction (x - y)
		MUL,					// multiplication (x * y)
		DIV,					// division (x / y)
		NEG						// unary minus (-x)
	};
	
	Expr(string n) {type = exprType::ID_ATOM; name = n; left = NULL; right = NULL; }	// atomic identifier
	Expr(double v) {type = exprType::NUM_ATOM; val = v; left = NULL; right = NULL; }	// atomic number
//	Expr(Expr *l, Expr *r, exprType t) { type = t; left = l; right = r; }							// SUM, SUB or MUL
//	Expr(Expr *e) { type = exprType::NEG; left = e; right = NULL; }										// NEG
	
	Expr *mul(Expr *e);
	Expr *div(Expr *e);
	Expr *sum(Expr *e);
	Expr *sub(Expr *e);
	Expr *neg();
	
	~Expr() { delete(left); delete(right); }
	
	exprType getType() { return type; }
	double getVal() { return val; }
	string getName() { return name; }
	Expr *getLeftOp() { return left; }
	Expr *getRightOp() { return right; }
	
	Expr *copy();			// deep copy of expression
	
	bool isNumeric(InputData *im);		// checks if the expression contains only numbers
	double evaluate(InputData *im);	// evaluates the value of a numeric expression
	
	ex toEx(InputData &m, const lst& vars, const lst& params);			// converts an Expr to a GiNaC ex
	
protected:
	Expr() { left = NULL; right = NULL; };
	
	exprType type;		// type of expression
	string name;			// name of ident (if ATOM)
	double val;				// value of number (if ATOM)
	Expr *left;				// left operand (if not ATOM)
	Expr *right;			// right operand (if not ATOM)
};



/*
 ***********************
 *       FORMULA       *
 ***********************
 */

class Formula
{
	friend ostream& operator<<(ostream& os, Formula& f);

public:
	enum formulaType
	{
		ATOM,					// atomic formula (ex <= 0)
		CONJ,					// conjunction (f1 && f2)
		DISJ,					// disjunction (f1 || f2)
		F_NEG,				// negation (! f1)
		ALW,					// always (G[i] f1)
		EVENT,				// eventually (F[i] f1)
		UNTIL					// until (f1 U[i] f2)
	};
	
	Formula(Expr *e) { type = formulaType::ATOM; ex = e; f1 = NULL; f2 = NULL; }							// atomic formula
//	Formula(Formula *l, Formula *r, formulaType t) { type = t; ex = NULL; f1 = l; f2 = r; }		// state formula (CONJ, DISJ, F_NEG
//	Formula(Formula *l, Formula *r, formulaType t, pair<int,int> in) { type = t; ex = NULL; f1 = l; f2 = r; i = in; }		// path formula (other types)
	
	Formula *conj(Formula *f);
	Formula *disj(Formula *f);
	Formula *neg();
	Formula *always(pair<int,int> in);
	Formula *eventually(pair<int,int> in);
	Formula *until(pair<int,int> in, Formula *f);
	
	~Formula() { delete(ex); delete(f1); delete(f2); }
	
	Expr *getEx() { return ex; }
	Formula *getLeft() { return f1; }
	Formula *getRight() { return f2; }
	pair<int, int> getInterval() { return i; }
	
	/*
	 * simplification, removes negations
	 */
	bool simplify();
	
	std::shared_ptr<STL> toSTL(InputData& m, const lst& vars, const lst& params);		// transforms a Formula into a SAPO STL formula
	
protected:
	Formula() { ex = NULL; f1 = NULL; f2 = NULL; }
	
	formulaType type;			// type of formula
	Expr *ex;							// expression of atomic formula (ex >= 0)
	Formula *f1;					// left formula
	Formula *f2;					// right formula
	pair<int, int> i;			// interval of application
	
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
	friend ostream& operator<<(ostream& os, Variable& v) { return os << v.name; }

public:
	Variable(string n) { name = n; dynamic = NULL; }
	~Variable() { delete(dynamic); }
	
	string getName() { return name; }
	Expr *getDynamic() { return dynamic; }
	
	void setDynamic(Expr *e) { dynamic = e; }
	
	// checks if dynamic has already been set
	bool isDynamicDefined() { return dynamic != NULL; }
	
protected:
	string name;					// name of the variable
	Expr *dynamic;				// dynamic associated with variable
};



/*
 *************************
 *       PARAMETER       *
 *************************
 */

class Parameter
{
	friend ostream& operator<<(ostream& os, Parameter& p) { return os << p.name; }

public:
	Parameter(string n) { name = n; }
	~Parameter() {}
	
	string getName() { return name; }
	
protected:
	string name;				// name of the parameter
};



/*
 ************************
 *       CONSTANT       *
 ************************
 */

class Constant
{
	friend ostream& operator<<(ostream& os, Constant& c) { return os << c.name; }

public:
	Constant(string n, double v) { name = n; val = v; }
	~Constant() {}
	
	string getName() { return name; }
	double getValue() { return val; }
	
protected:
	string name;				// name of the constant
	double val;					// value
};

class Definition
{
	friend ostream& operator<<(ostream& os, Definition& d) { return os << d.name; }
	
public:
	Definition(string id, Expr *e) { name = id; value = e; }
	~Definition() {}
	
	string getName() { return name; }
	Expr *getValue() { return value; }
	
private:
	string name;
	Expr *value;
};


/*
 ***********************
 *        MODEL        *
 ***********************
 */

class InputData
{
	friend ostream& operator<<(ostream& os, InputData& m);

public:
	InputData();
	
	~InputData() {}
	
	bool isProblemDefined() { return problem != problemType::P_UNDEF; }
	problemType getProblem() { return problem; }
	void setProblem(problemType t) { problem = t; }
	
	bool isVarModeDefined() { return varMode != modeType::M_UNDEF; }
	modeType getVarMode() { return varMode; }
	void setVarMode(modeType t) { varMode = t; }
	
	bool isParamModeDefined() { return paramMode != modeType::M_UNDEF; }
	modeType getParamMode() { return paramMode; }
	void setParamMode(modeType t) { paramMode = t; }
	
	int getIterations() { return iterations; }
	void setIterations(int n) { iterations = n; }
	
	unsigned getVarNum() { return vars.size(); }
	unsigned getParamNum() { return params.size(); }
	unsigned getConstNum() { return consts.size(); }
	unsigned getDefNumber() { return defs.size(); }
	
	bool isVarDefined(string name);			// checks if a variable named 'name' already exists
	bool isParamDefined(string name);		// checks if a parameter named 'name' already exists
	bool isConstDefined(string name);		// checks if a constant named 'name' already exists
	bool isDefDefined(string name);			// checks if a definition named 'name' already exists
	bool isSymbolDefined(string name);	// checks if a symbol (var, param, const or def) named 'name' already exists
	
	Variable *getVar(int i) { return vars[i]; }
	Variable *getVar(string name);				// return the variable named 'name', which must exist
	int getVarPos(string name);		// return an index such that vars[i] has name 'name'
	Parameter *getParam(int i) { return params[i]; }
	Parameter *getParam(string name);			// return the parameter named 'name', which must exist
	int getParamPos(string name);	// return an index i such that params[i] has name 'name'
	Constant *getConst(int i) { return consts[i]; }
	Constant *getConst(string name);			// return the constant named 'name', which must exist
	Definition *getDef(int i) { return defs[i]; }
	Definition *getDef(string name);			// return the definition named 'name', which must exist
	int getDefPos(string name);	// return an index i such that defs[i] has name 'name'
	
	void addVariable(Variable *v) { vars.push_back(v); }				// adds a new variable, which name is not already used
	void addParameter(Parameter *p) { params.push_back(p); }		// adds a new parameter, which name is not already used
	void addConstant(Constant *c) { consts.push_back(c); }			// adds a new constant, which name is not already used
	void addDefinition(Definition *d) { defs.push_back(d); }			// adds a new definition, which name is not already used
	
	bool isSpecDefined() { return spec != NULL; }
	void addSpec(Formula *f) { if (spec == NULL) spec = f; else spec = spec->conj(f); }
	Formula *getSpec() { return spec; }
	
	unsigned directionsNum() { return directions.size(); }
	vector<vector<double>> getDirections() { return directions; }
	void addDirection(vector<double> d, double LB, double UB);
	void addDirection(vector<double> d) { directions.push_back(d); }
	void addBounds(double LB, double UB) { LBoffsets.push_back(LB); UBoffsets.push_back(UB); }
	void defaultDirections();
	vector<double> getLB() { return LBoffsets; }
	vector<double> getUB() { return UBoffsets; }
	
	unsigned templateRows() { return templateMatrix.size(); }
	unsigned templateCols() { return templateMatrix[0].size(); }
	void setTemplate(vector<vector<int>> m) { templateMatrix = m; }
	void defaultTemplate();
	vector<vector<int>> getTemplate() { return templateMatrix; }
	
	unsigned paramDirectionsNum() { return paramDirections.size(); }
	void addParamDirection(vector<double> d, double LB, double UB);
	void addParamDirection(vector<double> d) { paramDirections.push_back(d); }
	void addParamBounds(double LB, double UB) { paramLBoffsets.push_back(LB); paramUBoffsets.push_back(UB); }
	void defaultParamDirections();
	vector<vector<double>> getParameterDirections() { return paramDirections; }
	vector<double> getParamDirection(int i) { return paramDirections[i]; }
	vector<double> getParamLB() { return paramLBoffsets; }
	vector<double> getParamUB() { return paramUBoffsets; }
	
	bool isTransModeDefined() { return trans != transType::T_UNDEF; }
	void setTransMode(transType t) { trans = t; }
	transType getTransMode() { return trans; }
	int getTransValue();			// returns int value used by sapo (AFO = 1, OFO = 0)
	
	bool isDecompositionDefined() { return decomp_defined; }
	void setDecomposition() { decomp = true; }
	bool getDecomposition() { return decomp; }
	
	bool isAlphaDefined() { return alpha >= 0; }
	void setAlpha(double a) { alpha = a; }
	double getAlpha() { return alpha; }
	
	bool check();				// checks for errors in model
	
protected:
	problemType problem;
	
	modeType varMode;
	modeType paramMode;
	
	int iterations;
	
	vector<Variable *> vars;
	vector<Parameter *> params;
	vector<Constant *> consts;
	vector<Definition *> defs;
	
	Formula *spec;
	
	vector<vector<double>> directions;
	vector<double> LBoffsets;
	vector<double> UBoffsets;
	
	vector<vector<int>> templateMatrix;
	
	vector<vector<double>> paramDirections;
	vector<double> paramLBoffsets;
	vector<double> paramUBoffsets;
	
	// SAPO options
	transType trans;
	bool decomp, decomp_defined;
	double alpha;
};

}

namespace std
{
	ostream& operator<<(ostream& os, const AbsSyn::problemType t);
	ostream& operator<<(ostream& os, const AbsSyn::modeType t);

	ostream& operator<<(ostream& os, const pair<int,int> p);
	ostream& operator<<(ostream& os, const pair<double,double> p);

	ostream& operator<<(ostream& os, const vector<double> v);
	ostream& operator<<(ostream& os, const vector<string> v);
	ostream& operator<<(ostream& os, const vector<int> v);
	ostream& operator<<(ostream& os, const vector<vector<double>> v);
	ostream& operator<<(ostream& os, const vector<vector<int>> v);
}
#endif
