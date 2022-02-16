%skeleton "lalr1.cc" // -*- C++ -*-
%require "3.5.1"
%defines

%define api.token.constructor
%define api.value.type variant
%define parse.assert

%code requires {
	#include <string>
	#include <fstream>
	#include <cmath>
	
	#include "locations.h"
	
	#include "AbsSynIO.h"
	
	#include "Expr.h"
	
	#include "STL.h"
	#include "Atom.h"
	#include "Conjunction.h"
	#include "Disjunction.h"
	#include "Negation.h"
	#include "Always.h"
	#include "Eventually.h"
	#include "Until.h"
	
	#include "types.h"
	#include "Variable.h"
	#include "Parameter.h"
	#include "Constant.h"
	#include "Definition.h"
	#include "Direction.h"
	#include "InputData.h"

	class driver;
	
	inline void warning(const yy::location &l, const std::string &m)
	{
		std::cerr << "Warning at line " << l.end.line << ", column " << l.end.column << ": " << m << '\n';
	}
	
	inline void print(const yy::location &l, const std::string filename)
	{
		using namespace std;
		
		// open file
		ifstream file;
		if (filename.empty() || filename == "-") {
			// TODO: check if we can read from stdin
			return;
		} else {
			file.open(filename);
		}
		
		// skip to correct line
		string line;
		for (int i = 0; i < l.begin.line - 1; i++) {
			getline(file, line);
		}
		
		// print useful lines
		int digits = floor(log10(l.begin.line)) + 2;
		for (int i = l.begin.line+1; i <= l.end.line + 1; i++) {
			getline(file, line);
			int current_digits = floor(log10(i)) + 1;
			
			cerr << "  ";
			for (int j = 0; j < digits - current_digits; j++) {
				cerr << " ";
			}
			cerr << i << " | " << line << endl;
		}
		
		if (l.begin.line == l.end.line) {
			cerr << "  ";
			for (int i = 0; i < digits; i++) {
				cerr << " ";
			}
			cerr << " |";
			for (int i = 0; i < l.begin.column; i++) {
				cerr << " ";
			}
			for (int i = l.begin.column; i < l.end.column; i++) {
				cerr << "^";
			}
			cerr << endl;
		}
		
		// close file
		file.close();
		
		return;
	}
	
	// macro for errors
	#define ERROR(loc, msg)\
		yy::parser::error(loc, msg);\
		print(loc, drv.file);\
		drv.errors = true;
	
	// macro for warnings
	#define WARNING(loc, msg)\
		warning(loc, msg);\
		print(loc, drv.file);
}

// The parsing context.
%param { driver& drv }

%locations
%define api.location.file "../include/locations.h"

%define parse.trace
%define parse.error verbose

%code {
#include "driver.h"
}

%define api.token.prefix {TOK_}
%token
	END 0
	PROB 
	REACH
	SYNTH
	VARMODE
	PARAMMODE
	BOX
	PARAL
	POLY
	VAR
	PARAM
	CONST
	DEFINE
	IN
	DYN
	SPEC
	ASSUME
	ITER
	PSPLITS
	PRESPLITS
	MAX_MAGNITUDE
	DIR
	TEMPL
	PDIR
	OPT
	TRANS
	AFO
	OFO
	DECOMP
	ALPHA
	AND					"&&"
	OR					"||"
	NOT					"!"
	LE					"<="
	GE					">="
	EQ					"="
	ALWAYS			"G"
	EVENTUALLY	"F"
	UNTIL				"U"
	MINUS				"-"
	PLUS				"+"
	STAR				"*"
	POW					"^"
	DIV					"/"
	LPAREN			"("
	RPAREN			")"
	LSQUARE			"["
	RSQUARE			"]"
	LBRACE			"{"
	RBRACE			"}"
	LANGLE			"<"
	RANGLE			">"
	COLON				":"
	SEMICOLON		";"
	COMMA				","
;

%token <std::string> IDENT
%token <std::string> ON
%token <std::string> OFF
%token <int> INTEGER
%token <double> DOUBLE

%left "&&"
%left "||"
%left "!"
%left "F"
%left "G"
%left "U"

%left "+" "-"
%left "*" "/"
%left UMINUS
%right "^"

%nterm <AbsSyn::problemType> problemType
%nterm <AbsSyn::modeType> modeType
%nterm <std::pair<double, double>> doubleInterval
%nterm <std::pair<int, int>> intInterval
%nterm <double> number
%nterm <std::vector<std::string>> identList
%nterm <SymbolicAlgebra::Expression<>> expr
%nterm <std::vector<int>> matrixRow
%nterm <std::vector<std::vector<int>>> rowList
%nterm <std::shared_ptr<STL> > formula
%nterm <AbsSyn::transType> transType
%nterm <AbsSyn::Direction *> direction
%nterm <AbsSyn::Direction::Type> directionType

%printer { yyo << $$; } <*>;

%start s

%%

s		: statement
		| s statement
		| END
		{
			ERROR(@1, "Empty file");
			YYERROR;
		}

statement		: header {}
						| symbol {}
						| matrices {}
						| option {}

header			: PROB ":" problemType ";"
						{
							if (drv.data.isProblemDefined()) {
								ERROR(@4, "Problem has already been defined");
							} else {
								drv.data.setProblem($3);
							}
						}
						| PARAMMODE ":" modeType ";"
						{
							WARNING(@$, "Parameter modality is deprecated and will be ignored");
						}
						| VARMODE ":" modeType ";"
						{
							WARNING(@$, "Variable modality is deprecated and will be ignored");
						}
						| ITER ":" INTEGER ";"
						{
							if (drv.data.isIterationSet()) {
								WARNING(@4, "Iteration number already defined");
							} else {
								drv.data.setIterations($3);
							}
						}
						| PSPLITS ":" INTEGER ";"
						{
							if (drv.data.getMaxParameterSplits() > 0) {
								WARNING(@4, "The maximum number of parameter splits has been already defined");
							} else {
								drv.data.setMaxParameterSplits($3);
							}
						}
						| PRESPLITS ":" ON ";"
						{
							drv.data.setPreSplits(true);
						}
						| PRESPLITS ":" OFF ";"
						{
							drv.data.setPreSplits(false);
						}
						| MAX_MAGNITUDE ":" DOUBLE ";"
						{
							drv.data.setMaxVersorMagnitude($3);
						}

symbol			: VAR identList IN doubleInterval ";"
						{
							for (unsigned i = 0; i < $2.size(); i++)
							{
								if (drv.data.isSymbolDefined($2[i])) {
									ERROR(@2, "Symbol '" + $2[i] + "' already defined");
								} else {
									SymbolicAlgebra::Symbol<> s($2[i]);
									drv.data.addVariable(new AbsSyn::Variable(s));
									drv.ctx.addVariable(s);
									
									AbsSyn::Direction *d = new AbsSyn::Direction(s, 0, AbsSyn::Direction::Type::INT, $4.first, $4.second, "default_" + $2[i]);
									drv.data.addDirectionConstraint(d, drv.ctx);
								}
							}
						}
						| VAR identList ";"
						{
							for (unsigned i = 0; i < $2.size(); i++)
							{
								if (drv.data.isSymbolDefined($2[i])) {
									ERROR(@2, "Symbol '" + $2[i] + "' already defined");
								} else {
									SymbolicAlgebra::Symbol<> s($2[i]);
									drv.data.addVariable(new AbsSyn::Variable(s));
									drv.ctx.addVariable(s);
								}
							}
						}
						| PARAM identList IN doubleInterval ";"
						{
							for (unsigned i = 0; i < $2.size(); i++)
							{
								if (drv.data.isSymbolDefined($2[i])) {
									ERROR(@2, "Symbol '" + $2[i] + "' already defined");
								} else {
									SymbolicAlgebra::Symbol<> s($2[i]);
									drv.data.addParameter(new AbsSyn::Parameter(s));
									drv.ctx.addParameter(s);
									
									AbsSyn::Direction *d = new AbsSyn::Direction(s, 0, AbsSyn::Direction::Type::INT, $4.first, $4.second, "default_" + $2[i]);
									drv.data.addParamDirectionConstraint(d, drv.ctx);
								}
							}
						}
						| PARAM identList ";"
						{
							for (unsigned i = 0; i < $2.size(); i++)
							{
								if (drv.data.isSymbolDefined($2[i])) {
									ERROR(@2, "Symbol '" + $2[i] + "' already defined");
								} else {
									SymbolicAlgebra::Symbol<> s($2[i]);
									drv.data.addParameter(new AbsSyn::Parameter(s));
									drv.ctx.addParameter(s);
								}
							}
						}
						| CONST IDENT "=" expr ";"
						{
							if (!isNumeric($4, drv.ctx)) {
								ERROR(@3, "Expression defining constant must be numeric");
								$4 = 0;
							}
							
							if (drv.data.isSymbolDefined($2)) {
								ERROR(@2, "Symbol '" + $2 + "' already defined");
							} else {
								SymbolicAlgebra::Symbol<> s($2);
								drv.data.addConstant(new AbsSyn::Constant(s, evaluate($4, drv.ctx)));
								drv.ctx.addConstant(s, evaluate($4, drv.ctx));
							}
						}
						| DEFINE IDENT "=" expr ";"
						{
							if (drv.data.isSymbolDefined($2)) {
								ERROR(@2, "Symbol '" + $2 + "' already defined");
							} else {
								SymbolicAlgebra::Symbol<> s($2);
								drv.data.addDefinition(new AbsSyn::Definition(s, $4));
								drv.ctx.addDefinition(s, $4);
							}
						}
						| DYN "(" IDENT ")" "=" expr ";"
						{
							if (getParamDegree($6, drv.ctx) > 1) {
								ERROR(@6, "Expression in dynamic must be at most linear w.r.t. parameters");
								$6 = 0;
							}
							
							if (!drv.data.isVarDefined($3)) {
								ERROR(@3, "Symbol '" + $3 + "' is not a variable");
							} else {
								AbsSyn::Variable *v = drv.data.getVar($3);
								if (v->isDynamicDefined()) {
									WARNING(@$, "Redefinition of dynamic for variable '" + $3 + "'");
								} else {
									v->setDynamic($6);
								}
							}
						}
						| SPEC ":" formula ";"
						{
							std::shared_ptr<STL> f;
							try {
								f = $3->simplify();
							} catch (std::logic_error &e) {
								ERROR(@3, "Negations of UNTIL are not allowed");
								f = std::make_shared<Atom>(0);
							}
							
							$3.reset();
							drv.data.addSpec(f);
						}
						| ASSUME direction ";"
						{
							if (hasParams($2->getLHS(), drv.ctx) || hasParams($2->getRHS(), drv.ctx)) {
								ERROR(@2, "Expressions in assumptions cannot contain parameters");
								$2 = 0;
							}
							
							if (getVarDegree($2->getLHS(), drv.ctx) > 1 || getVarDegree($2->getRHS(), drv.ctx) > 1) {
								ERROR(@2, "Assumption must be linear");
								$2 = 0;
							}
							
							if (drv.data.getProblem() == AbsSyn::problemType::SYNTH) {
								ERROR(@1, "Assumptions are not supported for synthesis yet");
							} else if ($2->getType() == AbsSyn::Direction::Type::EQ) {
								ERROR(@2, "Directions with \"=\" are not supported yet in assumptions");
							} else if ($2->getType() == AbsSyn::Direction::Type::INT) {
								ERROR(@2, "Directions with intervals are not supported yet in assumptions");
							} else {
								drv.data.addAssumption($2);
							}
						}

matrices		: var_direction {}
						| template {}
						| param_direction {}

direction	: expr directionType expr
						{
							if (getVarDegree($1, drv.ctx) > 1) {
								ERROR(@1, "Expression in directions must be at most linear");
								$1 = 0;
							}
							if (getVarDegree($3, drv.ctx) > 1) {
								ERROR(@3, "Expression in directions must be at most linear");
								$3 = 0;
							}
							if (getParamDegree($1, drv.ctx) > 0) {
								ERROR(@1, "Expression in directions cannot contain parameters");
								$1 = 0;
							}
							if (getParamDegree($3, drv.ctx) > 0) {
								ERROR(@3, "Expression in directions cannot contain parameters");
								$3 = 0;
							}
							$$ = new AbsSyn::Direction($1, $3, $2);
						}
						| expr IN doubleInterval
						{
							if (getVarDegree($1, drv.ctx) > 1) {
								ERROR(@1, "Expression in directions must be at most linear");
							}
							if (getParamDegree($1, drv.ctx) > 0) {
								ERROR(@1, "Expression in directions cannot contain parameters");
								$$ = new AbsSyn::Direction(0, 0, AbsSyn::Direction::Type::INT, $3.first, $3.second);
							} else {
								$$ = new AbsSyn::Direction($1, 0, AbsSyn::Direction::Type::INT, $3.first, $3.second);
							}
						}

directionType		: "<"
							{
								$$ = AbsSyn::Direction::Type::LT;
							}
							| "<="
							{
								$$ = AbsSyn::Direction::Type::LE;
							}
							| ">"
							{
								$$ = AbsSyn::Direction::Type::GT;
							}
							| ">="
							{
								$$ = AbsSyn::Direction::Type::GE;
							}
							| "="
							{
								$$ = AbsSyn::Direction::Type::EQ;
							}

var_direction	: DIR direction ";"
							{
								if ($2->hasParams(drv.ctx)) {
									ERROR(@2, "Variable directions cannot contain parameters");
								}
								drv.data.addDirectionConstraint($2, drv.ctx);
							}
							| DIR IDENT ":" direction ";"
							{
								if ($4->hasParams(drv.ctx)) {
									ERROR(@2, "Variable directions cannot contain parameters");
								}
								if (drv.data.isSymbolDefined($2)) {
									ERROR(@2, "Symbol '" + $2 + "' already defined");
								}
								$4->setName($2);
								drv.data.addDirectionConstraint($4, drv.ctx);
						}

template		: TEMPL "=" "{" rowList "}"
						{
							if ($4[0].size() != drv.data.getVarNum()) {
								ERROR(@3, "template matrix must have as many columns as the number of variables");
							}
							
							if (drv.data.templateRows() != 0) {
								ERROR(@$, "Redefinition of template matrix");
							} else {
								drv.data.setTemplate($4);
							}
						}

rowList			: "{" matrixRow "}" { $$ = std::vector<std::vector<int>>{$2}; }
						| rowList "," "{" matrixRow "}" { $1.push_back($4); $$ = $1; }

matrixRow	: matrixRow "," IDENT
					{
						if (!drv.data.isDirectionDefined($3)) {
							ERROR(@3, "Direction " + $3 + " is not defined");
							$1.push_back(0);
							$$ = $1;
						} else {
							$1.push_back(drv.data.findDirectionPos($3));
							$$ = $1;
						}
					}
					| IDENT
					{
						if (!drv.data.isDirectionDefined($1)) {
							ERROR(@1, "Direction " + $1 + " is not defined");
							$$ = {0};
						} else {
							$$ = {drv.data.findDirectionPos($1)};
						}
					}
					| matrixRow "," INTEGER
					{
						if ($3 < 0 || (unsigned int)$3 >= drv.data.getDirectionsNum()) {
							ERROR(@3, "Unknown direction " + std::to_string($3));
						}
						
						$1.push_back($3);
						$$ = $1;
					}
					| INTEGER
					{
						if ($1 < 0 || (unsigned int)$1 >= drv.data.getDirectionsNum()) {
							ERROR(@1, "Unknown direction " + std::to_string($1));
						}
						
						$$ = {$1};
					}

param_direction	: PDIR direction ";"
								{
									if ($2->hasVars(drv.ctx)) {
										ERROR(@2, "Parameter directions cannot contain variables");
									}
									
									drv.data.addParamDirectionConstraint($2, drv.ctx);
								}
								| PDIR IDENT ":" direction ";"
								{
									if ($4->hasVars(drv.ctx)) {
										ERROR(@2, "Parameter directions cannot contain variables");
										
									}
									
									if (drv.data.isSymbolDefined($2)) {
										ERROR(@2, "Symbol '" + $2 + "' already defined");
									}
									$4->setName($2);
									drv.data.addParamDirectionConstraint($4, drv.ctx);
								}

problemType	: REACH { $$ = AbsSyn::problemType::REACH; }
						| SYNTH { $$ = AbsSyn::problemType::SYNTH; }

modeType	: BOX { $$ = AbsSyn::modeType::BOX; }
					| PARAL { $$ = AbsSyn::modeType::PARAL; }
					| POLY { $$ = AbsSyn::modeType::POLY; }

identList	: IDENT { $$ = std::vector<std::string>{$1}; }
					| identList "," IDENT { $1.push_back($3); $$ = $1; }

intInterval			: "[" expr "," expr "]"
								{
									double x1, x2;
									
									if (!isNumeric($2, drv.ctx)) {
										ERROR(@2, "Intervals require numeric expressions");
										x1 = 0;
									} else {
										x1 = evaluate($2, drv.ctx);
									}
									
									if (!isNumeric($4, drv.ctx)) {
										ERROR(@4, "Intervals require numeric expressions");
										x2 = 1;
									} else {
										x2 = evaluate($4, drv.ctx);
									}
									
									if ((int) x1 != x1) {
										WARNING(@2, "Left endopoint in interval is not integer and will be truncated");
									}
									
									if ((int) x2 != x2) {
										WARNING(@4, "Left endopoint in interval is not integer and will be truncated");
									}
									
									if ((int) x2 < (int) x1) {
										ERROR(@$, "Right endpoint must be greater than or equal to the left one");
										x1 = 0;
										x2 = 1;
									}
									
									$$ = std::pair<int, int>((int) x1, (int) x2);
								}

doubleInterval	: "[" expr "," expr "]"
								{
									double x1, x2;
									
									if (!isNumeric($2, drv.ctx)) {
										ERROR(@2, "Intervals require numeric expressions");
										x1 = 0;
									} else {
										x1 = evaluate($2, drv.ctx);
									}
									
									if (!isNumeric($4, drv.ctx)) {
										ERROR(@4, "Intervals require numeric expressions");
										x2 = 1;
									} else {
										x2 = evaluate($4, drv.ctx);
									}
									
									if (x2 < x1) {
										ERROR(@$, "Right endpoint must be greater than or equal to the left one");
										x1 = 0;
										x2 = 1;
									}
									
									$$ = std::pair<double, double>(x1, x2);
								}

number		: DOUBLE { $$ = $1; }
					| INTEGER { $$ = (double) $1; }

expr		: number	{ $$ = $1; }
				| "+" expr { $$ = $2; }
				| IDENT	
				{
					if (!drv.data.isSymbolDefined($1)) {
						ERROR(@1, "Symbol " + $1 + " is undefined");
						$$ = 0;
					} else {
						$$ = drv.ctx.getSymbol($1);
					}
				}
				| expr "*" expr { $$ = $1 * $3; }
				| expr "^" expr
				{
					double val;
					
					if (!isNumeric($3, drv.ctx)) {
						ERROR(@3, "Exponent must be numeric");
						val = 1;
					} else {
						val = $3.evaluate<double>();
					}
					
					if (val != (int) val) {
						ERROR(@3, "Exponent must be integer");
					}
					
					int exp = (int) val;
					if (exp < 0) {
						ERROR(@3, "Exponent must be non-negative");
					}
					
					if (exp == 0) {
						$$ = 1;
					} else {
						SymbolicAlgebra::Expression<> res = $1;
						for (int i = 1; i < exp; i++) {
							res *= $1;
						}
						$$ = res;
					}
				}
				| expr "/" expr
				{
					if (!isNumeric($3, drv.ctx)) {
						ERROR(@3, "cannot divide by non--numeric expression");
						$$ = $1 / 1;
					} else {
						$$ = $1 / $3;
					}
				}
				| expr "+" expr { $$ = $1 + $3; }
				| expr "-" expr { $$ = $1 - $3; }
				| "-" expr %prec UMINUS { $$ = -$2; }
				| "(" expr ")" { $$ = $2;}

formula	: expr ">" expr { $$ = std::make_shared<Atom>($3 - $1); }
				| expr ">=" expr { $$ = std::make_shared<Atom>($3 - $1); }
				| expr "<" expr { $$ = std::make_shared<Atom>($1 - $3); }
				| expr "<=" expr { $$ = std::make_shared<Atom>($1 - $3); }
				| expr "=" expr
				{
					std::shared_ptr<STL> f1 = std::make_shared<Atom>($1 - $3);
					std::shared_ptr<STL> f2 = std::make_shared<Atom>($3 - $1);
					$$ = std::make_shared<Conjunction>(f1, f2);
				}
				| formula AND formula		{ $$ = std::make_shared<Conjunction>($1, $3); }
				| formula OR formula		{ $$ = std::make_shared<Disjunction>($1, $3); }
				| NOT formula									{ $$ = std::make_shared<Negation>($2); }
				| "(" formula ")" { $$ = $2; }
				| "G" intInterval formula %prec "G"	{ $$ = std::make_shared<Always>($2.first, $2.second, $3); }
				| "F" intInterval formula %prec "F"	{ $$ = std::make_shared<Eventually>($2.first, $2.second, $3); }
				| formula "U" intInterval formula %prec "U"	{ $$ = std::make_shared<Until>($1, $3.first, $3.second, $4); }

option	: OPT TRANS transType ";"
				{
					if (drv.data.isTransModeDefined()) {
						WARNING(@$, "Transformation type already defined");
					} else {
						drv.data.setTransMode($3);
					}
				}
				| OPT DECOMP
				{
					if (drv.data.isDecompositionDefined()) {
						WARNING(@$, "Decomposition option already defined");
					} else {
						drv.data.setDecomposition();
					}
				}
				| OPT ALPHA DOUBLE ";"
				{
					if (drv.data.isAlphaDefined()) {
						WARNING(@$, "Alpha already defined");
					}
					
					if ($3 > 1) {
						ERROR(@3, "Alpha must be between 0 and 1");
					}
					
					drv.data.setAlpha($3);
				}

transType : AFO { $$ = AbsSyn::transType::AFO; }
					| OFO { $$ = AbsSyn::transType::OFO; }

%%

void
yy::parser::error (const location_type& l, const std::string& m)
{
	std::cerr << "Error at line " << l.end.line << ", column " << l.end.column << ": " << m << '\n';
}
