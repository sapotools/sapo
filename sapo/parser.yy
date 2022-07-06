%skeleton "lalr1.cc" // -*- C++ -*-
//%require "3.8.1"
%require "3.5.1"
%language "c++"
%defines

%define api.token.constructor
%define api.value.type variant
%define parse.assert

%code requires {
	#include <string>
	#include <fstream>
	#include <cmath>
	
	#include <LinearAlgebraIO.h>

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
	
	#define YYDEBUG 1
	
	class driver;
	
	// macro for errors
	#define ERROR(loc, msg)\
		drv.error(loc, msg, drv.file);
	
	// macro for warnings
	#define WARNING(loc, msg)\
		drv.warning(loc, msg, drv.file);
	
	#define MISSING_SC(loc)\
		drv.missingSemicolon(loc, drv.file);\
		drv.errors = true;
	
	std::string possibleStatements(std::string s);
}

// The parsing context.
%param { driver& drv }

%locations
%define api.location.file "../include/locations.h"

%define parse.trace
//%define parse.error custom
%define parse.error detailed

%code {
#include "driver.h"
}

%define api.token.prefix {TOK_}
%token
	END 0
	PROB 
	REACH
	INVCHECK
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
	ON
	OFF
	MAX_MAGNITUDE
	DIR
	TEMPL
	PDIR
	OPT
	K_IND_JOIN
	LISTING
	PACKAGING
	MERGING
	NO_CACHE
	TRANS
	AFO
	OFO
	DECOMP
	ALPHA
	COMPOSE
	INVARIANT
	EXP_DELTA
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
%token <unsigned int> NATURAL
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
%nterm <std::vector<unsigned int>> matrixRow
%nterm <std::vector<std::vector<unsigned int>>> rowList
%nterm <std::shared_ptr<STL::STL> > formula
%nterm <AbsSyn::transType> transType
%nterm <AbsSyn::Direction *> direction
%nterm <AbsSyn::Direction::Type> directionType
%nterm <bool> presplits_value

%printer { yyo << $$; } <*>;

%start s

%%

s		: statement {}
		| s statement {}
		| END
		{
			ERROR(@1, "Empty file");
			YYERROR;
		}

statement		: header {}
						| symbol {}
						| matrices {}
						| footer {}
						| ";" { /* empty statement */ }
						| IDENT error ";"
						{
							ERROR(@1, "\"" + $1 + "\" is not a valid statement." + possibleStatements($1));
						}
						| error ";"
						{
							ERROR(@1, "Syntax error in statement");
						}

header			: PROB ":" problemType ";"
						{
							if (drv.data.isProblemDefined()) {
								ERROR(@$, "Problem has already been defined");
							} else {
								drv.data.setProblem($3);
							}
						}
						| PROB ":" error ";"
						{
							ERROR(@4, "Problem type must be \"reachability\" or \"synthesis\"");
							YYERROR;
						}
						| PROB ":" problemType error
						{
							MISSING_SC(@3);
						}
						| PARAMMODE ":" modeType ";"
						{
							WARNING(@$, "Parameter modality is deprecated and will be ignored");
						}
						| PARAMMODE ":" modeType error
						{
							MISSING_SC(@3);
						}
						| VARMODE ":" modeType ";"
						{
							WARNING(@$, "Variable modality is deprecated and will be ignored");
						}
						| VARMODE ":" modeType error
						{
							MISSING_SC(@3);
						}
						| ITER ":" NATURAL ";"
						{
							if (drv.data.isIterationSet()) {
								WARNING(@$, "Iteration number already defined");
							} else {
								drv.data.setIterations($3);
							}
						}
						| ITER ":" NATURAL error
						{
							MISSING_SC(@3);
						}
						| PSPLITS ":" NATURAL ";"
						{
							if (drv.data.getMaxParameterSplits() > 0) {
								WARNING(@$, "The maximum number of parameter splits has been already defined");
							} else {
								drv.data.setMaxParameterSplits($3);
							}
						}
						| PSPLITS ":" NATURAL error
						{
							MISSING_SC(@3);
						}
						| PRESPLITS ":" presplits_value ";"
						{
							drv.data.setPreSplits($3);
						}
						| PRESPLITS ":" presplits_value error
						{
							MISSING_SC(@3);
						}
						| PRESPLITS ":" error ";"
						{
							ERROR(@3, "Value for presplits is either \"on\" or \"off\"");
							yyerrok;
						}
						| MAX_MAGNITUDE ":" DOUBLE ";"
						{
							drv.data.setMaxVersorMagnitude($3);
						}
						| MAX_MAGNITUDE ":" DOUBLE error
						{
							MISSING_SC(@3);
						}

presplits_value	: ON { $$ = true; }
								| OFF { $$ = false; }

symbol			: VAR identList IN doubleInterval ";"
						{
							for (unsigned i = 0; i < $2.size(); i++)
							{
								if (drv.data.isSymbolDefined($2[i])) {
									ERROR(@2, "Symbol '" + $2[i] + "' already defined");
								} else {
									SymbolicAlgebra::Symbol<> s($2[i]);
									drv.data.addVariable(new AbsSyn::Variable(s));
									
									std::string str = "default_" + $2[i];
									SymbolicAlgebra::Symbol<> *sym = new SymbolicAlgebra::Symbol<>(str);
									AbsSyn::Direction *d = new AbsSyn::Direction(s, $4.first, $4.second, sym);
									drv.data.addVarDirectionConstraint(d);
								}
							}
						}
						| VAR identList IN doubleInterval error
						{
							MISSING_SC(@4);
							for (unsigned i = 0; i < $2.size(); i++)
							{
								if (drv.data.isSymbolDefined($2[i])) {
									ERROR(@2, "Symbol '" + $2[i] + "' already defined");
								} else {
									SymbolicAlgebra::Symbol<> s($2[i]);
									drv.data.addVariable(new AbsSyn::Variable(s));
									
									std::string str = "default_" + $2[i];
									SymbolicAlgebra::Symbol<> *sym = new SymbolicAlgebra::Symbol<>(str);
									AbsSyn::Direction *d = new AbsSyn::Direction(s, $4.first, $4.second, sym);
									drv.data.addVarDirectionConstraint(d);
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
								}
							}
						}
						| VAR identList error
						{
							MISSING_SC(@2);
							for (unsigned i = 0; i < $2.size(); i++)
							{
								if (drv.data.isSymbolDefined($2[i])) {
									ERROR(@2, "Symbol '" + $2[i] + "' already defined");
								} else {
									SymbolicAlgebra::Symbol<> s($2[i]);
									drv.data.addVariable(new AbsSyn::Variable(s));
								}
							}
						}
						| VAR error ";"
						{
							ERROR(@2, "Error in variable declaration");
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
									
									std::string str = "default_" + $2[i];
									SymbolicAlgebra::Symbol<> *sym = new SymbolicAlgebra::Symbol<>(str);
									AbsSyn::Direction *d = new AbsSyn::Direction(s, $4.first, $4.second, sym);
									drv.data.addParamDirectionConstraint(d);
								}
							}
						}
						| PARAM identList IN doubleInterval error
						{
							MISSING_SC(@4);
							for (unsigned i = 0; i < $2.size(); i++)
							{
								if (drv.data.isSymbolDefined($2[i])) {
									ERROR(@2, "Symbol '" + $2[i] + "' already defined");
								} else {
									SymbolicAlgebra::Symbol<> s($2[i]);
									drv.data.addParameter(new AbsSyn::Parameter(s));
									
									std::string str = "default_" + $2[i];
									SymbolicAlgebra::Symbol<> *sym = new SymbolicAlgebra::Symbol<>(str);
									AbsSyn::Direction *d = new AbsSyn::Direction(s, $4.first, $4.second, sym);
									drv.data.addParamDirectionConstraint(d);
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
								}
							}
						}
						| PARAM identList error
						{
							MISSING_SC(@2);
							for (unsigned i = 0; i < $2.size(); i++)
							{
								if (drv.data.isSymbolDefined($2[i])) {
									ERROR(@2, "Symbol '" + $2[i] + "' already defined");
								} else {
									SymbolicAlgebra::Symbol<> s($2[i]);
									drv.data.addParameter(new AbsSyn::Parameter(s));
								}
							}
						}
						| CONST IDENT "=" expr ";"
						{
							if (!AbsSyn::isNumeric($4)) {
								ERROR(@3, "Expression defining constant must be numeric");
								$4 = 0;
							}
							
							if (drv.data.isSymbolDefined($2)) {
								ERROR(@2, "Symbol '" + $2 + "' already defined");
							} else {
								SymbolicAlgebra::Symbol<> s($2);
								drv.data.addConstant(new AbsSyn::Constant(s, $4.evaluate()));
							}
						}
						| CONST IDENT "=" expr error
						{
							MISSING_SC(@4);
							if (!AbsSyn::isNumeric($4)) {
								ERROR(@3, "Expression defining constant must be numeric");
								$4 = 0;
							}
							
							if (drv.data.isSymbolDefined($2)) {
								ERROR(@2, "Symbol '" + $2 + "' already defined");
							} else {
								SymbolicAlgebra::Symbol<> s($2);
								drv.data.addConstant(new AbsSyn::Constant(s, $4.evaluate()));
							}
						}
						| DEFINE IDENT "=" expr ";"
						{
							if (drv.data.isSymbolDefined($2)) {
								ERROR(@2, "Symbol '" + $2 + "' already defined");
							} else {
								SymbolicAlgebra::Symbol<> s($2);
								drv.data.addDefinition(new AbsSyn::Definition(s, $4));
							}
						}
						| DEFINE IDENT "=" expr error
						{
							MISSING_SC(@4);
							if (drv.data.isSymbolDefined($2)) {
								ERROR(@2, "Symbol '" + $2 + "' already defined");
							} else {
								SymbolicAlgebra::Symbol<> s($2);
								drv.data.addDefinition(new AbsSyn::Definition(s, $4));
							}
						}
						| DYN "(" IDENT ")" "=" expr ";"
						{
							if (AbsSyn::getDegree($6, drv.data.getParamSymbols()) > 1) {
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
						| DYN "(" IDENT ")" "=" expr error
						{
							MISSING_SC(@6);
						}
						| DYN "(" IDENT error "=" expr ";"
						{
							ERROR(@3, "Missing \")\"");
						}
						| SPEC ":" formula ";"
						{
							std::shared_ptr<STL::STL> f;
							try {
								f = $3->get_PNF();
							} catch (std::logic_error &e) {
								ERROR(@3, "Negations of UNTIL are not allowed");
								f = std::make_shared<STL::Atom>(0);
							}
							
							$3.reset();
							drv.data.addSpec(f);
						}
						| SPEC ":" formula error
						{
							MISSING_SC(@3);
						}
						| INVARIANT ":" invariant ";"
						{}
						| INVARIANT ":" invariant error
						{
							MISSING_SC(@3);
						}
						| ASSUME direction ";"
						{
							if ($2->contains(drv.data.getParamSymbols())) {
								ERROR(@2, "Expressions in assumptions cannot contain parameters");
								$2 = 0;
							}
							
							if (AbsSyn::getDegree($2->getLHS() - $2->getRHS(), drv.data.getVarSymbols()) > 1) {
								ERROR(@2, "Assumption must be linear");
								$2 = 0;
							}
							
							if (drv.data.getProblem() == AbsSyn::problemType::SYNTH) {
								ERROR(@1, "Assumptions are not supported for synthesis yet");
							} else if ($2->getType() == AbsSyn::Direction::Type::EQ) {
								ERROR(@2, "Directions with \"=\" are not supported yet in assumptions");
							} else if ($2->getType() == AbsSyn::Direction::Type::IN) {
								ERROR(@2, "Directions with intervals are not supported yet in assumptions");
							} else {
								drv.data.addAssumption($2);
							}
						}
						| ASSUME direction error
						{
							MISSING_SC(@2);
						}

matrices		: var_direction {}
				| template {}
				| param_direction {}

direction	: expr directionType expr
						{
							if (AbsSyn::getDegree($1, drv.data.getVarSymbols()) > 1) {
								ERROR(@1, "Expression in directions must be at most linear");
								$1 = 0;
							}
							if (AbsSyn::getDegree($3, drv.data.getVarSymbols()) > 1) {
								ERROR(@3, "Expression in directions must be at most linear");
								$3 = 0;
							}
							if (AbsSyn::getDegree($1, drv.data.getParamSymbols()) > 0) {
								ERROR(@1, "Expression in directions cannot contain parameters");
								$1 = 0;
							}
							if (AbsSyn::getDegree($3, drv.data.getParamSymbols()) > 0) {
								ERROR(@3, "Expression in directions cannot contain parameters");
								$3 = 0;
							}
							$$ = new AbsSyn::Direction($1, $3, $2);
						}
						| expr IN doubleInterval
						{
							if (AbsSyn::getDegree($1, drv.data.getVarSymbols()) > 1) {
								ERROR(@1, "Expression in directions must be at most linear");
							}
							if (AbsSyn::getDegree($1, drv.data.getParamSymbols()) > 0) {
								ERROR(@1, "Expression in directions cannot contain parameters");
								$$ = new AbsSyn::Direction(0, $3.first, $3.second);
							} else {
								$$ = new AbsSyn::Direction($1, $3.first, $3.second);
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

invariant 	: 	direction
				{
					if ($1->contains(drv.data.getParamSymbols())) {
						ERROR(@1, "Invariant formula cannot contain parameters");
					}
					drv.data.addInvariantConstraint($1);
				}
			|	invariant AND direction
				{
					if ($3->contains(drv.data.getParamSymbols())) {
						ERROR(@3, "Invariant formula cannot contain parameters");
					}
					drv.data.addInvariantConstraint($3);
				}

var_direction	: DIR direction ";"
							{
								if ($2->contains(drv.data.getParamSymbols())) {
									ERROR(@2, "Variable directions cannot contain parameters");
								}
								drv.data.addVarDirectionConstraint($2);
							}
							| DIR direction error
							{
								MISSING_SC(@2);
								if ($2->contains(drv.data.getParamSymbols())) {
									ERROR(@2, "Variable directions cannot contain parameters");
								}
								drv.data.addVarDirectionConstraint($2);
							}
							| DIR error ";"
							{
								ERROR(@2, "Error in direction");
							}
							| DIR IDENT ":" direction ";"
							{
								if ($4->contains(drv.data.getParamSymbols())) {
									ERROR(@2, "Variable directions cannot contain parameters");
								}
								if (drv.data.isSymbolDefined($2)) {
									ERROR(@2, "Symbol '" + $2 + "' already defined");
								}
								SymbolicAlgebra::Symbol<> *s = new SymbolicAlgebra::Symbol<>($2);
								$4->setSymbol(s);
								drv.data.addVarDirectionConstraint($4);
							}
							| DIR IDENT ":" direction error
							{
								MISSING_SC(@3);
								if ($4->contains(drv.data.getParamSymbols())) {
									ERROR(@2, "Variable directions cannot contain parameters");
								}
								if (drv.data.isSymbolDefined($2)) {
									ERROR(@2, "Symbol '" + $2 + "' already defined");
								}
								SymbolicAlgebra::Symbol<> *s = new SymbolicAlgebra::Symbol<>($2);
								$4->setSymbol(s);
								drv.data.addVarDirectionConstraint($4);
							}
							| DIR IDENT error direction ";"
							{
								ERROR(@2, "Missing \":\"");
							}
							| DIR IDENT error ";"
							{
								ERROR(@3, "Error in direction");
							}

template		: TEMPL "=" "{" rowList "}" ";"
						{
							if ($4[0].size() != drv.data.getVarNum()) {
								ERROR(@4, "template matrix must have as many columns as the number of variables");
							}
							
							if (drv.data.templateRows() != 0) {
								ERROR(@$, "Redefinition of template matrix");
							} else {
								drv.data.setTemplate($4);
							}
						}
						| TEMPL "=" "{" rowList "}" error
						{
							MISSING_SC(@5);
							if ($4[0].size() != drv.data.getVarNum()) {
								ERROR(@4, "template matrix must have as many columns as the number of variables");
							}
							
							if (drv.data.templateRows() != 0) {
								ERROR(@$, "Redefinition of template matrix");
							} else {
								drv.data.setTemplate($4);
							}
						}
						| TEMPL "=" "{" rowList error ";"
						{
							ERROR(@4, "Missing \"}\"");
							if ($4[0].size() != drv.data.getVarNum()) {
								ERROR(@4, "template matrix must have as many columns as the number of variables");
							}
							
							if (drv.data.templateRows() != 0) {
								ERROR(@$, "Redefinition of template matrix");
							} else {
								drv.data.setTemplate($4);
							}
						}
						| TEMPL "=" "{" "}" ";"
						{
							ERROR(@$, "Template cannot be empty");
						}
						| TEMPL error ";"
						{
							ERROR(@2, "Error in definition of template");
						}

rowList			: "{" matrixRow "}" { $$ = std::vector<std::vector<unsigned int>>{$2}; }
						| "{" matrixRow error
						{
							ERROR(@2, "Missing \"}\""); $$ = std::vector<std::vector<unsigned int>>{$2};
						}
						| rowList "," "{" matrixRow "}" { $1.push_back($4); $$ = $1; }
						| rowList "," "{" matrixRow error
						{
							ERROR(@4, "Missing \"}\"");
							$1.push_back($4);
							$$ = $1;
						}
						| rowList "," error
						{
							ERROR(@2, "Unexpected \",\"");
							$$ = $1;
						}

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
					| matrixRow "," NATURAL
					{
						if ($3 >= drv.data.getDirectionsNum()) {
							ERROR(@3, "Unknown direction " + std::to_string($3));
						}
						
						$1.push_back($3);
						$$ = $1;
					}
					| NATURAL
					{
						if ($1 >= drv.data.getDirectionsNum()) {
							ERROR(@1, "Unknown direction " + std::to_string($1));
						}
						
						$$ = {$1};
					}
					| matrixRow "," error
					{
						ERROR(@2, "Unexpected \",\"");
						$$ = $1;
					}
					| %empty
					{
						$$ = {};
					}

param_direction	: PDIR direction ";"
								{
									if ($2->contains(drv.data.getVarSymbols())) {
										ERROR(@2, "Parameter directions cannot contain variables");
									}
									
									drv.data.addParamDirectionConstraint($2);
								}
								| PDIR direction error
								{
									MISSING_SC(@2);
									if ($2->contains(drv.data.getVarSymbols())) {
										ERROR(@2, "Parameter directions cannot contain variables");
									}
									
									drv.data.addParamDirectionConstraint($2);
								}
								| PDIR error ";"
								{
									ERROR(@2, "Error in parameter direction");
								}
								| PDIR IDENT ":" direction ";"
								{
									if ($4->contains(drv.data.getVarSymbols())) {
										ERROR(@2, "Parameter directions cannot contain variables");
										
									}
									
									if (drv.data.isSymbolDefined($2)) {
										ERROR(@2, "Symbol '" + $2 + "' already defined");
									}
									SymbolicAlgebra::Symbol<> *s = new SymbolicAlgebra::Symbol<>($2);
									$4->setSymbol(s);
									drv.data.addParamDirectionConstraint($4);
								}
								| PDIR IDENT ":" direction error
								{
									MISSING_SC(@4);
									if ($4->contains(drv.data.getVarSymbols())) {
										ERROR(@2, "Parameter directions cannot contain variables");
										
									}
									
									if (drv.data.isSymbolDefined($2)) {
										ERROR(@2, "Symbol '" + $2 + "' already defined");
									}
									SymbolicAlgebra::Symbol<> *s = new SymbolicAlgebra::Symbol<>($2);
									$4->setSymbol(s);
									drv.data.addParamDirectionConstraint($4);
								}
								| PDIR IDENT error direction ";"
								{
									ERROR(@2, "Missing \":\"");
								}
								| PDIR IDENT error ";"
								{
									ERROR(@3, "Error in parameter direction");
								}

problemType	: REACH { $$ = AbsSyn::problemType::REACH; }
						| SYNTH { $$ = AbsSyn::problemType::SYNTH; }
						| INVCHECK { $$ = AbsSyn::problemType::INVARIANT; }

modeType	: BOX { $$ = AbsSyn::modeType::BOX; }
					| PARAL { $$ = AbsSyn::modeType::PARAL; }
					| POLY { $$ = AbsSyn::modeType::POLY; }

identList	: IDENT { $$ = std::vector<std::string>{$1}; }
					| identList "," IDENT { $1.push_back($3); $$ = $1; }

intInterval			: "[" expr "," expr "]"
								{
									double x1, x2;
									
									if (!AbsSyn::isNumeric($2)) {
										ERROR(@2, "Intervals require numeric expressions");
										x1 = 0;
									} else {
										x1 = $2.evaluate();
									}
									
									if (!AbsSyn::isNumeric($4)) {
										ERROR(@4, "Intervals require numeric expressions");
										x2 = 1;
									} else {
										x2 = $4.evaluate();
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
								| "[" expr ";" expr error
								{
									ERROR(@4, "Missing \"]\"");
									$$ = std::pair<int, int>(0, 1);
								}

doubleInterval	: "[" expr "," expr "]"
								{
									double x1, x2;
									
									if (!AbsSyn::isNumeric($2)) {
										ERROR(@2, "Intervals require numeric expressions");
										x1 = 0;
									} else {
										x1 = $2.evaluate();
									}
									
									if (!AbsSyn::isNumeric($4)) {
										ERROR(@4, "Intervals require numeric expressions");
										x2 = 1;
									} else {
										x2 = $4.evaluate();
									}
									
									if (x2 < x1) {
										ERROR(@$, "Right endpoint must be greater than or equal to the left one");
										x1 = 0;
										x2 = 1;
									}
									
									$$ = std::pair<double, double>(x1, x2);
								}
								| "[" expr "," expr error
								{
									ERROR(@4, "Missing \"]\"");
									$$ = std::pair<double, double>(0, 1);
								}

number		: DOUBLE { $$ = $1; }
				| NATURAL { $$ = (double) $1; }

expr		: number	{ $$ = $1; }
				| "+" expr { $$ = $2; }
				| IDENT	
				{
					if (!drv.data.isSymbolDefined($1)) {
						ERROR(@1, "Symbol " + $1 + " is undefined");
						$$ = 0;
					} else if (drv.data.isConstDefined($1)) {
						$$ = drv.data.getConst($1)->getValue();
					} else if (drv.data.isDefDefined($1)) {
						$$ = drv.data.getDef($1)->getValue();
					} else {
						$$ = drv.data.getSymbol($1);
					}
				}
				| expr "*" expr { $$ = $1 * $3; }
				| expr "^" expr
				{
					double val;
					
					if (!AbsSyn::isNumeric($3)) {
						ERROR(@3, "Exponent must be numeric");
						val = 1;
					} else {
						val = $3.evaluate();
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
					if (!AbsSyn::isNumeric($3)) {
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
				| "(" expr error
				{
					ERROR(@2, "Missing \")\"");
					$$ = $2;
				}

formula	: expr ">" expr { $$ = std::make_shared<STL::Atom>($3 - $1); }
		| expr ">=" expr { $$ = std::make_shared<STL::Atom>($3 - $1); }
		| expr "<" expr { $$ = std::make_shared<STL::Atom>($1 - $3); }
		| expr "<=" expr { $$ = std::make_shared<STL::Atom>($1 - $3); }
		| expr "=" expr
		{
			std::shared_ptr<STL::STL> f1 = std::make_shared<STL::Atom>($1 - $3);
			std::shared_ptr<STL::STL> f2 = std::make_shared<STL::Atom>($3 - $1);
			$$ = std::make_shared<STL::Conjunction>(f1, f2);
		}
		| formula AND formula		{ $$ = std::make_shared<STL::Conjunction>($1, $3); }
		| formula OR formula		{ $$ = std::make_shared<STL::Disjunction>($1, $3); }
		| NOT formula									{ $$ = std::make_shared<STL::Negation>($2); }
		| "(" formula ")" { $$ = $2; }
		| "(" formula error
		{
			ERROR(@2, "Missing \"(\"");
			$$ = $2;
		}
		| "G" intInterval formula %prec "G"	{ $$ = std::make_shared<STL::Always>($2.first, $2.second, $3); }
		| "F" intInterval formula %prec "F"	{ $$ = std::make_shared<STL::Eventually>($2.first, $2.second, $3); }
		| formula "U" intInterval formula %prec "U"	{ $$ = std::make_shared<STL::Until>($1, $3.first, $3.second, $4); }

footer	: OPT option {}
		| OPT error ";"
		{
			ERROR(@2, "Syntax error in option");
		}

option	: TRANS transType ";"
		{
			if (drv.data.isTransModeDefined()) {
				WARNING(@$, "Transformation type already defined");
			} else {
				drv.data.setTransMode($2);
			}
		}
		| TRANS transType error
		{
			MISSING_SC(@2);
			if (drv.data.isTransModeDefined()) {
				WARNING(@$, "Transformation type already defined");
			} else {
				drv.data.setTransMode($2);
			}
		}
		| TRANS error ";"
		{
			ERROR(@2, "Unknown transformation type");
		}
		| DECOMP ";"
		{
			if (drv.data.isDecompositionDefined()) {
				WARNING(@$, "Decomposition option already defined");
			} else {
				drv.data.setDecomposition();
			}
		}
		| DECOMP error
		{
			MISSING_SC(@1);
			if (drv.data.isDecompositionDefined()) {
				WARNING(@$, "Decomposition option already defined");
			} else {
				drv.data.setDecomposition();
			}
		}
		| ALPHA DOUBLE ";"
		{
			if (drv.data.isAlphaDefined()) {
				WARNING(@$, "Alpha already defined");
			}
			
			if ($2 > 1) {
				ERROR(@2, "Alpha must be between 0 and 1");
			}
			
			drv.data.setAlpha($2);
		}
		| ALPHA DOUBLE error
		{
			MISSING_SC(@2);
			if (drv.data.isAlphaDefined()) {
				WARNING(@$, "Alpha already defined");
			}
			
			if ($2 > 1) {
				ERROR(@2, "Alpha must be between 0 and 1");
			}
			
			drv.data.setAlpha($2);
		}
		| COMPOSE NATURAL ";"
		{
			if ($2 < 1) {
				yy::parser::error(@2, "Degree of composing dynamic must be at least 1");
			}
			drv.data.setDynamicDegree($2);
		}
		| COMPOSE NATURAL error
		{
			MISSING_SC(@2);
			if ($2 < 1) {
				yy::parser::error(@2, "Degree of composing dynamic must be at least 1");
			}
		}
		| K_IND_JOIN LISTING ";"
		{
			drv.data.setApproxType(Sapo::NO_APPROX);
		}
		| K_IND_JOIN PACKAGING ";"
		{
			drv.data.setApproxType(Sapo::FULL_JOIN);
		}
		| K_IND_JOIN MERGING ";"
		{
			drv.data.setApproxType(Sapo::CHAIN_JOIN);
		}
		| NO_CACHE ";"
		{
			drv.data.setBernsteinCaching(false);
		}

transType : AFO { $$ = AbsSyn::transType::AFO; }
					| OFO { $$ = AbsSyn::transType::OFO; }

%%

void yy::parser::error (const location_type& l, const std::string& m)
{
	drv.error(l, m, drv.file);
}

/*void yy::parser::report_syntax_error (const yy::parser::context &ctx) const
{
	(void) ctx;
	
	return;
	
#define N 30
	
	std::cerr << "Syntax error at " << ctx.location() << std::endl;
	std::cerr << "Unexpected " << ctx.lookahead().name() << ", ";
	yy::parser::symbol_kind_type toks[N] = {};
	ctx.expected_tokens(toks, N);
	std::cerr << "expected tokens: ";
	for (unsigned i = 0; i < N; i++) {
		std::cerr << toks[i] << ", ";
	}
	std::cerr << std::endl;
}*/


inline unsigned editDistance(std::string s1, std::string s2)
{
	using namespace std;
	
	vector<vector<unsigned>> dists(s1.size(), vector<unsigned>(s2.size()));
	
	for (unsigned i = 0; i < s1.size(); i++) {
		dists[i][0] = i;
	}
	for (unsigned i = 0; i < s2.size(); i++) {
		dists[0][i] = i;
	}
	
	for (unsigned i = 1; i < s1.size(); i++) {
		for (unsigned j = 1; j < s2.size(); j++) {
			if (s1[i] == s2[j]) {
					dists[i][j] = dists[i-1][j-1];
			} else {
				dists[i][j] = 1 + min({
					dists[i-1][j],
					dists[i][j-1],
					dists[i-1][j-1]
				});
			}
		}
	}
	
	return dists[s1.size() - 1][s2.size() - 1];
}

std::string possibleStatements(std::string s)
{
	using namespace std;
	
	vector<string> statements = {
		"problem",
		"variable_mode",
		"parameter_mode",
		"var",
		"param",
		"const",
		"define",
		"dynamic",
		"invariant",
		"spec",
		"assume",
		"iterations",
		"max_parameter_splits",
		"presplit_parameters",
		"max_bundle_magnitude",
		"direction",
		"parameter_direction",
		"template",
		"option"
	};
	
	vector<string> good {};
	for (unsigned i = 0; i < statements.size(); i++) {
		double l = max(s.size(), statements[i].size());
		if (editDistance(s, statements[i]) / l <= 0.34) {
			good.push_back(statements[i]);
		}
	}
	
	if (good.size() == 0) {
		return "";
	}
	
	string res = " Did you mean ";
	for (unsigned i = 0; i < good.size(); i++) {
		res += good[i] + (i == good.size() - 1 ? "" : (i == good.size() - 2 ? " or " : ", "));
	}
	res += "?";
	
	return res;
}
