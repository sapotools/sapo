%skeleton "lalr1.cc" // -*- C++ -*-
%require "3.5.1"
%defines

%define api.token.constructor
%define api.value.type variant
%define parse.assert

%code requires {
	# include <string>
	#include "AbsSyn.h"
	class driver;
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

%nterm <AbsSyn::problemType> problemType
%nterm <AbsSyn::modeType> modeType
%nterm <std::pair<double, double>> doubleInterval
%nterm <std::pair<int, int>> intInterval
%nterm <double> number
%nterm <std::vector<std::string>> identList
%nterm <AbsSyn::Expr *> expr
%nterm <std::vector<double>> numList
%nterm <std::vector<int>> matrixRow intList
%nterm <std::vector<std::vector<int>>> _matrix rowList
%nterm <AbsSyn::Formula *> path_formula state_formula
%nterm <AbsSyn::transType> transType

%printer { yyo << $$; } <*>;

%start s

%%

s		: headerList
				{
					if (!drv.data.isProblemDefined())
					{
						yy::parser::error(@1, "Problem type must be defined");
						YYERROR;
					}
					
					if (!drv.data.isVarModeDefined())
						drv.data.setVarMode(AbsSyn::modeType::BOX);
					
					if (!drv.data.isParamModeDefined())
						drv.data.setParamMode(AbsSyn::modeType::BOX);
					
					if (!drv.data.isIterationSet())
					{
						yy::parser::error(@1, "Iteration number must be defined");
						YYERROR;
					}
				}
			symbolList matricesList footerList END
		{
			if (drv.data.getVarMode() == AbsSyn::modeType::BOX)
						drv.data.defaultDirections();
					
			if (drv.data.getParamMode() == AbsSyn::modeType::BOX)
				drv.data.defaultParamDirections();
			
			if (drv.data.getVarMode() != AbsSyn::modeType::POLY)
							drv.data.defaultTemplate();
			
			if (!drv.data.isTransModeDefined())
				drv.data.setTransMode(AbsSyn::transType::AFO);
			
			if (!drv.data.isAlphaDefined())
				drv.data.setAlpha(0.5);
		}
		| END { yy::parser::error(@1, "Empty file"); YYERROR; }

headerList	: header {}
						| headerList header {}

header			: PROB ":" problemType ";"
						{
							if (drv.data.isProblemDefined())
							{
								yy::parser::error(@4, "Problem has already been defined");
								YYERROR;
							}
							drv.data.setProblem($3);
						}
						| VARMODE ":" modeType ";"
						{
							if (drv.data.isVarModeDefined())
							{
								yy::parser::error(@4, "Variable modality has already been defined");
								YYERROR;
							}
							drv.data.setVarMode($3);
						}
						| PARAMMODE ":" modeType ";"
						{
							if (drv.data.isParamModeDefined())
							{
								yy::parser::error(@4, "Parameter modality has already been defined");
								YYERROR;
							}
							drv.data.setParamMode($3);
						}
						| ITER ":" INTEGER ";"
						{
							if (drv.data.isIterationSet())
							{
								yy::parser::error(@4, "Iteration number already defined");
								YYERROR;
							}
							
							drv.data.setIterations($3);
						}
						| PSPLITS ":" INTEGER ";"
						{
							if (drv.data.getMaxParameterSplits() > 0)
							{
								yy::parser::error(@4, "The maximum number of parameter splits has been already defined");
								YYERROR;
							}

							drv.data.setMaxParameterSplits($3);
						}
symbolList	: symbol {}
						| symbolList symbol {}

symbol			: VAR identList IN doubleInterval ";"
						{
							if (drv.data.getVarMode() != AbsSyn::modeType::BOX)
							{
								yy::parser::error(@$, "Cannot define variable bounds if modality is not 'boxes'");
								YYERROR;
							}
							
							for (unsigned i = 0; i < $2.size(); i++)
							{
								if (drv.data.isSymbolDefined($2[i]))
								{
									yy::parser::error(@2, "Symbol '" + $2[i] + "' already defined");
									YYERROR;
								}
								drv.data.addVariable(new AbsSyn::Variable($2[i]));
								drv.data.addBounds($4.first, $4.second);
							}
						}
						| VAR identList ";"
						{
							if (drv.data.getVarMode() == AbsSyn::modeType::BOX)
							{
								yy::parser::error(@$, "If variable modality is 'boxes', bounds must be provided");
								YYERROR;
							}
							
							for (unsigned i = 0; i < $2.size(); i++)
							{
								if (drv.data.isSymbolDefined($2[i]))
								{
									yy::parser::error(@2, "Symbol '" + $2[i] + "' already defined");
									YYERROR;
								}
								drv.data.addVariable(new AbsSyn::Variable($2[i]));
							}
						}
						| PARAM identList IN doubleInterval ";"
						{
							if (drv.data.getParamMode() != AbsSyn::modeType::BOX)
							{
								yy::parser::error(@$, "Cannot define parameter bounds if modality is not 'boxes'");
								YYERROR;
							}
							
							for (unsigned i = 0; i < $2.size(); i++)
							{
								if (drv.data.isSymbolDefined($2[i]))
								{
									yy::parser::error(@2, "Symbol '" + $2[i] + "' already defined");
									YYERROR;
								}
								drv.data.addParameter(new AbsSyn::Parameter($2[i]));
								drv.data.addParamBounds($4.first, $4.second);
							}
						}
						| PARAM identList ";"
						{
							for (unsigned i = 0; i < $2.size(); i++)
							{
								if (drv.data.isSymbolDefined($2[i]))
								{
									yy::parser::error(@2, "Symbol '" + $2[i] + "' already defined");
									YYERROR;
								}
								drv.data.addParameter(new AbsSyn::Parameter($2[i]));
							}
						}
						| CONST IDENT "=" expr ";"
						{
							if (!$4->isNumeric(drv.data))
							{
								yy::parser::error(@3, "Expression defining constant must be numeric");
								YYERROR;
							}
							if (drv.data.isSymbolDefined($2))
							{
								yy::parser::error(@2, "Symbol '" + $2 + "' already defined");
								YYERROR;
							}
							drv.data.addConstant(new AbsSyn::Constant($2, $4->evaluate(drv.data)));
							delete $4;
						}
						| DEFINE IDENT "=" expr ";"
						{
							if (drv.data.isSymbolDefined($2))
							{
								yy::parser::error(@2, "Symbol '" + $2 + "' already defined");
								YYERROR;
							}
							
							drv.data.addDefinition(new AbsSyn::Definition($2, $4));
						}
						| DYN "(" IDENT ")" "=" expr ";"
						{
							if (!drv.data.isSymbolDefined($3))
							{
								yy::parser::error(@3, "Unknown symbol name '" + $3 + "'");
								YYERROR;
							}
							
							if (!drv.data.isVarDefined($3))
							{
								yy::parser::error(@3, "'" + $3 + "' is not a variable");
								YYERROR;
							}
							
							AbsSyn::Variable *v = drv.data.getVar($3);
							if (v->isDynamicDefined())
							{
								yy::parser::error(@$, "Redefinition of dynamic for variable '" + $3 + "'");
								YYERROR;
							}
							v->setDynamic($6);
						}
						| SPEC ":" path_formula ";"
						{
							/*if (drv.data.getProblem() != AbsSyn::problemType::SYNTH)
							{
								yy::parser::error(@$, "Specification not required when problem is 'reachability'");
								YYERROR;
							}*/
							
							if (!$3->simplify())
							{
								yy::parser::error(@3, "Negations of UNTIL are not allowed");
								YYERROR;
							}
							
							drv.data.addSpec($3);
						}
						| SPEC ":" state_formula ";"
						{
							/*if (drv.data.getProblem() != AbsSyn::problemType::SYNTH)
							{
								yy::parser::error(@$, "Specification not required when problem is 'reachability'");
								YYERROR;
							}*/
							
							if (!$3->simplify())
							{
								yy::parser::error(@3, "Negations of UNTIL are not allowed");
								YYERROR;
							}
							
							drv.data.addSpec($3);
						}
						| ASSUME state_formula ";"
						{
							if (drv.data.getProblem() == AbsSyn::problemType::SYNTH) {
								yy::parser::error(@1, "Assumptions are not supported for synthesis yet");
								YYERROR;
							}
							
							if ($2->getType() != AbsSyn::Formula::formulaType::ATOM) {
								yy::parser::error(@2, "Only atomic formulas are supported");
								YYERROR;
							}
							if ($2->getEx()->hasParams(drv.data)) {
								yy::parser::error(@2, "Expressions in assertions cannot contain parameters");
								YYERROR;
							}
							if (!$2->isLinear(drv.data)) {
								yy::parser::error(@2, "Assumption must be linear");
								YYERROR;
							}
							
							AbsSyn::Assumption *a = new AbsSyn::Assumption($2->getEx());
							drv.data.addAssumption(a);
						}

matricesList	: %empty {}
							| matricesList matrices {}

matrices		: direction
						{
							/*if (drv.data.getVarMode() == AbsSyn::modeType::BOX)
							{
								yy::parser::error(@$, "Cannot define directions if variable modality is 'boxes'");
								YYERROR;
							}*/
							
							if (drv.data.templateRows() != 0)
							{
								yy::parser::error(@$, "Cannot define directions after defining template matrix");
								YYERROR;
							}
						}
						| template
						{
							/*if (drv.data.getVarMode() != AbsSyn::modeType::POLY)
							{
								yy::parser::error(@$, "Template matrix can be provided only if variable modality is 'polytopes'");
								YYERROR;
							}*/
							
							if (drv.data.directionsNum() == 0)
							{
								yy::parser::error(@$, "Template matrix cannot be provided before directions");
								YYERROR;
							}
						}
						| paramDir
						{
							/*if (drv.data.getParamMode() == AbsSyn::modeType::BOX)
							{
								yy::parser::error(@$, "Cannot define parameter directions if parameter modality is 'boxes'");
								YYERROR;
							}*/
						}

direction		: DIR "<" numList ">" IN doubleInterval ";" 
						{
							if (drv.data.getVarNum() != $3.size())
							{
								yy::parser::error(@3, "A direction must have as many components as the number of variables");
								YYERROR;
							}
							
							drv.data.addDirection($3, $6.first, $6.second);
						}

template		: TEMPL "=" _matrix
						{
							if ($3[0].size() != drv.data.getVarNum())
							{
								yy::parser::error(@3, "template matrix must have as many columns as the number of variables");
								YYERROR;
							}
							
							if (drv.data.templateRows() != 0)
							{
								yy::parser::error(@$, "Redefinition of template matrix");
								YYERROR;
							}
							
							for (unsigned i = 0; i < $3.size(); i++)
							{
								for (unsigned j = 0; j < $3[i].size(); j++)
								{
									if ($3[i][j] < 0 || (unsigned int)$3[i][j] >= drv.data.directionsNum())
									{
										yy::parser::error(@3, "Unknown direction");
										YYERROR;
									}
								}
							}
							
							drv.data.setTemplate($3);
						}

paramDir		: PDIR "<" numList ">" IN doubleInterval ";"
						{
							if (drv.data.getParamNum() != $3.size())
							{
								yy::parser::error(@3, "A paramter direction must have as many components as the number of parameters");
								YYERROR;
							}
							
							drv.data.addParamDirection($3, $6.first, $6.second);
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
									if (!$2->isNumeric(drv.data))
									{
										yy::parser::error(@2, "Intervals require numeric expressions");
										YYERROR;
									}
									
									if (!$4->isNumeric(drv.data))
									{
										yy::parser::error(@4, "Intervals require numeric expressions");
										YYERROR;
									}
									
									int x1 = (int) $2->evaluate(drv.data);
									int x2 = (int) $4->evaluate(drv.data);

									delete $2;
									delete $4;

									if (x2 < x1)
									{
										yy::parser::error(@$, "Right endpoint must be greater than or equal to the left one");
										YYERROR;
									}
									
									$$ = std::pair<int, int>(x1, x2);
								}

doubleInterval	: "[" expr "," expr "]"
								{
									if (!$2->isNumeric(drv.data))
									{
										yy::parser::error(@2, "Intervals require numeric expressions");
										YYERROR;
									}
									
									if (!$4->isNumeric(drv.data))
									{
										yy::parser::error(@4, "Intervals require numeric expressions");
										YYERROR;
									}

									double x1 = $2->evaluate(drv.data);
									double x2 = $4->evaluate(drv.data);

									delete $2;
									delete $4;

									if (x2 < x1)
									{
										yy::parser::error(@$, "Right endpoint must be greater than or equal to the left one");
										YYERROR;
									}
									
									$$ = std::pair<double, double>(x1, x2);
								}

number		: DOUBLE { $$ = $1; }
					| INTEGER { $$ = (double) $1; }

expr		: number	{ $$ = new AbsSyn::Expr($1); }
				| IDENT	
				{
					if (!drv.data.isSymbolDefined($1))
					{
						yy::parser::error(@1, "unknown symbol name");
						YYERROR;
					}
				
					if (drv.data.isConstDefined($1))
						$$ = new AbsSyn::Expr(drv.data.getConst($1)->getValue());
//					else if (drv.data.isDefDefined($1))
//						$$ = drv.data.getDef($1)->getValue()->copy();
					else
						$$ = new AbsSyn::Expr($1);
				}
				| expr "*" expr { $$ = $1->mul($3); }
				| expr "/" expr
				{
					if (!$3->isNumeric(drv.data))
					{
						yy::parser::error(@3, "cannot divide by non--numeric expression");
						YYERROR;
					}
					$$ = $1->div($3);
				}
				| expr "+" expr { $$ = $1->sum($3); }
				| expr "-" expr { $$ = $1->sub($3); }
				| "-" expr %prec UMINUS { $$ = $2->neg(); }
				| "(" expr ")" { $$ = $2;}

matrixRow	: "{" intList "}" { $$ = $2; }
					| "{" "}" { yy::parser::error(@$, "Matrix row cannot be empty"); YYERROR; }

numList		: expr
					{
						if (!$1->isNumeric(drv.data))
						{
							yy::parser::error(@1, "Expression must be numeric only");
							YYERROR;
						}
						$$ = std::vector<double>{$1->evaluate(drv.data)};
						delete $1;
					}
					| numList "," expr
					{
						if (!$3->isNumeric(drv.data))
						{
							yy::parser::error(@3, "Expression must be numeric only");
							YYERROR;
						}
						$1.push_back($3->evaluate(drv.data));
						$$ = $1;
						delete $3;
					}

intList		: INTEGER
					{
						if ($1 < 0 || (unsigned int)$1 >= drv.data.directionsNum())
						{
							yy::parser::error(@1, "Unknown direction");
							YYERROR;
						}
						$$ = std::vector<int>{$1};
					}
					| intList "," INTEGER
					{
						if ($3 < 0 || (unsigned int)$3 >= drv.data.directionsNum())
						{
							yy::parser::error(@3, "Unknown direction");
							YYERROR;
						}
						$1.push_back($3);
						$$ = $1;
					}

_matrix			: "{" rowList "}" { $$ = $2; }
						| "{" "}" { yy::parser::error(@$, "Matrix cannot be empty"); YYERROR; }

rowList			: matrixRow { $$ = std::vector<std::vector<int>>{$1}; }
						| rowList "," matrixRow { $1.push_back($3); $$ = $1; }

state_formula	: expr ">" expr { $$ = new AbsSyn::Formula($3->sub($1)); }
							| expr ">=" expr { $$ = new AbsSyn::Formula($3->sub($1)); }
							| expr "<" expr { $$ = new AbsSyn::Formula($1->sub($3)); }
							| expr "<=" expr { $$ = new AbsSyn::Formula($1->sub($3)); }
							| expr "=" expr
							{
								AbsSyn::Formula *f1 = new AbsSyn::Formula($1->sub($3));
								AbsSyn::Formula *f2 = new AbsSyn::Formula($3->sub($1));
								$$ = f1->conj(f2);
							}
							| state_formula AND state_formula		{ $$ = $1->conj($3); }
							| state_formula OR state_formula		{ $$ = $1->disj($3); }
							| "!" state_formula									{ $$ = $2->neg(); }
							| "(" state_formula ")" { $$ = $2; }

path_formula	: "G" intInterval state_formula %prec "G"	{ $$ = $3->always($2); }
							| "F" intInterval state_formula %prec "F"	{ $$ = $3->eventually($2); }
							| state_formula "U" intInterval state_formula %prec "U"	{ $$ = $1->until($3, $4); }
							| path_formula "&&" path_formula					{ $$ = $1->conj($3); }
							| path_formula "||" path_formula					{ $$ = $1->disj($3); }
							| "!" path_formula												{ $$ = $2->neg(); }
							| "(" path_formula ")" { $$ = $2; }

footerList	: %empty {}
						| footerList footer {}

footer	: OPT TRANS transType ";"
				{
					if (drv.data.isTransModeDefined())
					{
						yy::parser::error(@$, "Transformation type already defined");
						YYERROR;
					}
					drv.data.setTransMode($3);
				}
				| OPT DECOMP
				{
					if (drv.data.isDecompositionDefined())
					{
						yy::parser::error(@$, "Decomposition option already defined");
						YYERROR;
					}
					
					drv.data.setDecomposition();
				}
				| OPT ALPHA DOUBLE ";"
				{
					if (drv.data.isAlphaDefined())
					{
						yy::parser::error(@$, "Alpha already defined");
						YYERROR;
					}
					
					if ($3 > 1)
					{
						yy::parser::error(@3, "Alpha must be between 0 and 1");
						YYERROR;
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
