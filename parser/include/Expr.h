#ifndef __EXPR_H__
#define __EXPR_H__

#include <vector>
#include <algorithm>
#include <set>

#include "SymbolicAlgebra.h"
#include "Context.h"

using namespace SymbolicAlgebra;
using namespace std;

namespace AbsSyn
{

// returns an equivalent expression, where each definition name is replaced by its expression
inline Expression<> expandDefinitions(const Expression<> &e, Context ctx)
{
	Expression<> res(e);
	res.replace(ctx.getDefines());
	res.replace(ctx.getConsts());
	res.expand();
	return res;
}

// checks if the expression contains at least one of the symbols
inline bool contains(const Expression<> &e, Context ctx, const vector<Symbol<>> &symbols)
{
	std::set<Symbol<>> ids = e.get_symbols();
	
	for (auto it = ids.begin(); it != ids.end(); it++) {
		if (ctx.isDefDefined(*it) && contains(ctx.getDefineVal(*it), ctx, symbols)) {
			return true;
		} else {
			for (auto sym = symbols.begin(); sym != symbols.end(); sym++) {
				if (it->get_symbol_name(it->get_id()) == sym->get_symbol_name(sym->get_id())) {
					return true;
				}
			}
		}
	}
	return false;
}

// checks if the expression contains variables
inline bool hasVars(const Expression<> &e, Context ctx)
{
	return contains(e, ctx, ctx.getVars());
}

// checks if the expression contains parameters
inline bool hasParams(const Expression<> &e, Context ctx)
{
	return contains(e, ctx, ctx.getParams());
}

// returns the degree of the polynomial expression, considering only variables
inline unsigned getVarDegree(const Expression<> &e, Context ctx)
{
	Expression<> temp = expandDefinitions(e, ctx);
	
	unsigned degree = 0;
	
	for (unsigned i = 0; i < ctx.getVarNum(); i++) {
		unsigned d = temp.degree(ctx.getVar(i));
		if (d > degree) {
			degree = d;
		}
	}
	
	return degree;
}

// returns the degree of the polynomial expression, considering only parameters
inline unsigned getParamDegree(const Expression<> &e, Context ctx)
{
	Expression<> temp = expandDefinitions(e, ctx);
	
	unsigned degree = 0;
	
	for (unsigned i = 0; i < ctx.getParamNum(); i++) {
		unsigned d = temp.degree(ctx.getParam(i));
		if (d > degree) {
			degree = d;
		}
	}
	
	return degree;
}

// checks if the expression is numeric (no vars or parameters)
inline bool isNumeric(const Expression<> &e, Context ctx)
{
	std::set<Symbol<>> ids = e.get_symbols();
	
	for (auto it = ids.begin(); it != ids.end(); it++) {
		if (ctx.isVarDefined(*it) || ctx.isParamDefined(*it)) {
			return false;
		} else if (ctx.isDefDefined(*it) && !isNumeric(ctx.getDefineVal(*it), ctx)) {
			return false;
		}
	}
	
	return true;
}

// evaluates a numeric expression
inline double evaluate(const Expression<> &e, Context ctx)
{
	Expression<> temp = expandDefinitions(e, ctx);
	
	return temp.evaluate<double>();
}

/*
 * returns the coefficient of the symbol "name"
 * in the simplified expression (in which each symbol
 * apperas only one time)
 * Applies only to linear expressions w.r.t the symbol specified
 */
inline double getCoefficient(const Expression<> &e, Context ctx, const Symbol<> &s)
{
	Expression<> temp = expandDefinitions(e, ctx);
	
	return temp.get_coeff(s, 1).evaluate<double>();
}

/*
 * returns the numerical term
 * in the simplified expression (in which each symbol
 * apperas only one time)
 */
inline double getOffset(const Expression<> &e, Context ctx)
{
	Expression<> temp = expandDefinitions(e, ctx);

	Expression<>::replacement_type rep{};
	
	for (unsigned i = 0; i < ctx.getVarNum(); i++) {
		rep[ctx.getVar(i)] = 0;
	}
	for (unsigned i = 0; i < ctx.getParamNum(); i++) {
		rep[ctx.getParam(i)] = 0;
	}
	
	temp.replace(rep);
	
	return temp.evaluate<double>();
}


}

#endif
