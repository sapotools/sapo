#ifndef __CONTEXT_H__
#define __CONTEXT_H__

#include <vector>
#include <map>

#include "SymbolicAlgebra.h"

namespace AbsSyn
{

class Context
{
public:
	Context():
		vars(0), params(0), consts{}, defines{}
		{}
	
	~Context() {}
	
	bool isVarDefined(std::string name) const;
	inline bool isVarDefined(SymbolicAlgebra::Symbol<> s) const
	{
		return isVarDefined(s.get_symbol_name(s.get_id()));
	}
	
	bool isParamDefined(std::string name) const;
	inline bool isParamDefined(SymbolicAlgebra::Symbol<> s) const
	{
		return isParamDefined(s.get_symbol_name(s.get_id()));
	}
	
	bool isConstDefined(std::string name) const;
	inline bool isConstDefined(SymbolicAlgebra::Symbol<> s) const
	{
		return isConstDefined(s.get_symbol_name(s.get_id()));
	}
	
	bool isDefDefined(std::string name) const;
	inline bool isDefDefined(SymbolicAlgebra::Symbol<> s) const
	{
		return isDefDefined(s.get_symbol_name(s.get_id()));
	}
	
	bool isSymbolDefined(std::string name) const;
	inline bool isSymbolDefined(SymbolicAlgebra::Symbol<> s) const
	{
		return isSymbolDefined(s.get_symbol_name(s.get_id()));
	}
	
	SymbolicAlgebra::Symbol<> getSymbol(std::string name) const;
	
	unsigned getVarNum() const
	{
		return vars.size();
	}
	SymbolicAlgebra::Symbol<> getVar(unsigned i) const
	{
		return vars[i];
	}
	void addVariable(SymbolicAlgebra::Symbol<> s)
	{
		vars.push_back(s);
	}
	std::vector<SymbolicAlgebra::Symbol<>> getVars() const
	{
		return vars;
	}
	
	unsigned getParamNum() const
	{
		return params.size();
	}
	SymbolicAlgebra::Symbol<> getParam(unsigned i) const
	{
		return params[i];
	}
	void addParameter(SymbolicAlgebra::Symbol<> s)
	{
		params.push_back(s);
	}
	std::vector<SymbolicAlgebra::Symbol<>> getParams() const
	{
		return params;
	}
	
	unsigned getConstantNum() const
	{
		return consts.size();
	}
	double getConst(SymbolicAlgebra::Symbol<> s) const
	{
		return consts.at(s);
	}
	void addConstant(SymbolicAlgebra::Symbol<> s, double val)
	{
		consts[s] = val;
	}
	SymbolicAlgebra::Expression<>::replacement_type getConsts();
	
	unsigned getDefineNum() const
	{
		return defines.size();
	}
	void addDefinition(SymbolicAlgebra::Symbol<> s, SymbolicAlgebra::Expression<> e)
	{
		defines[s] = e;
	}
	
	SymbolicAlgebra::Expression<> getDefineVal(SymbolicAlgebra::Symbol<> s)
	{
		return defines[s];
	}
	
	SymbolicAlgebra::Expression<>::replacement_type getDefines()
	{
		return defines;
	}
	
protected:
	std::vector<SymbolicAlgebra::Symbol<>> vars;
	std::vector<SymbolicAlgebra::Symbol<>> params;
	std::map<SymbolicAlgebra::Symbol<>, double> consts;
	std::map<SymbolicAlgebra::Symbol<>, SymbolicAlgebra::Expression<>> defines;
};

}

#endif
