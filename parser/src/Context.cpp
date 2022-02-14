#include "Context.h"

using namespace AbsSyn;

bool Context::isVarDefined(std::string name) const
{
	for (unsigned i = 0; i < vars.size(); i++) {
		if (vars[i].get_symbol_name(vars[i].get_id()) == name) {
			return true;
		}
	}
	return false;
}

bool Context::isParamDefined(std::string name) const
{
	for (unsigned i = 0; i < params.size(); i++) {
		if (params[i].get_symbol_name(params[i].get_id()) == name) {
			return true;
		}
	}
	return false;
}

bool Context::isConstDefined(std::string name) const
{
	for (auto it = consts.begin(); it != consts.end(); it++) {
		if (it->first.get_symbol_name(it->first.get_id()) == name) {
			return true;
		}
	}
	return false;
}

bool Context::isDefDefined(std::string name) const
{
	for (auto it = defines.begin(); it != defines.end(); it++) {
		if (it->first.get_symbol_name(it->first.get_id()) == name) {
			return true;
		}
	}
	return false;
}

bool Context::isSymbolDefined(std::string name) const
{
	return isVarDefined(name) || isParamDefined(name) || isConstDefined(name) || isDefDefined(name);
}


SymbolicAlgebra::Symbol<> Context::getSymbol(std::string name) const
{
	for (unsigned i = 0; i < vars.size(); i++) {
		if (vars[i].get_symbol_name(vars[i].get_id()) == name) {
			return vars[i];
		}
	}
	
	for (unsigned i = 0; i < params.size(); i++) {
		if (params[i].get_symbol_name(params[i].get_id()) == name) {
			return params[i];
		}
	}
	
	for (auto it = consts.begin(); it != consts.end(); it++) {
		if (it->first.get_symbol_name(it->first.get_id()) == name) {
			return it->first;
		}
	}
	
	for (auto it = defines.begin(); it != defines.end(); it++) {
		if (it->first.get_symbol_name(it->first.get_id()) == name) {
			return it->first;
		}
	}
	
	throw std::logic_error("symbol " + name + " is not defined");
}



SymbolicAlgebra::Expression<>::replacement_type Context::getConsts()
{
	SymbolicAlgebra::Expression<>::replacement_type rep{};
	
	for (auto it = consts.begin(); it != consts.end(); it++) {
		rep[it->first] = it->second;
	}
	
	return rep;
}
