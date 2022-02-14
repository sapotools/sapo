#ifndef __VARIABLE_H__
#define __VARIABLE_H__

#include "SymbolicAlgebra.h"

namespace AbsSyn
{

class Variable
{
	friend std::ostream &operator<<(std::ostream &os, Variable &v)
	{
		return os << v.s;
	}

public:
	Variable(const SymbolicAlgebra::Symbol<> &n): s(n), dynamic(0), covered(false), dynamicDefined(false) {}

	~Variable() {}

	const std::string &getName() const
	{
		return s.get_symbol_name(s.get_id());
	}
	
	SymbolicAlgebra::Symbol<> getSymbol() const
	{
		return s;
	}

	SymbolicAlgebra::Expression<> getDynamic() const
	{
		return dynamic;
	}

	void setDynamic(SymbolicAlgebra::Expression<> e)
	{
		dynamic = e;
		dynamicDefined = true;
	}

	// checks if dynamic has already been set
	bool isDynamicDefined()
	{
		return dynamicDefined;
	}
	
	bool isCovered()
	{
		return covered;
	}
	
	void setCovered()
	{
		this->covered = true;
	}

protected:
	SymbolicAlgebra::Symbol<> s;
	SymbolicAlgebra::Expression<> dynamic;
	bool covered;
	bool dynamicDefined;
};

}

#endif
