#ifndef __DEFINITION_H__
#define __DEFINITION_H__

#include "SymbolicAlgebra.h"

namespace AbsSyn
{

class Definition
{
	friend std::ostream &operator<<(std::ostream &os, const Definition &d)
	{
		return os << d.s;
	}

public:
	Definition(SymbolicAlgebra::Symbol<> sym, SymbolicAlgebra::Expression<> e): s(sym), value(e) {}

	~Definition() {}

	std::string getName()
	{
		return s.get_symbol_name(s.get_id());
	}
	
	SymbolicAlgebra::Symbol<> getSymbol() const
	{
		return s;
	}
	
	SymbolicAlgebra::Expression<> getValue()
	{
		return value;
	}
	const SymbolicAlgebra::Expression<> getValue() const
	{
		return value;
	}

private:
	SymbolicAlgebra::Symbol<> s;
	SymbolicAlgebra::Expression<> value;
};

}

#endif
