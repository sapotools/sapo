#ifndef __DIRECTION_H__
#define __DIRECTION_H__

#include <cmath>

#include "Context.h"
#include "SymbolicAlgebra.h"
#include "Expr.h"

namespace AbsSyn
{

class Direction
{
  friend std::ostream &operator<<(std::ostream &os, const Direction &d);

public:
	enum Type
	{
		LT,		// <
		LE,		// <=
		GT,		// >
		GE,		// >=
		EQ,		// =
		INT		// lhs in [a,b]
	};
	
	Direction(SymbolicAlgebra::Expression<> e1, SymbolicAlgebra::Expression<> e2, Type t, double lb = -std::numeric_limits<double>::infinity(),
						double ub = std::numeric_limits<double>::infinity(), std::string dirName = ""):
				lhs(e1), rhs(e2), type(t), LB(lb), UB(ub), name(dirName) {}
	
	~Direction() {}
	
	std::vector<double> getDirectionVector(const Context &ctx, bool variables)
				const;	// returns the vector representing the direction
	
	double getOffset(const Context &ctx)
				const;		// return the offset of the direction, if type is not INT
	
	std::string getName() const
	{
		return name;
	}
	void setName(std::string n)
	{
		name = n;
	}
	
	double getLB(const Context &ctx) const;
	double getUB(const Context &ctx) const;
	
	bool hasLB() const
	{
		return type == Type::INT || type == Type::EQ || LB != -std::numeric_limits<double>::infinity();
	}
	bool hasUB() const
	{
		return type == Type::INT || type == Type::EQ || UB != std::numeric_limits<double>::infinity();
	}
	
	void setLB(const Context &ctx, double val);
	void setUB(const Context &ctx, double val);
	
	SymbolicAlgebra::Expression<> getLHS()
	{
		return lhs;
	}
	SymbolicAlgebra::Expression<> getRHS()
	{
		return rhs;
	}
	
	Type getType() const
	{
		return type;
	}
	
	// checks if the direction contains variable names
	bool hasVars(const Context &ctx) const
	{
		return AbsSyn::hasVars(lhs, ctx) || (type != Type::INT && AbsSyn::hasVars(rhs, ctx));
	}
	
	// checks if the direction contains parameter names
	bool hasParams(const Context &ctx) const
	{
		return AbsSyn::hasParams(lhs, ctx) || (type != Type::INT && AbsSyn::hasParams(rhs, ctx));
	}
	
	Direction *copy() const;		// deep copy of direction
	
	Direction *getComplementary() const;	// returns the negated direction
	
	bool compare(Direction *d, const Context &ctx, bool variable = true) const;			// comparison between directions
	
	bool covers(const Context &ctx, const Symbol<> &s)
				const;	// checks if the symbol named "name" is present in the direction
	
protected:
	SymbolicAlgebra::Expression<> lhs, rhs;
	Type type;
	double LB, UB;		// used only if type is INT
	std::string name;
};

}

#endif
