/**
 * @file STL.h
 * STL formula
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#ifndef STL_H_
#define STL_H_

#include "Common.h"
#include <string.h>

class STL {
public:
	inline const formula_type& getType() const { return type; }

/*
	virtual inline const ex& getPredicate() const { throw std::logic_error("This object is not equipped with this method."); };
	virtual inline STL * getSubFormula() { throw std::logic_error("This object is not equipped with this method."); }
	virtual inline STL * getLeftSubFormula() { throw std::logic_error("This object is not equipped with this method."); }
	virtual inline STL * getRightSubFormula() { throw std::logic_error("This object is not equipped with this method."); }
	virtual inline void setPredicateControlPts(const vector<lst>& predicateControlPts) { throw std::logic_error("This object is not equipped with this method."); }
	virtual inline const vector<lst>& getPredicateControlPts() const { throw std::logic_error("This object is not equipped with this method."); }

	virtual inline const int getA() const { throw std::logic_error("This object is not equipped with this method."); }
	virtual inline const int getB() const { throw std::logic_error("This object is not equipped with this method."); }
	virtual inline void setA(const int a) { throw std::logic_error("This object is not equipped with this method."); }
	virtual inline void setB(const int b) { throw std::logic_error("This object is not equipped with this method."); }
	virtual inline int getID() const { throw std::logic_error("This object is not equipped with this method."); }
*/

	virtual void print() const = 0;

protected:
	STL(formula_type type);

	formula_type type;
};

#endif /* STL_H_ */
