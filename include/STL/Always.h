/**
 * @file Always.h
 * Always STL formula
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#ifndef ALWAYS_H_
#define ALWAYS_H_

#include "STL.h"
#include "TimeInterval.h"

class Always : public STL {

private:

	const std::shared_ptr<STL> f;	// subformula
	TimeInterval t_itvl;			// temporal formula bounds

public:
	Always(int a, int b, const std::shared_ptr<STL> f);

	inline const std::shared_ptr<STL> getSubFormula() const { return this->f; }

	inline const int& getA() const { return t_itvl.begin(); }
	inline const int& getB() const { return t_itvl.end(); }

	void print() const;

	virtual ~Always();
};

#endif /* ALWAYS_H */
