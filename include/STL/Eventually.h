/**
 * @file Eventually.h
 * Eventually STL formula
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#ifndef EVENTUALLY_H_
#define EVENTUALLY_H_

#include "STL.h"
#include "TimeInterval.h"

class Eventually : public STL {

private:

	const std::shared_ptr<STL> f;	// subformula
	TimeInterval t_itvl;			// temporal formula bounds

public:

	Eventually(int a, int b, const std::shared_ptr<STL> f);

	inline const std::shared_ptr<STL> getSubFormula() const { return f; }

	inline const int& getA() const { return t_itvl.begin(); }
	inline const int& getB() const { return t_itvl.end(); }

	void print() const;

	~Eventually();
};

#endif /* Eventually_H */
