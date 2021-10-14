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

class Eventually : public STL {

private:

	const std::shared_ptr<STL> f;			// subformula
	int a, b;			// interval bounds

public:

	Eventually(int a, int b, const std::shared_ptr<STL> f);

	inline const std::shared_ptr<STL> getSubFormula() const { return f; }

	inline int getA() const { return a; }
	inline int getB() const { return b; }

	void print() const;

	~Eventually();
};

#endif /* Eventually_H */
