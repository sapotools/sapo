/**
 * @file Until.h
 * Until STL formula
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#ifndef UNTIL_H_
#define UNTIL_H_

#include "STL.h"
#include "TimeInterval.h"

class Until : public STL {

private:

	const std::shared_ptr<STL> f1, f2;		// subformulas
	TimeInterval t_itvl;					// temporal formula bounds

public:

	Until(const std::shared_ptr<STL> f1, int a, int b, const std::shared_ptr<STL> f2);

	inline const std::shared_ptr<STL> getLeftSubFormula() const {return f1;}
	inline const std::shared_ptr<STL> getRightSubFormula() const {return f2;}

	inline const int& getA() const { return t_itvl.begin(); }
	inline const int& getB() const { return t_itvl.end(); }

	void print() const;

	~Until();
};

#endif /* UNTIL_H */
