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

#define DO_NOT_DELETE_SUBFORMULAS 1<<0

class STL {
	const formula_type type;
protected:

	int options;

	STL(const formula_type type, const int options=0);

	inline bool delete_subformulas() const
	{
		return !(this->options&DO_NOT_DELETE_SUBFORMULAS);
	}
public:

	const STL& set_options(const int options);

	inline const formula_type& getType() const 
	{ 
		return type; 
	}

	virtual void print() const = 0;
	virtual ~STL();
};

#endif /* STL_H_ */
