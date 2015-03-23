/*
 * Box.h
 *
 *  Created on: Nov 3, 2014
 *      Author: dreossi
 */

#ifndef BOX_H_
#define BOX_H_

#include "Common.h"
#include "Polyhedron.h"

class Box : public Polyhedron{
private:

public:

	Box(vector<lst> vars);
	LinearSystem* gen2const(vector<double> q, vector<double> beta);
	poly_values const2gen(LinearSystem *constr);
	virtual ~Box();
};

#endif /* BOX_H_ */
