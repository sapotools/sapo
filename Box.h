/*
 * Box.h
 *
 *  Created on: Nov 3, 2014
 *      Author: Tommaso Dreossi
 */

#ifndef BOX_H_
#define BOX_H_

#include "Common.h"
#include "Polyhedron.h"
#include "LinearSystem.h"


class Box : public Polyhedron{

	private:

	vector< vector<double> > A;
	vector< double > b;

	public:

		Box(vector< vector<double> > A, vector< double > b);

		double volume();

		void print(){ };
		virtual ~Box();
	};

#endif /* BOX_H_ */
