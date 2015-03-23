/*
 * LinearSystem.h
 *
 *  Created on: Oct 24, 2014
 *      Author: dreossi
 */

#ifndef LINEARSYSTEM_H_
#define LINEARSYSTEM_H_

#include "Common.h"
#include <glpk.h>

class LinearSystem {
private:
	int n_vars;			// number of variables
	lst vars;			// 	list of variables
	lst constraints;	//	list of constraints
	vector< vector<double> > A; // matrix A
	vector< double > b; 		// vector b

	void initLS();				// initialize A and b
	double solveLinearSystem(vector< vector< double > > A, vector< double > b, vector< double > obj_fun, int min_max);


public:
	LinearSystem(vector< vector<double> > A, vector< double > b);
	LinearSystem(lst vars, lst constraints);

	vector< vector<double> > getA(); //return A
	vector<double> getb(); 			//return b
	double getA(int i, int j); 		//return A[i][j]
	double getb(int i); 			//return b[i]

	double minLinearSystem(lst vars, ex obj_fun);
	double maxLinearSystem(lst vars, ex obj_fun);
	bool isEmpty();									// determine this LS is empty
	LinearSystem* appendLinearSystem(LinearSystem *LS);

	void print();


	virtual ~LinearSystem();
};

#endif /* LINEARSYSTEM_H_ */
