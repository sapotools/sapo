/**
 * @file LinearSystem.h
 * Represent and manipulate a linear system
 * It can be used to represent polytopes (reached states, parameters, etc.)
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#ifndef LINEARSYSTEM_H_
#define LINEARSYSTEM_H_

#include "Common.h"
#include <glpk.h>

#include <iostream>
#include <fstream>

class LinearSystem {

private:
	int n_vars;					// number of variables
	lst vars;					// 	list of variables
	lst constraints;			//	list of constraints
	vector< vector<double> > A; // matrix A
	vector< double > b; 		// vector b

	bool isIn(vector< double > Ai, double bi);	// check if a constraint is already in
	void initLS();								// initialize A and b
	double solveLinearSystem(vector< vector< double > > A, vector< double > b, vector< double > obj_fun, int min_max);
	bool zeroLine(vector<double> line);


public:

	LinearSystem();
	LinearSystem(vector< vector<double> > A, vector< double > b);
	LinearSystem(lst vars, lst constraints);

	vector< vector<double> > getA(); 	//return A
	vector<double> getb(); 				//return b
	double getA(int i, int j); 			//return A[i][j]
	double getb(int i); 				//return b[i]

	// optimization functions
	double minLinearSystem(lst vars, ex obj_fun);
	double maxLinearSystem(lst vars, ex obj_fun);
	double maxLinearSystem(vector< double > obj_fun_coeffs);
	bool isEmpty();						// determine this LS is empty

	// operations on linear system
	LinearSystem* appendLinearSystem(LinearSystem *LS);
	vector<bool> redundantCons();
	LinearSystem *simplify();

	int dim(){ if(!this->isEmpty()){ return this->A[0].size(); }else{ return 0;}	};
	int size(){ return this->b.size(); };

	double volBoundingBox();

	void print();
	void plotRegion();
	void plotRegionToFile(char *file_name, char color);
	void plotRegionT(double t);
	void plotRegion(vector<int> rows, vector<int> cols);

	virtual ~LinearSystem();
};

#endif /* LINEARSYSTEM_H_ */
