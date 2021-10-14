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

	bool isIn(vector< double > Ai, const double bi) const;	// check if a constraint is already in
	void initLS();								// initialize A and b

public:

	LinearSystem();
	LinearSystem(vector< vector<double> > A, vector< double > b);
	LinearSystem(lst vars, lst constraints);

	/**
	 * Return the template matrix
	 *
	 * @return template matrix
	 */
	inline const vector< vector<double> >& getA() const { return this->A; }

	/**
	 * Return the offset vector
	 *
	 * @return offset vector
	 */
	inline const vector<double>& getb() const { return this->b; }

	const double& getA(int i, int j) const;
	const double& getb(int i) const;

	// optimization functions
	double minLinearSystem(const lst& vars, const ex& obj_fun) const;
	double maxLinearSystem(const lst& vars, const ex& obj_fun) const;
	double maxLinearSystem(const vector< double >& obj_fun_coeffs) const;
	bool isEmpty(const bool strict_inequality=false) const;						// determine this LS is empty
	bool operator==(const LinearSystem& LS) const;  // decide whether two systems are equivalent
	bool solutionsAlsoSatisfy(const LinearSystem& LS) const;

	// operations on linear system
	LinearSystem* intersectWith(const LinearSystem *LS) const;
	vector<bool> redundantCons() const;
	LinearSystem *simplify();

	/**
	 * Return the number of variables
	 */
	inline int dim() const { if(this->isEmpty()){ return 0; } return this->A[0].size(); };

	/**
	 * Return the number of inequalities
	 */
	inline int size() const { return this->b.size(); };

	double volBoundingBox();

	void print() const;
	void plotRegion() const;
	void plotRegionToFile(const char *file_name, const char color) const;
	void plotRegionT(const double t) const;
	void plotRegion(const vector<int>& rows, const vector<int>& cols) const;

	~LinearSystem();
};

#endif /* LINEARSYSTEM_H_ */
