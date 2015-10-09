/*
 * MultiParallelotope.h
 *
 *  Created on: Oct 30, 2014
 *      Author: Tommaso Dreossi
 */

#ifndef MULTIPARALLELOTOPE_H_
#define MULTIPARALLELOTOPE_H_

#include "float.h"
#include "Common.h"
#include "Polyhedron.h"
#include "Parallelotope.h"
#include "LinearSystem.h"

class MultiParallelotope : public Polyhedron{

private:

	lst qs, as, bs;							// parallelotope variables
	vector< lst > us;
	lst generator_funs;						// generator functions

	vector< vector< double > > directions;	// direction of facets
	vector< double > upper_bounds;			// upper and lower bounds
	vector< double > lower_bounds;
	vector< vector< int > > templates;		// parallelotope templates

	vector< vector< double > > angles;

	pair< vector<int>,  vector<bool> > pickNdirections(vector< vector<double> > directions, vector<bool> marks);

	bool allMarked(vector<bool> v);
	pair< int,int > minAngleBothNotMarked( vector<double> angles, vector< vector< int > > angle_idx, vector<bool> marks);
	int minAngleBothNotMarkedWith( vector<double> angles, vector< vector< int > > angle_idx, vector<bool> marks, vector<int> tmp_lambda);
	vector<int> fillLambda( vector<double> angles, vector< vector< int > > angle_idx, vector<int> tmp_lambda, int dim_lambda);


	double norm(vector<double> v);
	double prod(vector<double> v1, vector<double> v2);
	double angle(vector<double> v1, vector<double> v2);
	bool isIn(int e, vector<int> v);

	void initAngles();

public:

	MultiParallelotope(lst qs, lst as, lst bs, vector< lst > us);
	MultiParallelotope(vector< vector< double > > directions, vector< double > upper_bounds, vector< double > lower_bounds, vector< vector< int > > templates);

	vector< vector< int > > coupleDirectionsAngle(vector< vector< double > > dirs);		// couple the directions
	pair< vector< double >, vector< double > > shrink(LinearSystem *H);					// shrink current set of parallelotopes
	void decompose(LinearSystem *H, vector< vector< double > > directions);
	void decompose(vector<int> hp_idx, vector<bool> offset);

	LinearSystem* parallelotopesIntersection();

	Parallelotope* getPar(int p);
	LinearSystem* getParallelotope(int p);
	lst getParallelotopeGenFun(int p);

	lst getBaseVertexVars(){ return this->qs; };
	lst getFreeVars(){ return this->as; };
	lst getLengthsVars(){ return this->bs; };
	vector< lst > getGenVars(){ return this->us; };
	lst getGenVars(int i){ return this->us[i]; };
	int getDim(){ return this->directions[0].size(); };
	int getCard();
	int getNumDirs();
	vector< vector< double > > getDirections(){ return this->directions; };
	vector< vector< int > > getTemplates(){ return this->templates; };
	vector< double > getDirection(int i);
	vector< int > getTemplate(int i);
	ex getGenFuns(int i){ return this->generator_funs[i]; }
	lst getGenFuns(){ return this->generator_funs; }

	void setUpperBounds(vector< double > upper_bounds);
	void setLowerBounds(vector< double > lower_bounds);
	void setBounds(vector< double > upper_bounds, vector< double > lower_bounds);

	vector< double > negate(vector< double > v);

	void print();
	void printParallelotopes();
	void printParallelotopesIntersection();

	virtual ~MultiParallelotope();
};

#endif /* PARALLELOTOPE_H_ */
