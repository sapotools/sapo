/*
 * MultiReacher.h
 *
 *  Created on: Jun 3, 2015
 *      Author: dreossi
 */

#ifndef MULTIREACHER_H_
#define MULTIREACHER_H_

#include "MultiParallelotope.h"
#include "LinearSystem.h"
#include "Parallelotope.h"
#include "BaseConverter.h"
#include "float.h"
#include <algorithm>



class MultiReacher {

private:
	lst vars, dynamics;
	lst qs, as, bs, ls;
	vector< lst > us;
	lst dyn_bern_coeffs;
	vector< vector<lst> > static_bern_coeffs;
	vector< vector<ex> > Lambda_fogs;

public:
	MultiReacher(lst vars, lst dynamics, MultiParallelotope *D, lst ls);


	MultiParallelotope* reachStepStaticTemplate(MultiParallelotope *D);
	MultiParallelotope* reachStepDynTemplate(MultiParallelotope *D);
	MultiParallelotope* reachStepNumerical(MultiParallelotope *D);

	double maxBernCoeffDirs(vector<double> base_vertex, vector<double> lenghts, vector< vector<double> > generators, vector<double> temp_direction);
	double maxBernCoeffStatic(int parall_idx, int dir_idx, vector<double> base_vertex, vector<double> lenghts);
	double maxBernCoeffNumerical(int parall_idx, int dir_idx, vector<double> base_vertex, vector<double> lenghts);


	void reach(MultiParallelotope *D, int k);
	void staticReach(MultiParallelotope *D, int k);
	void numericalReach(MultiParallelotope *D, int k);

	vector<double> generateRndDir(int n);
	vector< vector<double> > generateRndDirs(int m, int n);

	double norm(vector<double> v);
	vector<double> normalize(vector<double> v);
	vector< double > negate(vector< double > v);

	virtual ~MultiReacher();
};

#endif /* MULTIREACHER_H_ */
