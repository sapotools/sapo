/*
 * Polyhedron.h
 *
 *  Created on: Oct 30, 2014
 *      Author: dreossi
 */

#ifndef POLYHEDRON_H_
#define POLYHEDRON_H_

#include "Common.h"
#include "LinearSystem.h"

struct poly_values{
	vector<double> base_vertex;
	vector< vector<double> > versors;
	vector<double> lenghts;
};

class Polyhedron {

protected:

	int dim;							// dimension of the parallelotope
	vector<lst> vars;					// variables appearing in generato function
										// vars[0] q: base vertex
										// vars[1] alpha : free variables \in [0,1]
										// vars[2] beta : generator amplitudes
	lst generator_function;				// generator function
	vector< vector<double> > u;			// versors
	vector< vector<double> > template_matrix;	// Template matrix

public:

	Polyhedron();
	lst getGeneratorFunction();
	lst getQ();
	lst getAlpha();
	lst getBeta();
	int getDim();
	vector< vector< double > > getTemplate();
	virtual LinearSystem* gen2const(vector<double> q, vector<double> beta){ return 0; }
	virtual poly_values const2gen(LinearSystem *constr){ poly_values p; return p; };

	virtual int getMaxSize(){ return 0; };
	virtual vector<ex> getConvCombs(int i){ vector<ex> c; return c; };

	virtual ~Polyhedron();
};

#endif /* POLYHEDRON_H_ */
