/**
 * @file Bundle.h
 * Represent and manipulate bundles of parallelotopes whose intersection
 * represents a polytope
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#ifndef BUNDLE_H_
#define BUNDLE_H_

#include "float.h"
#include "Common.h"
#include "BaseConverter.h"
#include "Parallelotope.h"
#include "LinearSystem.h"
#include "VarsGenerator.h"
#include <cmath>

class Bundle {

private:
	int dim;							// dimension
	vector< vector< double > > L;		// direction matrix
	vector< double > offp;				// superior offset
	vector< double > offm;				// inferior offset
	vector< vector< int > > T;			// templates matrix
	vector< vector< double > > Theta;	// matrix of orthogonal proximity
	vector<lst> vars;					// variables appearing in generato function
										// vars[0] q: base vertex
										// vars[1] alpha : free variables \in [0,1]
										// vars[2] beta : generator amplitudes

	// map with Bernstein coefficients
	map< vector<int>, lst > bernCoeffs;

	double initTheta();
	vector< double > offsetDistances();

	// operations on vectors
	double norm(vector<double> v);
	double prod(vector<double> v1, vector<double> v2);
	double angle(vector<double> v1, vector<double> v2);
	double orthProx(vector<double> v1, vector<double> v2);
	double maxOrthProx(int vIdx, vector<int> dirsIdx);
	double maxOrthProx(vector<int> dirsIdx);
	double maxOrthProx(vector< vector<int> > T);
	double maxOffsetDist(int vIdx, vector<int> dirsIdx, vector<double> dists);
	double maxOffsetDist(vector<int> dirsIdx, vector<double> dists);
	double maxOffsetDist(vector< vector<int> > T, vector<double> dists);
	
	//vector< double > negate(vector< double > v);
	
	bool isIn(int n, vector<int> v);
	bool isIn(vector<int> v, vector< vector< int > > vlist);
	bool isPermutation(vector<int> v1, vector<int> v2);


	bool validTemp(vector< vector<int> > T, int card, vector<int> dirs);	// check if a template is valid
	vector<lst> transformContrPts(lst vars, lst f, int mode);

public:

	// constructors
	Bundle(vector< vector< double > > L, vector< double > offp, vector< double > offm, 	vector< vector< int > > T);
	Bundle(vector<lst> vars, vector< vector< double > > L, vector< double > offp, vector< double > offm, 	vector< vector< int > > T);

	int getDim(){ return this->dim; };
	int getSize(){ return L.size(); };
	int getCard(){ return T.size(); };
	int getNumDirs(){ return this->L.size(); };

	vector<int> getTemplate(int i){ return this->T[i]; };
	double getOffp(int i){ return this->offp[i]; };
	double getOffm(int i){ return this->offm[i]; };
	LinearSystem getBundle();
	Parallelotope* getParallelotope(int i);

	void setTemplate(vector< vector< int > > T);
	void setOffsetP(vector< double > offp){ this->offp = offp; }
	void setOffsetM(vector< double > offm){ this->offm = offm; }

	// operations on bundles
	Bundle* canonize();
	Bundle* decompose(double alpha, int max_iters);
	Bundle* transform(lst vars, lst f, map< vector<int>,pair<lst,lst> > &controlPts, int mode);
	Bundle* transform(lst vars, lst params, lst f, LinearSystem *paraSet, map< vector<int>,pair<lst,lst> > &controlPts, int mode);

	virtual ~Bundle();
};

#endif /* BUNDLE_H_ */
