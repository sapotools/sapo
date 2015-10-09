/*
 * Bundle.h
 *
 *  Created on: Sep 18, 2015
 *      Author: dreossi
 */

#ifndef BUNDLE_H_
#define BUNDLE_H_

#include "float.h"
#include "Common.h"
#include "BaseConverter.h"
#include "Parallelotope.h"
#include "LinearSystem.h"
#include "CoeffMap.h"

class Bundle {

private:
	int dim;							// dimension
	vector< vector< double > > L;		// direction matrix
	vector< double > offp;				// superior offset
	vector< double > offm;				// inferior offset
	vector< vector< int > > T;		// templates matrix
	vector< vector< double > > Theta;	// matrix of orthogonal proximity
	vector<lst> vars;					// variables appearing in generato function
										// vars[0] q: base vertex
										// vars[1] alpha : free variables \in [0,1]
										// vars[2] beta : generator amplitudes
	//CoeffMap *bernCoeffs;	// map with Bernstein coefficients
	map< vector<int>, lst, classcomp > bernCoeffs;
	CoeffMap *coeffMap;

	double initTheta();

	vector< double > offsetDistances();

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
	vector< double > negate(vector< double > v);

	bool isIn(int n, vector<int> v);
	bool isIn(vector<int> v, vector< vector< int > > vlist);
	bool isPermutation(vector<int> v1, vector<int> v2);

	bool validTemp(vector< vector<int> > T, int card, vector<int> dirs);




public:

	Bundle(vector<lst> vars, vector< vector< double > > L, vector< double > offp, vector< double > offm, 	vector< vector< int > > T);

	int getDim(){ return this->dim; };
	int getSize(){ return L.size(); };
	int getCard(){ return T.size(); };
	int getNumDirs(){ return this->L.size(); };

	double getOffp(int i){ return this->offp[i]; };
	double getOffm(int i){ return this->offm[i]; };

	void setTemplate(vector< vector< int > > T);
	void setOffsetP(vector< double > offp){ this->offp = offp; }
	void setOffsetM(vector< double > offm){ this->offm = offm; }

	void canonize();
	vector< vector<int> > decompose();
	vector< vector<int> > decomposeRand();
	vector< vector<int> > decomposeTotalRand();
	pair< vector<double>, vector<double> > transform(lst vars, lst f, bool mode);

	LinearSystem *getBundle();
	Parallelotope* getParallelotope(int i);

	virtual ~Bundle();
};

#endif /* BUNDLE_H_ */
