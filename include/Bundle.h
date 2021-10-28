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
	using Vector = vector<double>;
	using Matrix = vector<Vector>;
	using CtrlPointType = map< vector<int>,pair<lst,lst> >;

private:
	unsigned int dim;					// dimension
	Matrix L;		// direction matrix
	Vector offp;				// superior offset
	Vector offm;				// inferior offset
	vector< vector< int > > T;			// templates matrix
	Matrix Theta;	// matrix of orthogonal proximity
	vector<lst> vars;					// variables appearing in generato function
										// vars[0] q: base vertex
										// vars[1] alpha : free variables \in [0,1]
										// vars[2] beta : generator amplitudes

	// map with Bernstein coefficients
	map< vector<int>, lst > bernCoeffs;

	double initTheta();
	Vector offsetDistances();

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
	
	bool isIn(int n, vector<int> v);
	bool isIn(vector<int> v, vector< vector< int > > vlist);
	bool isPermutation(vector<int> v1, vector<int> v2);


	bool validTemp(vector< vector<int> > T, unsigned int card, vector<int> dirs);	// check if a template is valid
	vector<lst> transformContrPts(lst vars, lst f, int mode);

public:

	// constructors
	Bundle(const Bundle &orig);
	Bundle(Bundle&& orig);
	Bundle(const Matrix& L, const Vector& offp, const Vector& offm,
		   const vector< vector< int > >& T);
	Bundle(const vector<lst>& vars, const Matrix& L, const Vector& offp,
		   const Vector& offm, const vector< vector< int > >& T);

	inline unsigned int getDim() const
	{
		return this->dim;
	}

	inline unsigned int getSize() const
	{
		return L.size();
	}

	inline unsigned int getCard() const
	{
		return T.size();
	}

	inline unsigned int getNumDirs() const
	{
		return this->L.size();
	}

	inline const vector<int>& getTemplate(long unsigned int i) const
	{
		return this->T[i];
	}

	inline const std::vector< std::vector<double> >& getDirectionMatrix() const
	{ 
		return this->L;
	}

	inline const double& getOffp(long unsigned int i) const
	{
		return this->offp[i];
	}

	inline const double& getOffm(long unsigned int i) const
	{
		return this->offm[i];
	}

	LinearSystem getLinearSystem() const;
	Parallelotope getParallelotope(unsigned int i) const;

	void setTemplate(vector< vector< int > > T);
	inline void setOffsetP(Vector offp)
	{
		this->offp = offp;
	}

	void setOffsetM(Vector offm)
	{
		this->offm = offm;
	}

	// operations on bundles
	Bundle get_canonical() const;
	Bundle decompose(double alpha, int max_iters);
	Bundle transform(const lst& vars, const lst& f, 
					 CtrlPointType &controlPts, int mode) const;
	Bundle transform(const lst& vars, const lst& params, 
					 const lst& f, const LinearSystem& paraSet,
					 CtrlPointType &controlPts, int mode) const;

	Bundle& operator=(Bundle&&);

	virtual ~Bundle();

	friend void swap(Bundle& A, Bundle& B);
};

void swap(Bundle& A, Bundle& B);

#endif /* BUNDLE_H_ */
