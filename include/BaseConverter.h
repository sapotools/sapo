/**
 * @file BaseConverter.h
 * Convert the basis of the given polynomial
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#ifndef BASECONVERTER_H_
#define BASECONVERTER_H_

#include "Common.h"
#include <math.h>

class BaseConverter {
private:
	lst vars;				// list of variables
	ex polynomial;			// polynomial in symbolic form
	ex num, denom;			// numerator and denominator for rational polynomials
	vector<int> degrees;	// collection of degrees
	vector<int> shifts;		// degree shifts
	vector<ex> coeffs;		// list of the coefficients

	void initShifts();														// initialize degrees shifts
	int multi_index2pos(vector<int> multi_index);							// convert a multi-index to a position
	vector<int> pos2multi_index(unsigned int position);								// convert a position to a multi-index
	void extractCoeffs(ex polynomial, unsigned int var_idx,vector<int> multi_index);	// extract the coefficients of the polynomial

	int nChoosek(int n, int k);								// binomial coefficient
	int multi_index_nChoosek(vector<int> n, vector<int> k);	// binomial coefficient of multi-indices
	bool multi_index_leq(vector<int> a, vector<int> b);		// check whether b dominates a

	// auxiliary operations
	int prod( vector<int> v, int a, int b );											// productory of element v[a]*...*v[b]
	int nchoosek(int n, int k);															// n choose k
	vector< vector<ex> > matrixProd(vector< vector<ex> > A, vector< vector<ex> > B);	// matrix product

	void print( vector< vector< ex > > M);

public:
	// constructors
	BaseConverter(lst vars, ex polynomial);
	BaseConverter(lst vars, ex polynomial, vector<int> degrees);
	BaseConverter(lst vars, ex num, ex denom);

	// get Bernstein coefficients
	ex bernCoeff(vector<int> mi);
	lst getBernCoeffs();
	lst getRationalBernCoeffs();
	lst getBernCoeffsMatrix();

	// operations on multi-indices
	vector< int > n2t( vector<int> a, vector<int> degs );
	vector< int > t2n( vector<int> a, vector<int> degs );

	//operations on mutlidimensional matrices
	vector< vector< ex > > genUtilde(int dim);
	vector< int > transp( vector<int> b, vector<int> degs, int degs_prod );
	vector< vector< ex > > transp( vector< vector<ex> > M, vector<int> degs );
	vector< int > transp_naive( vector<int> b, vector<int> degs );
	pair< vector<ex>, vector< vector<int> > > compressZeroCoeffs();
	void implicitMaxIndex();
	vector<int> shift(vector<int> v);

	// operations of split
	vector< vector< int > > getMultiIdxList();
	void split(long unsigned int direction, double split_point);
	void print();

	virtual ~BaseConverter();
};

#endif /* BASECONVERTER_H_ */
