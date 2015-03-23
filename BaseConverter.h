/*
 * BaseConverter.h
 *
 *  Created on: Oct 23, 2014
 *      Author: dreossi
 */

#ifndef BASECONVERTER_H_
#define BASECONVERTER_H_

#include "Common.h"

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
	vector<int> pos2multi_index(int position);								// convert a position to a multi-index
	void extractCoeffs(ex polynomial,int var_idx,vector<int> multi_index);	// extract the coefficients of the polynomial

	int nChoosek(int n, int k);								// binomial coefficient
	int multi_index_nChoosek(vector<int> n, vector<int> k);	// binomial coefficient of multi-indices
	bool multi_index_leq(vector<int> a, vector<int> b);		// check whether b dominates a

public:
	BaseConverter(lst vars, ex polynomial);
	BaseConverter(lst vars, ex polynomial, vector<int> degrees);
	BaseConverter(lst vars, ex num, ex denom);
	ex bernCoeff(vector<int> mi);
	lst getBernCoeffs();
	lst getRationalBernCoeffs();
	vector< vector< int > > getMultiIdxList();
	void split(int direction, double split_point);

	virtual ~BaseConverter();
};

#endif /* BASECONVERTER_H_ */
