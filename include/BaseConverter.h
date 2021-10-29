/**
 * @file BaseConverter.h
 * Convert the basis of the given polynomial
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#ifndef BASECONVERTER_H_
#define BASECONVERTER_H_

#include <math.h>

#include "Common.h"

class BaseConverter
{
private:
  const GiNaC::lst &vars; // list of variables
  GiNaC::ex polynomial;   // polynomial in symbolic form
  GiNaC::ex num, denom;   // numerator and denominator for rational polynomials
  std::vector<int> degrees;      // collection of degrees
  std::vector<int> shifts;       // degree shifts
  std::vector<GiNaC::ex> coeffs; // list of the coefficients

  void initShifts(); // initialize degrees shifts
  int multi_index2pos(
      std::vector<int> multi_index); // convert a multi-index to a position
  std::vector<int> pos2multi_index(
      unsigned int position); // convert a position to a multi-index
  void extractCoeffs(GiNaC::ex polynomial, unsigned int var_idx,
                     std::vector<int> multi_index); // extract the coefficients
                                                    // of the polynomial

  int nChoosek(int n, int k); // binomial coefficient
  int multi_index_nChoosek(
      std::vector<int> n,
      std::vector<int> k); // binomial coefficient of multi-indices
  bool multi_index_leq(std::vector<int> a,
                       std::vector<int> b); // check whether b dominates a

  // auxiliary operations
  int prod(std::vector<int> v, int a,
           int b);            // productory of element v[a]*...*v[b]
  int nchoosek(int n, int k); // n choose k
  std::vector<std::vector<GiNaC::ex>>
  matrixProd(std::vector<std::vector<GiNaC::ex>> A,
             std::vector<std::vector<GiNaC::ex>> B); // matrix product

  void print(std::vector<std::vector<GiNaC::ex>> M);

public:
  // constructors
  BaseConverter(const GiNaC::lst &vars, const GiNaC::ex &polynomial);
  BaseConverter(const GiNaC::lst &vars, const GiNaC::ex &polynomial,
                const std::vector<int> &degrees);
  BaseConverter(const GiNaC::lst &vars, const GiNaC::ex &num,
                const GiNaC::ex &denom);

  // get Bernstein coefficients
  GiNaC::ex bernCoeff(std::vector<int> mi);
  GiNaC::lst getBernCoeffs();
  GiNaC::lst getRationalBernCoeffs();
  GiNaC::lst getBernCoeffsMatrix();

  // operations on multi-indices
  std::vector<int> n2t(std::vector<int> a, std::vector<int> degs);
  std::vector<int> t2n(std::vector<int> a, std::vector<int> degs);

  // operations on mutlidimensional matrices
  std::vector<std::vector<GiNaC::ex>> genUtilde(int dim);
  std::vector<int> transp(std::vector<int> b, std::vector<int> degs,
                          int degs_prod);
  std::vector<std::vector<GiNaC::ex>>
  transp(std::vector<std::vector<GiNaC::ex>> M, std::vector<int> degs);
  std::vector<int> transp_naive(std::vector<int> b, std::vector<int> degs);
  std::pair<std::vector<GiNaC::ex>, std::vector<std::vector<int>>>
  compressZeroCoeffs();
  void implicitMaxIndex();
  std::vector<int> shift(std::vector<int> v);

  // operations of split
  std::vector<std::vector<int>> getMultiIdxList();
  void split(long unsigned int direction, double split_point);
  void print();

  virtual ~BaseConverter();
};

#endif /* BASECONVERTER_H_ */
