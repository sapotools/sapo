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
  std::vector<unsigned int> degrees; // collection of degrees
  std::vector<unsigned int> shifts;  // degree shifts
  std::vector<GiNaC::ex> coeffs;     // list of the coefficients

  // TODO: complete two the following methods if necessary. Otherwise, remove
  // them operations on multidimensional matrices
  std::pair<std::vector<GiNaC::ex>, std::vector<std::vector<unsigned int>>>
  compressZeroCoeffs() const;
  void implicitMaxIndex() const;

  /**
   * Encode a multi-index into a coefficient position index
   *
   * @param[in] multi_index multi-index to convert
   * @returns converted multi-index
   */
  unsigned int
  multi_index2pos(const std::vector<unsigned int> &multi_index) const;

  /**
   * Decode a position into a multi-index
   *
   * @param[in] position position to convert
   * @returns converted position
   */
  std::vector<unsigned int> pos2multi_index(unsigned int position) const;

  /**
   * Initialize the degree shift vector used to extract the multi-indices
   */
  void initShifts();

  /**
   * Recursively extract the coefficients of the polynomial and populate the
   * vector of coefficients
   *
   * @param[in] polynomial polynomial from which to extract the coefficients
   */
  void initCoeffs(const GiNaC::ex &polynomial)
  {
    initCoeffs(polynomial, 0, 0);
  }

  /**
   * Recursively extract the coefficients of the polynomial and populate the
   * vector of coefficients
   *
   * @param[in] polynomial polynomial from which to extract the coefficients
   * @param[in] var_idx index of the variable to be considered
   * @param[in] pos is the position of the next coefficients to be initialized
   */
  void initCoeffs(const GiNaC::ex &polynomial, unsigned int var_idx,
                  const unsigned int pos); // extract and initialize
                                           // the coefficients
                                           // of the polynomial

public:
  // constructors
  BaseConverter(const GiNaC::lst &vars, const GiNaC::ex &polynomial);
  BaseConverter(const GiNaC::lst &vars, const GiNaC::ex &polynomial,
                const std::vector<unsigned int> &degrees);
  BaseConverter(const GiNaC::lst &vars, const GiNaC::ex &num,
                const GiNaC::ex &denom);

  // get Bernstein coefficients
  GiNaC::ex bernCoeff(const std::vector<unsigned int> &mi) const;
  GiNaC::lst getBernCoeffs() const;
  GiNaC::lst getRationalBernCoeffs() const;
  GiNaC::lst getBernCoeffsMatrix() const;

  // operations of split
  std::vector<std::vector<unsigned int>> getMultiIdxList() const;
  void split(unsigned int direction, double split_point) const;
  void print() const;

  virtual ~BaseConverter();
};

#endif /* BASECONVERTER_H_ */
