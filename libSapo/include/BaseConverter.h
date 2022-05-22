/**
 * @file BaseConverter.h
 * Convert the basis of the given polynomial
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#ifndef BASECONVERTER_H_
#define BASECONVERTER_H_

#include <vector>

#include "SymbolicAlgebra.h"

/**
 * @brief A base converter for polynomials
 *
 * This class converts polynomials from their canonical base
 * to Bernstein base.
 */
class BaseConverter
{
private:
  const std::vector<SymbolicAlgebra::Symbol<>> &vars; // list of variables
  SymbolicAlgebra::Expression<> polynomial; // polynomial in symbolic form
  SymbolicAlgebra::Expression<> num,
      denom; // numerator and denominator for rational polynomials
  std::vector<unsigned int> degrees; // collection of degrees
  std::vector<unsigned int> shifts;  // degree shifts
  std::vector<SymbolicAlgebra::Expression<>>
      coeffs; // list of the coefficients

  // TODO: complete two the following methods if necessary. Otherwise, remove
  // them operations on multidimensional matrices
  std::pair<std::vector<SymbolicAlgebra::Expression<>>,
            std::vector<std::vector<unsigned int>>>
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
  void initCoeffs(const SymbolicAlgebra::Expression<> &polynomial)
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
  void initCoeffs(const SymbolicAlgebra::Expression<> &polynomial,
                  unsigned int var_idx,
                  const unsigned int pos); // extract and initialize
                                           // the coefficients
                                           // of the polynomial

public:
  // constructors
  BaseConverter(const std::vector<SymbolicAlgebra::Symbol<>> &vars,
                const SymbolicAlgebra::Expression<> &polynomial);
  BaseConverter(const std::vector<SymbolicAlgebra::Symbol<>> &vars,
                const SymbolicAlgebra::Expression<> &polynomial,
                const std::vector<unsigned int> &degrees);
  BaseConverter(const std::vector<SymbolicAlgebra::Symbol<>> &vars,
                const SymbolicAlgebra::Expression<> &num,
                const SymbolicAlgebra::Expression<> &denom);

  // get Bernstein coefficients
  SymbolicAlgebra::Expression<>
  bernCoeff(const std::vector<unsigned int> &mi) const;
  std::vector<SymbolicAlgebra::Expression<>> getBernCoeffs() const;
  // std::vector<SymbolicAlgebra::Expression<>> getRationalBernCoeffs() const;
  std::vector<SymbolicAlgebra::Expression<>> getBernCoeffsMatrix() const;

  // operations of split
  std::vector<std::vector<unsigned int>> getMultiIdxList() const;
  void split(unsigned int direction, double split_point) const;
  void print() const;

  virtual ~BaseConverter();
};

#endif /* BASECONVERTER_H_ */
