/**
 * @file LinearAlgebra.h
 * Contains linear algebra classes and functions code
 *
 * @author Alberto Casagrande <acasagrande@units.it>
 * @version 0.1
 */

#ifndef LINEAR_ALGEBRA_H
#define LINEAR_ALGEBRA_H

#include <iostream>
#include <vector>
#include <map>
#include <set>

#include <numeric> // due to std::iota
#include <math.h>  // due to sqrt
#include <cmath>   // due to floor
#include <limits>  // due to numeric_limits

/**
 * @brief Approximate all the values in a vector up to given decimal digit.
 *
 * @tparam T is a numerical type
 * @param v is a the vector
 * @param decimal is the number of decimal number in the result.
 * @return the given vector approximated up to `decimal` decimal digits.
 */
template<typename T>
inline std::vector<T> approx(const std::vector<T> &v,
                             const unsigned short decimal
                             = std::numeric_limits<T>::digits10);

/**
 * @brief Approximate all the values in a vector up to given decimal digit.
 *
 * This function can be mainly used to speed up the computation.
 *
 * @tparam T is a numerical type
 * @param v is a the vector
 * @param decimal is the number of decimal number in the result.
 * @return the given vector approximated up to `decimal` decimal digits.
 */
template<typename T>
inline std::vector<T> approx(std::vector<T> &&orig,
                             const unsigned short decimal
                             = std::numeric_limits<T>::digits10);

/**
 * @brief Approximate all the values in a vector up to given decimal digit.
 *
 * This function can be mainly used to speed up the computation.
 *
 * @param v is a the vector
 * @param decimal is the number of decimal number in the result.
 * @return the given vector approximated up to `decimal` decimal digits.
 */
template<>
inline std::vector<double> approx<double>(std::vector<double> &&orig,
                                          const unsigned short decimal)
{
  if (decimal == std::numeric_limits<double>::digits10) {
    return orig;
  }

  const double mult = std::pow(10, decimal);
  for (auto it = std::begin(orig); it != std::end(orig); ++it) {
    *it = floor(*it * mult) / mult;
  }

  return orig;
}

/**
 * @brief Approximate all the values in a vector up to given decimal digit.
 *
 * This function can be mainly used to speed up the computation.
 *
 * @param v is a the vector
 * @param decimal is the number of decimal number in the result.
 * @return the given vector approximated up to `decimal` decimal digits.
 */
template<>
inline std::vector<double> approx<double>(const std::vector<double> &orig,
                                          const unsigned short decimal)
{
  return approx<double>(std::vector<double>(orig), decimal);
}

/**
 * @brief Compute the Euclidean norm of a vector
 *
 * @tparam T is a numeric type supporting the `sqrt` function.
 * @param v is a vector
 * @return the Euclidean norm of the given vector
 */
template<typename T>
T norm_2(const std::vector<T> &v)
{
  T norm = 0;

  for (auto v_it = std::begin(v); v_it != std::end(v); ++v_it) {
    norm += *v_it * *v_it;
  }

  return sqrt(norm);
}

/**
 * @brief Compute the Manhattan norm of a vector
 *
 * @tparam T is a numeric type supporting the `abs` function.
 * @param v is a vector
 * @return the Manhattan norm of the given vector
 */
template<typename T>
T norm_1(const std::vector<T> &v)
{
  T norm = 0;

  for (auto v_it = std::begin(v); v_it != std::end(v); ++v_it) {
    norm += std::abs(*v_it);
  }

  return norm;
}

/**
 * @brief Compute the maximum norm of a vector
 *
 * @tparam T is a numeric type supporting the `abs` function.
 * @param v is a vector
 * @return the maximum norm of the given vector
 */
template<typename T>
T norm_infinity(const std::vector<T> &v)
{
  T norm = 0;

  for (auto v_it = std::begin(v); v_it != std::end(v); ++v_it) {
    if (std::abs(*v_it) > norm) {
      norm = std::abs(*v_it);
    }
  }

  return norm;
}

/**
 * Compute the complementary of a vector of values.
 *
 * @param[in] orig the vector of values to be complementated.
 * @return A vector containing the complementaries of the parameter values
 */
template<typename T>
std::vector<T> operator-(const std::vector<T> &orig)
{
  std::vector<T> res = orig;

  transform(res.begin(), res.end(), res.begin(), std::negate<T>());

  return res;
}

/**
 * Compute the complementary of a vector of values.
 *
 * @param[in] orig the vector of values to be complementated.
 * @return A vector containing the complementaries of the parameter values
 */
template<typename T>
std::vector<T> operator-(std::vector<T> &&v)
{
  transform(v.begin(), v.end(), v.begin(), std::negate<T>());

  return v;
}

/**
 * @brief Compute the element-wise sum of two vectors.
 *
 * @tparam T any numeric type
 * @param a is the first vector to be added.
 * @param b is the second vector to be added.
 * @return the element-wise sum of the two parameters.
 */
template<typename T>
std::vector<T> operator+(const std::vector<T> &a, const std::vector<T> &b)
{
  if (a.size() != b.size()) {
    throw std::domain_error("The two vectors have different dimensions");
  }

  std::vector<T> res(a.size());

  for (unsigned int i = 0; i < res.size(); ++i) {
    res[i] = a[i] + b[i];
  }

  return res;
}

/**
 * @brief Compute the element-wise sum of two vectors.
 *
 * @tparam T any numeric type
 * @param a is the first vector to be added.
 * @param b is the second vector to be added.
 * @return the element-wise sum of the two parameters.
 */
template<typename T>
std::vector<T> operator+(std::vector<T> &&a, const std::vector<T> &b)
{
  if (a.size() != b.size()) {
    throw std::domain_error("The two vectors have different dimensions");
  }

  for (unsigned int i = 0; i < a.size(); ++i) {
    a[i] += b[i];
  }

  return a;
}

/**
 * @brief Compute the element-wise sum of two vectors.
 *
 * @tparam T any numeric type
 * @param a is the first vector to be added.
 * @param b is the second vector to be added.
 * @return the element-wise sum of the two parameters.
 */
template<typename T>
inline std::vector<T> operator+(const std::vector<T> &a, std::vector<T> &&b)
{
  return operator+(b, a);
}

/**
 * @brief Compute the element-wise difference of two vectors.
 *
 * @tparam T any numeric type
 * @param a is the vector from which the second parameter must be subtracted.
 * @param b is the vector that must be subtracted.
 * @return the element-wise difference of the second parameter from the first
 * one.
 */
template<typename T>
std::vector<T> operator-(const std::vector<T> &a, const std::vector<T> &b)
{
  if (a.size() != b.size()) {
    throw std::domain_error("The two vectors have different dimensions");
  }

  std::vector<T> res(a.size());

  for (unsigned int i = 0; i < res.size(); ++i) {
    res[i] = a[i] - b[i];
  }

  return res;
}

/**
 * @brief Compute the element-wise difference of two vectors.
 *
 * @tparam T any numeric type
 * @param a is the vector from which the second parameter must be subtracted.
 * @param b is the vector that must be subtracted.
 * @return the element-wise difference of the second parameter from the first
 * one.
 */
template<typename T>
std::vector<T> operator-(std::vector<T> &&a, const std::vector<T> &b)
{
  if (a.size() != b.size()) {
    throw std::domain_error("The two vectors have different dimensions");
  }

  for (unsigned int i = 0; i < a.size(); ++i) {
    a[i] -= b[i];
  }

  return a;
}

/**
 * @brief Compute the element-wise difference of two vectors.
 *
 * @tparam T any numeric type
 * @param a is the vector from which the second parameter must be subtracted.
 * @param b is the vector that must be subtracted.
 * @return the element-wise difference of the second parameter from the first
 * one.
 */
template<typename T>
std::vector<T> operator-(const std::vector<T> &a, std::vector<T> &&b)
{
  if (a.size() != b.size()) {
    throw std::domain_error("The two vectors have different dimensions");
  }

  for (unsigned int i = 0; i < a.size(); ++i) {
    b[i] -= a[i];
    b[i] *= -1;
  }

  return a;
}

/**
 * @brief Compute the element-wise scalar product.
 *
 * @tparam T any numeric type.
 * @param s is a scalar value.
 * @param v is a vector.
 * @return the element-wise scalar product $s * v$.
 */
template<typename T>
std::vector<T> operator*(const T s, const std::vector<T> &v)
{
  std::vector<T> res(v.size());

  for (unsigned int i = 0; i < res.size(); ++i) {
    res[i] = s * v[i];
  }

  return res;
}

/**
 * @brief Compute the element-wise scalar product.
 *
 * @tparam T any numeric type.
 * @param v is a vector.
 * @param s is a scalar value.
 * @return the element-wise scalar product $v * s$.
 */
template<typename T>
inline std::vector<T> operator*(const std::vector<T> &v, const T s)
{
  return operator*(s, v);
}

/**
 * @brief Compute the vector product.
 *
 * @tparam T any numeric type.
 * @param v1 is a vector.
 * @param v2 is a vector.
 * @return the vector product $v1 \cdot v2$.
 */
template<typename T>
T operator*(const std::vector<T> &v1, const std::vector<T> &v2)
{
  if (v1.size() != v2.size()) {
    throw std::domain_error("The two vectors have not the same dimensions.");
  }

  T res = 0;
  auto it2 = std::begin(v2);
  for (auto it1 = std::begin(v1); it1 != std::end(v1); ++it1) {
    res += (*it1) * (*it2);

    ++it2;
  }

  return res;
}

/**
 * @brief Compute the Hadamard product of two vectors.
 *
 * This method computes the Hadamard product of two vectors
 * $v_1 \circ v_2$, i.e., it returns a vector whose elements
 * are obtained by element-wise multiplying of the
 * parameters.
 * $$(v_1 \circ v_2)[i] = v_1[i]*v_2[i]$$
 *
 * @tparam T any numeric type.
 * @param v1 is a vector.
 * @param v2 is a vector.
 * @return the vector product $v1 \circ v2$.
 */
template<typename T>
std::vector<T> H_prod(const std::vector<T> &v1, const std::vector<T> &v2)
{
  if (v1.size() != v2.size()) {
    throw std::domain_error("The two vectors have not the same dimensions.");
  }

  std::vector<T> res(v1.size());
  for (unsigned int i = 0; i < v1.size(); ++i) {
    res[i] = v1[i] * v2[i];
  }

  return res;
}

/**
 * @brief Compute the element-wise scalar product.
 *
 * @tparam T any numeric type.
 * @param s is a scalar value.
 * @param v is a vector.
 * @return the element-wise scalar product $s * v$.
 */
template<typename T>
std::vector<T> operator*(const T s, std::vector<T> &&v)
{
  for (auto v_it = std::begin(v); v_it != std::end(v); ++v_it) {
    *v_it *= s;
  }

  return v;
}

/**
 * @brief Compute the element-wise scalar product.
 *
 * @tparam T any numeric type.
 * @param v is a vector.
 * @param s is a scalar value.
 * @return the element-wise scalar product $v * s$.
 */
template<typename T>
inline std::vector<T> operator*(std::vector<T> &&v, const T s)
{
  return operator*(s, v);
}

/**
 * @brief Compute the element-wise scalar division.
 *
 * @tparam T any numeric type.
 * @param v is a vector.
 * @param s is a scalar value.
 * @return the element-wise scalar division $v/a$.
 */
template<typename T>
std::vector<T> operator/(const std::vector<T> &v, const T s)
{
  std::vector<T> res(v.size());

  for (unsigned int i = 0; i < res.size(); ++i) {
    res[i] = v[i] / s;
  }

  return res;
}

/**
 * @brief Compute the element-wise scalar division.
 *
 * @tparam T any numeric type.
 * @param v is a vector.
 * @param s is a scalar value.
 * @return the element-wise scalar division $v/a$.
 */
template<typename T>
std::vector<T> operator/(std::vector<T> &&v, const T s)
{
  for (auto v_it = std::begin(v); v_it != std::end(v); ++v_it) {
    *v_it /= s;
  }

  return v;
}

/**
 * Compute the angle between two vectors
 *
 * @tparam T is a numeric type
 * @param[in] v1 vector
 * @param[in] v2 vector
 * @returns angle between v1 and v2
 */
template<typename T>
inline T angle(const std::vector<T> &v1, const std::vector<T> &v2)
{
  return acos(v1 * v2 / (norm_2(v1) * norm_2(v2)));
}

/**
 * @brief Compute the row-column matrix-matrix multiplication
 *
 * @tparam T is a numeric type
 * @param A is a dense matrix
 * @param B is a dense matrix
 * @return the row-column matrix-matrix multiplication $A \cdot B$.
 */
template<typename T>
std::vector<std::vector<T>> operator*(const std::vector<std::vector<T>> &A,
                                      const std::vector<std::vector<T>> &B)
{
  if (((A.size() != 0 && A.front().size() == 0) || (A.size() == 0))
      && B.size() == 0) {
    return std::vector<std::vector<T>>();
  }

  if (A.front().size() != B.size()) {
    throw std::domain_error("The two matrices are not compatible for the "
                            "matrix-matrix multiplication.");
  }

  std::vector<std::vector<T>> res(A.size(),
                                  std::vector<T>(B.front().size(), 0));

  for (unsigned int row_idx = 0; row_idx < A.size(); ++row_idx) {
    for (unsigned int col_idx = 0; col_idx < B.front().size(); ++col_idx) {
      T value = 0;
      for (unsigned int k = 0; k < A[row_idx].size(); ++k) {
        value += A[row_idx][k] * B[k][col_idx];
      }
      res[row_idx][col_idx] = value;
    }
  }

  return res;
}

/**
 * @brief A class to represent pemutations.
 *
 * Permutation are square matrices whose row and column
 * sums are one. However, storing a permutation for
 * $n$-dim vectors by using this natural representation
 * requires $n^2$ integer numbers and its application
 * costs $O(n^2)$.
 * This class stores permutation is a more efficient
 * way representing only swapped element positions in
 * a map from integer to integer. The application cost
 * belongs to $O(n)$.
 */
class Permutation
{
  std::map<unsigned int, unsigned int> _swaps; //!< the swapped positions
public:
  /**
   * @brief Create a new Permutation object
   */
  Permutation(): _swaps() {}

  /**
   * @brief Copy constructor for Permutation object
   *
   * @param orig
   */
  Permutation(const Permutation &orig): _swaps(orig._swaps) {}

  /**
   * @brief Swaps two positions.
   *
   * This method swaps two positions in the permutation.
   *
   * @param i first position to be swapped.
   * @param j second position to be swapped.
   * @return The updated permutation.
   */
  Permutation &swap(const unsigned int i, const unsigned int j)
  {
    // for both the positions
    for (const unsigned int &pos: {i, j}) {
      // if the position is not stored among the swaps yet
      if (_swaps.find(pos) == std::end(_swaps)) {

        // store it
        _swaps[pos] = pos;
      }
    }

    // swap the two positions
    std::swap(_swaps[i], _swaps[j]);

    return *this;
  }

  template<typename T>
  std::vector<T> operator()(const std::vector<T> &v) const
  {
    std::vector<T> pv = v;

    for (auto it = std::begin(_swaps); it != std::end(_swaps); ++it) {
      if (it->first >= v.size() || it->second >= v.size()) {
        throw std::domain_error(
            "Permutation and vector dimensions are not compatible.");
      }
      pv[it->first] = v[it->second];
    }

    return pv;
  }

  Permutation &operator=(const Permutation &P)
  {
    _swaps = P._swaps;

    return *this;
  }

  friend inline void swap(Permutation &P1, Permutation &P2)
  {
    std::swap(P1._swaps, P2._swaps);
  }

  friend std::ostream &operator<<(std::ostream &os, const Permutation &P)
  {
    os << "Permutation[";
    for (auto it = std::begin(P._swaps); it != std::end(P._swaps); ++it) {
      if (it != std::begin(P._swaps)) {
        os << ",";
      }
      os << it->second << "->" << it->first;
    }
    os << "]";

    return os;
  }
};

/*!
 *  \addtogroup DenseLinearAlgebra
 *  @{
 */

//! Linear algebra for dense matrices
namespace DenseLinearAlgebra
{

/**
 * @brief A Permutation-LU factorization for dense matrices.
 *
 * This class represents a Permutation, Lower-triangular,
 * Upper-triangular factorization for dense matrices. It
 * also offers a method to solve linear systems
 * represented by dense matrices. The permutation is not
 * represented by a matrix, but instead by using a
 * vector of indexes swaps.
 *
 * @tparam T any numeric type.
 */
template<typename T>
class PLU_Factorization
{
  Permutation _P;                  //!< the permutation
  std::vector<std::vector<T>> _LU; //!< the L and U matrices

  /**
   * @brief Solve the linear system $L\cdot x = b$.
   *
   * @param[in] b is the known term of the linear system.
   * @return the vector $x = L^{-1} \cdot b$.
   */
  std::vector<T> solveL(const std::vector<T> &b) const
  {
    std::vector<T> x(b.size());

    for (unsigned int j = 0; j < _LU.size(); ++j) {
      T value = b[j];
      const std::vector<T> &row = _LU[j];
      for (unsigned int i = 0; i < j; ++i) {
        if (row[i] != 0) {
          value -= row[i] * x[i];
        }
      }
      x[j] = value;
    }

    return x;
  }

  /**
   * @brief Solve the linear system $U\cdot x = b$.
   *
   * @param[in] b is the known term of the linear system.
   * @return the vector $x = U^{-1} \cdot b$.
   */
  std::vector<T> solveU(const std::vector<T> &b) const
  {
    std::vector<T> x(b.size());

    for (unsigned int rj = 0; rj < _LU.size(); ++rj) {
      const unsigned int j = _LU.size() - rj - 1;
      T value = b[j];
      const std::vector<T> &row = _LU[j];

      const T &diag_value = row[j];
      if (diag_value == 0) {
        throw std::domain_error("The linear system is underdetermined.");
      }
      for (unsigned int i = j + 1; i < row.size(); ++i) {
        value -= row[i] * x[i];
      }
      x[j] = value / diag_value;
    }

    return x;
  }

public:
  /**
   * @brief Create a new PLU_Factorization object
   *
   */
  PLU_Factorization(): _P(), _LU() {}

  /**
   * @brief Create a new PLU_Factorization object
   *
   * @param orig is the model for the new object.
   */
  PLU_Factorization(const PLU_Factorization &orig): _P(orig._P), _LU(orig._LU)
  {
  }

  /**
   * @brief Create a new PLU_Factorization object
   *
   * This constructor performs the P-LU factorization of
   * a dense matrix.
   *
   * @param M is the matrix whose P-LU factorization must be computed.
   */
  PLU_Factorization(const std::vector<std::vector<T>> &M): _P(), _LU(M)
  {
    if (M.size() == 0) {
      throw std::domain_error("The parameter is an empty matrix.");
    }

    for (unsigned int j = 0; j < _LU.size(); ++j) {
      if (j == M.front().size()) {
        return;
      }

      // find the first non-null value below A[j-1][j]
      unsigned int k = j;
      while (k < _LU.size() && _LU[k][j] == 0) {
        k++;
      }

      // if it does not exist, skip to the next row/column
      if (k < _LU.size()) {
        // otherwise, swap rows
        if (j != k) {
          std::swap(_LU[j], _LU[k]);
          _P.swap(j, k);
        }

        // nullify all the values below A[j][j]
        const std::vector<T> &row_j = _LU[j];
        const T &v = row_j[j];
        k += 1;
        while (k < _LU.size()) {
          std::vector<T> &row_k = _LU[k];
          if (row_k[j] != 0) {
            const T ratio = -row_k[j] / v;
            for (unsigned int i = j + 1; i < row_j.size(); ++i) {
              row_k[i] += ratio * row_j[i];
            }
            row_k[j] = -ratio;
          }
          ++k;
        }
      }
    }
  }

  /**
   * @brief Compute the solution of the system $(L \cot U) \cdot x = P^{-1}
   * \cdot b$
   *
   * @param is the known term of the linear system.
   * @return the vector $x = (U^{-1} \cdot L^{-1} \cdot P^{-1}) \cdot b$.
   */
  std::vector<T> solve(const std::vector<T> &b) const
  {
    if (b.size() != _LU.size()) {
      throw std::domain_error("Wrong dimensions");
    }

    if (_LU.size() != _LU.front().size()) {
      throw std::domain_error("The factorization is not square and the linear "
                              "system cannot be solved.");
    }

    std::vector<T> Pb = _P(b);

    // Forward substitution for L
    Pb = solveL(Pb);

    // Backward substitution for U
    return solveU(Pb);
  }

  /**
   * @brief Copy a PLU_Factorization object.
   *
   * @param orig is the model that must be copied.
   * @return a reference to the updated object.
   */
  PLU_Factorization<T> &operator=(const PLU_Factorization<T> &orig)
  {
    _P = orig._P;
    _LU = orig._LU;

    return *this;
  }

  /**
   * @brief Copy a PLU_Factorization object.
   *
   * @param orig is the model that must be copied.
   * @return a reference to the updated object.
   */
  PLU_Factorization<T> &operator=(PLU_Factorization<T> &&orig)
  {
    std::swap(_P, orig._P);
    std::swap(_LU, orig._LU);

    return *this;
  }

  /**
   * @brief Get the factorization permutation
   *
   * @return the factorization permutation
   */
  const Permutation &P() const
  {
    return _P;
  }

  /**
   * @brief Return the factorization lower-triangular matrix
   *
   * @return the factorization lower-triangular matrix.
   */
  std::vector<std::vector<T>> L() const
  {
    std::vector<std::vector<T>> res(_LU.size(), std::vector<T>(_LU.size(), 0));

    for (unsigned int row_idx = 0; row_idx < _LU.size(); ++row_idx) {
      const std::vector<T> &row = _LU[row_idx];
      std::vector<T> &res_row = res[row_idx];
      for (unsigned int col_idx = 0; col_idx < row_idx && col_idx < row.size();
           ++col_idx) {
        res_row[col_idx] = row[col_idx];
      }
      res_row[row_idx] = 1;
    }

    return res;
  }

  /**
   * @brief Return the factorization upper-triangular matrix
   *
   * @return the factorization upper-triangular matrix.
   */
  std::vector<std::vector<T>> U() const
  {
    std::vector<std::vector<T>> res(_LU.size(),
                                    std::vector<T>(_LU.front().size(), 0));

    for (unsigned int row_idx = 0; row_idx < _LU.size(); ++row_idx) {
      const std::vector<T> &row = _LU[row_idx];
      std::vector<T> &res_row = res[row_idx];
      for (unsigned int col_idx = row_idx; col_idx < row.size(); ++col_idx) {
        res_row[col_idx] = row[col_idx];
      }
    }

    return res;
  }

  template<typename E>
  friend std::ostream &std::operator<<(std::ostream &os,
                                       const PLU_Factorization<E> &D);
};

}
/*! @} End of DenseLinearAlgebra group */

/*!
 *  \addtogroup SparseLinearAlgebra
 *  @{
 */

//! Linear algebra for sparse matrices
namespace SparseLinearAlgebra
{

/**
 * @brief A class to represent sparse matrices.
 *
 * @tparam T any numeric type.
 */
template<typename T>
class Matrix
{
public:
  typedef std::map<unsigned int, T> RowType;

protected:
  class _row_ref_type;

  /**
   * @brief  A reference class for sparse matrix elements
   */
  class _elem_ref_type
  {
    Matrix<T> &A;               //!< a reference to the matrix
    const unsigned int row_idx; //!< the index of the element row
    const unsigned int col_idx; //!< the index of the element column

    /**
     * @brief Create a new matrix element reference object
     *
     * @param row is a reference to the element row
     * @param col_idx is the index of the element column
     */
    _elem_ref_type(_row_ref_type &row, const unsigned int col_idx):
        A(row.A), row_idx(row.row_idx), col_idx(col_idx)
    {
      if (col_idx >= A.num_of_rows()) {
        throw std::out_of_range("Out-of_bound: the column does not exist.");
      }
    }

  public:
    /**
     * @brief The assignment method
     *
     * This method assign a value to the referenced matrix element. If the
     * element already exists, it reassigns it with a new value. If the
     * matrix element does not exists yet, it creates the referenced
     * element (and its row if necessary) and assigns it with `value`.
     *
     * @param value is value that must be assign to the referenced element.
     * @return a reference to the assigned matrix element.
     */
    T &operator=(const T &value)
    {
      // TODO: delete the row if it only contains zeros
      auto row_it = A._matrix.find(row_idx);

      if (row_it == std::end(A._matrix)) {
        A._matrix[row_idx] = Matrix<T>::RowType();

        row_it = A._matrix.find(row_idx);
      }

      row_it->second[col_idx] = value;

      return row_it->second[col_idx];
    }

    /**
     * @brief The typecast operator for the type T
     *
     * @return The value of the matrix element.
     */
    operator T() const
    {
      auto row_it = A._matrix.find(row_idx);

      if (row_it == std::end(A._matrix)) {
        return 0;
      }

      auto elem_it = row_it->second.find(col_idx);

      if (elem_it == std::end(row_it->second)) {
        return 0;
      }

      return elem_it->second;
    }

    friend class _row_ref_type;
  };

  /**
   * @brief A reference class for sparse matrix rows
   */
  class _row_ref_type
  {
    Matrix<T> &A;               //!< a reference to the matrix
    const unsigned int row_idx; //!< the index of the referenced row

    /**
     * @brief Create a new row reference object
     *
     * @param A
     * @param row_idx
     */
    _row_ref_type(Matrix<T> &A, const unsigned int row_idx):
        A(A), row_idx(row_idx)
    {
      if (row_idx >= A.num_of_rows()) {
        throw std::out_of_range("Out-of_bound: the row does not exist.");
      }
    }

  public:
    /**
     * @brief Copy constructor for row reference objects
     *
     * @param orig is the model for the new object.
     */
    _row_ref_type(const _row_ref_type &orig): A(orig.A), row_idx(orig.row_idx)
    {
    }

    /**
     * @brief Return a reference of the element in the `col_idx`-th position of
     * the row
     *
     * @param col_idx is the index of the aimed value.
     * @return a reference of the element in position `col_idx`.
     */
    _elem_ref_type operator[](const unsigned int col_idx)
    {
      return _elem_ref_type(*this, col_idx);
    }

    friend class _elem_ref_type;
    friend class Matrix<T>;
  };

  /**
   * @brief A constant reference class for sparse matrix rows
   */
  class _const_row_ref_type
  {
    const Matrix<T> &A;         //!< a reference to the matrix
    const unsigned int row_idx; //!< the index of the referenced row

    /**
     * @brief Create a new constant row reference object
     *
     * @param A
     * @param row_idx
     */
    _const_row_ref_type(const Matrix<T> &A, const unsigned int row_idx):
        A(A), row_idx(row_idx)
    {
      if (row_idx >= A.num_of_rows()) {
        throw std::out_of_range("Out-of_bound: the row does not exist.");
      }
    }

  public:
    /**
     * @brief Copy constructor for constant row reference objects
     *
     * @param orig is the model for the new object.
     */
    _const_row_ref_type(const _const_row_ref_type &orig):
        A(orig.A), row_idx(orig.row_idx)
    {
    }

    /**
     * @brief Return the value in the `col_idx`-th position of the row
     *
     * @param col_idx is the index of the aimed value.
     * @return a copy of the value in position `col_idx`.
     */
    T operator[](const unsigned int col_idx) const
    {
      if (col_idx >= A.num_of_cols()) {
        throw std::out_of_range("Out-of_bound: the column does not exist.");
      }

      auto row_it = A._matrix.find(row_idx);

      if (row_it == std::end(A._matrix)) {
        return 0;
      }

      auto elem_it = row_it->second.find(col_idx);

      if (elem_it == std::end(row_it->second)) {
        return 0;
      }

      return elem_it->second;
    }

    friend class Matrix<T>;
  };

  std::map<unsigned int, RowType> _matrix; //!< indexed list of rows
  size_t _num_of_rows;                     //!< number of rows
  size_t _num_of_cols;                     //!< number of columns

  /**
   * @brief Transform a vector in a RowType
   *
   * @param row to original vector
   * @return the RowType representing `row`.
   */
  static RowType get_RowType(const std::vector<T> &row)
  {
    RowType row_map;
    unsigned int elem_idx = 0;
    for (auto elem_it = std::begin(row); elem_it != std::end(row); ++elem_it) {
      if (*elem_it != 0) {
        row_map[elem_idx] = *elem_it;
      }
      ++elem_idx;
    }

    return row_map;
  }

public:
  /**
   * @brief Create a new Matrix object
   *
   */
  Matrix(): _matrix(), _num_of_rows(0), _num_of_cols(0) {}

  /**
   * @brief Create a new Matrix object
   *
   * This method creates a new sparse matrix having a given number
   * of rows and columns whose elements equal zero.
   *
   * @param num_of_rows is the number of rows in the new matrix.
   * @param num_of_cols is the number of columns in the new matrix.
   */
  Matrix(const unsigned int num_of_rows, const unsigned int num_of_cols):
      _matrix(), _num_of_rows(num_of_rows), _num_of_cols(num_of_cols)
  {
  }

  /**
   * @brief Create a new Matrix object
   *
   * @param orig is the model for the new matrix.
   * @param up_to_row is the number of rows to be considered in `orig`.
   */
  Matrix(const std::vector<std::vector<T>> &orig,
         const unsigned int up_to_row):
      _matrix(),
      _num_of_rows(0), _num_of_cols(0)
  {
    if (up_to_row == 0) {
      return;
    }

    _num_of_cols = orig.front().size();

    for (unsigned int j = 0; j < up_to_row; ++j) {
      RowType s_row;
      const std::vector<T> &row = orig[j];
      if (row.size() != _num_of_cols) {
        std::domain_error("The parameter is not a matrix");
      }
      for (unsigned int i = 0; i < row.size(); ++i) {
        if (row[i] != 0) {
          s_row[i] = row[i];
        }
      }
      add_row(j, s_row);
    }
  }

  /**
   * @brief Create a new Matrix object
   *
   * @param orig is the model for the new matrix.
   */
  Matrix(const std::vector<std::vector<T>> &orig): Matrix(orig, orig.size()) {}

  /**
   * @brief Clone Matrix object
   *
   * @param orig is the model for the new object.
   */
  Matrix(const Matrix<T> &orig):
      _matrix(orig._matrix), _num_of_rows(orig._num_of_rows),
      _num_of_cols(orig._num_of_cols)
  {
  }

  /**
   * @brief Add a row to the matrix.
   *
   * @param row_idx is the index of the new row.
   * @param row is the new row.
   * @return a reference to the updated object.
   */
  Matrix<T> &add_row(unsigned row_idx, const RowType &row)
  {
    if (_matrix.count(row_idx) > 0) {
      std::domain_error("The row already exists.");
    }

    if (row_idx + 1 > _num_of_rows) {
      _num_of_rows = row_idx + 1;
    }

    auto greatest = std::rbegin(row);

    if (greatest != std::rend(row)) {
      if (_num_of_cols < greatest->first + 1) {
        _num_of_cols = greatest->first + 1;
      }

      _matrix[row_idx] = row;
    }

    return *this;
  }

  /**
   * @brief Add a row to the matrix.
   *
   * @param row_idx is the index of the new row.
   * @param row is the new row.
   * @return a reference to the updated object.
   */
  inline Matrix<T> &add_row(unsigned row_idx, const std::vector<T> &row)
  {
    return add_row(row_idx, get_RowType(row));
  }

  /**
   * @brief Add a row as last row of the matrix.
   *
   * @param row is the new row.
   * @return a reference to the updated object.
   */
  inline Matrix<T> &add_row(const RowType &row)
  {
    return add_row(this->num_of_rows(), row);
  }

  /**
   * @brief Add a row as last row of the matrix.
   *
   * @param row is the new row.
   * @return a reference to the updated object.
   */
  inline Matrix<T> &add_row(const std::vector<T> &row)
  {
    return add_row(this->num_of_rows(), get_RowType(row));
  }

  /**
   * @brief Get a constant reference to one of the matrix rows
   *
   * @param row_idx is the index of the aimed row.
   * @return the constant reference of the `row_idx`-th row.
   */
  typename Matrix<T>::_const_row_ref_type
  operator[](const unsigned int row_idx) const
  {
    return _const_row_ref_type(*this, row_idx);
  }

  /**
   * @brief Get a reference to one of the matrix rows
   *
   * @param row_idx is the index of the aimed row.
   * @return the reference of the `row_idx`-th row.
   */
  typename Matrix<T>::_row_ref_type operator[](const unsigned int row_idx)
  {
    return _row_ref_type(*this, row_idx);
  }

  /**
   * @brief Return the number of columns in the matrix
   *
   * @return the number of columns in the matrix.
   */
  size_t num_of_cols() const
  {
    return _num_of_cols;
  }

  /**
   * @brief Return the number of rows in the matrix
   *
   * @return the number of rows in the matrix.
   */
  size_t num_of_rows() const
  {
    return _num_of_rows;
  }

  /**
   * @brief Matrix-vector product
   *
   * @param v is the vector that must the multiplied to the matrix.
   * @return the product between the current matrix and `v`.
   */
  std::vector<T> operator*(const std::vector<T> &v) const
  {
    if (v.size() != _num_of_cols) {
      throw std::domain_error(
          "The matrix and the vector have non-compatible dimensions.");
    }

    std::vector<T> res(_num_of_rows);

    for (auto row_it = std::cbegin(_matrix); row_it != std::cend(_matrix);
         ++row_it) {
      T value = 0;
      for (auto elem_it = std::cbegin(row_it->second);
           elem_it != std::cend(row_it->second); ++elem_it) {
        value += elem_it->second * v[elem_it->first];
      }
      res[row_it->first] = value;
    }

    return res;
  }

  /**
   * @brief Copy a matrix
   *
   * @param orig is the matrix that must be copied.
   * @return a reference to the updated object.
   */
  Matrix<T> &operator=(const Matrix<T> &orig)
  {
    _matrix = orig._matrix;
    _num_of_cols = orig._num_of_cols;
    _num_of_rows = orig._num_of_rows;

    return *this;
  }

  /**
   * @brief Copy a matrix
   *
   * @param orig is the matrix that must be copied.
   * @return a reference to the updated object.
   */
  Matrix<T> &operator=(Matrix<T> &&orig)
  {
    std::swap(_matrix, orig._matrix);
    std::swap(_num_of_cols, orig._num_of_cols);
    std::swap(_num_of_rows, orig._num_of_rows);

    return *this;
  }

  /**
   * @brief Compute the row-column matrix-matrix multiplication
   *
   * @param A is a matrix
   * @return The row-column multiplication between this object and `A`.
   */
  Matrix<T> operator*(const Matrix<T> &A) const
  {
    if (num_of_cols() != A.num_of_rows()) {
      std::domain_error("The two matrices are not compatible for the "
                        "matrix-matrix multiplication.");
    }

    Matrix<T> res(num_of_rows(), A.num_of_cols());

    for (auto row_it = std::cbegin(_matrix); row_it != std::cend(_matrix);
         ++row_it) {
      const unsigned int row_idx = row_it->first;
      for (unsigned int col_idx = 0; col_idx < A.num_of_cols(); ++col_idx) {
        T value = 0;
        for (auto elem_it = std::cbegin(row_it->second);
             elem_it != std::cend(row_it->second); ++elem_it) {
          value += elem_it->second * A[elem_it->first][col_idx];
        }
        if (value != 0) {
          res._matrix[row_idx][col_idx] = value;
        }
      }
    }

    return res;
  }

  template<typename E>
  friend class PLU_Factorization;

  template<typename E>
  friend std::ostream &std::operator<<(std::ostream &os, const Matrix<E> &M);
};

/**
 * @brief A Permutation-LU factorization for sparse matrices.
 *
 * This class represents a Permutation, Lower-triangular,
 * Upper-triangular factorization for sparse matrices. It
 * also offers a method to solve linear systems
 * represented by sparse matrices. The permutation is not
 * represented by a matrix, but instead by using a
 * vector of indexes swaps.
 *
 * @tparam T any numeric type.
 */
template<typename T>
class PLU_Factorization
{
  Permutation _P; //!< the permutation
  Matrix<T> _L;   //!< the lower-triangular sparse matrix
  Matrix<T> _U;   //!< the upper-triangular sparse matrix

  /**
   * @brief Build the vector of non-zero rows below the diagonal
   *
   * This method builds a vector `V` such that `V[i]` is the set
   * of row indexes `j` greater than `i` and having non-zero
   * value in column `i`.
   *
   * @param M is a sparse matrix.
   * @return the vector of non-zero rows below the diagonal.
   */
  static std::vector<std::set<unsigned int>>
  get_non_zero_below_diag(const Matrix<T> &M)
  {
    std::vector<std::set<unsigned int>> non_zero_below_diag(M.num_of_cols());

    for (auto row_it = std::cbegin(M._matrix); row_it != std::cend(M._matrix);
         ++row_it) {
      for (auto elem_it = std::cbegin(row_it->second);
           elem_it != std::cend(row_it->second); ++elem_it) {

        if (elem_it->first < row_it->first) {
          non_zero_below_diag[elem_it->first].insert(row_it->first);
        }
      }
    }

    return non_zero_below_diag;
  }

  /**
   * @brief Swap two matrix rows.
   *
   * This method swaps the `row_idx`-th row and the first row
   * below the diagonal whose value in column `row_idx` is not
   * zero.
   *
   * @param non_zero_below_diag is the vector such that
   *     `non_zero_below_diag[i]` is the set of row indexes `j`
   *     greater than `i` and having non-zero value in column `i`.
   * @param row_idx is the index of the row that must be swaped down.
   */
  void swap_with_the_leastest_non_zero_row_in_column(
      std::vector<std::set<unsigned int>> &non_zero_below_diag,
      const unsigned int &row_idx)
  {
    // if no value below the diagonal is non-empty
    if (non_zero_below_diag[row_idx].empty()) {

      // the matrix is singular.
      std::domain_error("The parameter is a singular matrix.");
    }

    // get the index of the leastest row having a non-empty value in this
    // column
    unsigned int new_row_idx = *(non_zero_below_diag[row_idx].begin());

    // otherwise swap the current row and the row of the leastest row having a
    // non-empty value in this column
    std::swap(_U._matrix[row_idx], _U._matrix[new_row_idx]);
    std::swap(_L._matrix[row_idx], _L._matrix[new_row_idx]);
    _P.swap(row_idx, new_row_idx);

    // delete the new_row_idx from the sets of non-zero rows below the diagonal
    for (auto elem_it = std::cbegin(_U._matrix[row_idx]);
         elem_it != std::cend(_U._matrix[row_idx]); ++elem_it) {
      non_zero_below_diag[elem_it->first].erase(new_row_idx);
    }

    // add the row_idx to the sets of non-zero rows below the diagonal
    for (auto elem_it = std::cbegin(_U._matrix[new_row_idx]);
         elem_it != std::cend(_U._matrix[new_row_idx]); ++elem_it) {
      if (elem_it->first < new_row_idx) {
        non_zero_below_diag[elem_it->first].insert(new_row_idx);
      }
    }
  }

  /**
   * @brief Solve the linear system $L\cdot x = b$.
   *
   * @param[in] b is the known term of the linear system.
   * @return the vector $x = L^{-1} \cdot b$.
   */
  std::vector<T> solveL(const std::vector<T> &b) const
  {
    std::vector<T> x(b.size());

    for (auto row_it = std::begin(_L._matrix); row_it != std::end(_L._matrix);
         ++row_it) {
      const T &diag_value = std::rbegin(row_it->second)->second;
      if (diag_value == 0) {
        throw std::domain_error("The linear system is underdetermined.");
      }
      T value = b[row_it->first];
      for (auto elem_it = ++std::rbegin(row_it->second);
           elem_it != std::rend(row_it->second); ++elem_it) {
        value -= elem_it->second * x[elem_it->first];
      }
      x[row_it->first] = value / diag_value;
    }

    return x;
  }

  /**
   * @brief Solve the linear system $U\cdot x = b$.
   *
   * @param[in] b is the known term of the linear system.
   * @return the vector $x = U^{-1} \cdot b$.
   */
  std::vector<T> solveU(const std::vector<T> &b) const
  {
    std::vector<T> x(b.size());

    for (auto row_it = std::rbegin(_U._matrix);
         row_it != std::rend(_U._matrix); ++row_it) {
      const T &diag_value = std::begin(row_it->second)->second;
      if (diag_value == 0) {
        throw std::domain_error("The linear system is underdetermined.");
      }
      T value = b[row_it->first];
      for (auto elem_it = ++std::begin(row_it->second);
           elem_it != std::end(row_it->second); ++elem_it) {
        value -= elem_it->second * x[elem_it->first];
      }
      x[row_it->first] = value / diag_value;
    }

    return x;
  }

public:
  /**
   * @brief Create a new PLU_Factorization object
   */
  PLU_Factorization(): _P(), _L(), _U() {}

  /**
   * @brief Copy construct a new PLU_Factorization object
   *
   * @param orig is the model for the new object.
   */
  PLU_Factorization(const PLU_Factorization &orig):
      _P(orig._P), _L(orig._L), _U(orig._U)
  {
  }

  /**
   * @brief Create a new PLU_Factorization object
   *
   * This constructor performs the P-LU factorization of
   * a sparse matrix.
   *
   * @param M is the matrix whose P-LU factorization must be computed.
   */
  PLU_Factorization(const Matrix<T> &M):
      _P(), _L(M.num_of_rows(), M.num_of_rows()), _U(M)
  {
    if (M.num_of_rows() == 0 || M.num_of_cols() == 0) {
      throw std::domain_error("The parameter is an empty matrix.");
    }

    std::vector<std::set<unsigned int>> non_zero_below_diag
        = get_non_zero_below_diag(M);

    for (unsigned int row_idx = 0; row_idx < _L.num_of_rows(); ++row_idx) {
      // set the L diagonal to 1
      _L._matrix[row_idx][row_idx] = 1;
    }

    for (auto row_it = std::begin(_U._matrix); row_it != std::end(_U._matrix);
         ++row_it) {
      const unsigned int &row_idx = row_it->first;

      if (row_idx == M.num_of_cols()) {
        return;
      }

      // get the diagonal element on the row
      auto diag_elem = row_it->second.find(row_idx);

      if (diag_elem == std::end(row_it->second)) {
        // the diagonal element is empty

        if (!non_zero_below_diag[row_idx].empty()) {
          swap_with_the_leastest_non_zero_row_in_column(non_zero_below_diag,
                                                        row_idx);

          diag_elem = _U._matrix[row_it->first].find(row_idx);
        }
      }

      // after the previous conditional statement either the
      // diagonal element is no more zero or all the values
      // below the diagonal in column row_idx are 0 and we
      // can skip to the next row/column

      if (diag_elem != std::end(row_it->second)) {
        // the diagonal element is non-empty
        auto &U_row = _U._matrix[row_idx];

        // for any non-zero row in this column
        for (auto nz_row_it = std::cbegin(non_zero_below_diag[row_idx]);
             nz_row_it != std::cend(non_zero_below_diag[row_idx]);
             ++nz_row_it) {

          auto &nz_row = _U._matrix[*nz_row_it];

          // compute the ratio between the element on the U diagonal and the
          // non-null value below it
          auto ratio = -nz_row[row_idx] / diag_elem->second;

          // for any value in the row row_idx
          for (auto elem_it = std::cbegin(U_row); elem_it != std::cend(U_row);
               ++elem_it) {

            // search for the corresponding element in the considered non-zero
            // row
            auto nz_elem_it = nz_row.find(elem_it->first);

            // if it does not exist, initialize it to the opportune value
            if (nz_elem_it == std::end(nz_row)) {
              nz_row[elem_it->first] = ratio * elem_it->second;
              non_zero_below_diag[elem_it->first].insert(*nz_row_it);
            } else { // otherwise, it exists
              // increase it by the the opportune value
              nz_elem_it->second += ratio * elem_it->second;

              // if now is 0 or it is the element that we aimed to to remove
              // which is different from 0 due to approximation errors
              if (nz_elem_it->second == 0 || nz_elem_it->first == row_idx) {
                // remove the row from the set of non-zero rows in column if
                // necessary
                if (elem_it->first != row_idx) {
                  // postponing row_idx remova to avoid messing up the iterator
                  // on it
                  non_zero_below_diag[elem_it->first].erase(*nz_row_it);
                }

                // delete it from the matrix
                nz_row.erase(nz_elem_it);
              }
            }
          }
          _L._matrix[*nz_row_it][row_idx] = -ratio;
        }
        non_zero_below_diag[row_idx].clear();
      }
    }
  }

  /**
   * @brief Compute the solution of the system $(L \cot U) \cdot x = P^{-1}
   * \cdot b$
   *
   * @param b is the known term of the linear system.
   * @return the vector $x = (U^{-1} \cdot L^{-1} \cdot P^{-1}) \cdot b$.
   */
  std::vector<T> solve(const std::vector<T> &b) const
  {
    if (b.size() != _U.num_of_rows()) {
      throw std::domain_error("Wrong dimensions");
    }

    if (_U.num_of_rows() != _U.num_of_cols()) {
      throw std::domain_error("The factorization is not square and the linear "
                              "system cannot be solved.");
    }

    std::vector<T> Pb = _P(b);

    // Forward substitution for L
    Pb = solveL(Pb);

    // Backward substitution for U
    return solveU(Pb);
  }

  /**
   * @brief Get the permutation
   *
   * @return the factorization permutation
   */
  const Permutation &P() const
  {
    return _P;
  }

  /**
   * @brief Get the lower-triangular matrix
   *
   * @return the factorization lower-triangular matrix
   */
  const Matrix<T> &L() const
  {
    return _L;
  }

  /**
   * @brief Get the upper-triangular matrix
   *
   * @return the factorization upper-triangular matrix
   */
  const Matrix<T> &U() const
  {
    return _U;
  }

  /**
   * @brief Copy a PLU_Factorization object.
   *
   * @param orig is the model that must be copied.
   * @return a reference to the updated object.
   */
  PLU_Factorization<T> &operator=(const PLU_Factorization<T> &orig)
  {
    _P = orig._P;
    _L = orig._L;
    _U = orig._U;

    return *this;
  }

  /**
   * @brief Copy a PLU_Factorization object.
   *
   * @param orig is the model that must be copied.
   * @return a reference to the updated object.
   */
  PLU_Factorization<T> &operator=(PLU_Factorization<T> &&orig)
  {
    std::swap(_P, orig._P);
    std::swap(_L, orig._L);
    std::swap(_U, orig._U);

    return *this;
  }
};

}
/*! @} End of SparseLinearAlgebra group */

namespace std
{
template<typename T>
std::ostream &operator<<(std::ostream &out, const std::vector<T> &v)
{
  out << "[";
  for (auto el_it = std::begin(v); el_it != std::end(v); ++el_it) {
    if (el_it != std::begin(v)) {
      out << ",";
    }
    out << *el_it;
  }
  out << "]";

  return out;
}

template<typename T>
std::ostream &operator<<(std::ostream &out, const std::set<T> &v)
{
  out << "{";
  for (auto el_it = std::begin(v); el_it != std::end(v); ++el_it) {
    if (el_it != std::begin(v)) {
      out << ",";
    }
    out << *el_it;
  }
  out << "}";

  return out;
}

template<typename T>
std::ostream &operator<<(std::ostream &out,
                         const std::vector<std::vector<T>> &A)
{
  for (auto row_it = std::begin(A); row_it != std::end(A); ++row_it) {
    if (row_it != std::begin(A)) {
      out << std::endl;
    }

    out << *row_it;
  }

  return out;
}

template<typename T>
std::ostream &operator<<(std::ostream &os,
                         const SparseLinearAlgebra::Matrix<T> &M)
{
  os << "[ #rows: " << M.num_of_rows() << " #cols: " << M.num_of_cols() << " ";
  for (auto row_it = std::begin(M._matrix); row_it != std::end(M._matrix);
       ++row_it) {
    os << std::endl << " row " << row_it->first << ": [";
    for (auto elem_it = std::begin(row_it->second);
         elem_it != std::end(row_it->second); ++elem_it) {
      if (elem_it != std::begin(row_it->second)) {
        os << ", ";
      }
      os << "col " << elem_it->first << ": " << elem_it->second;
    }
    os << "]";
  }
  os << "]";

  return os;
}

template<typename T>
std::ostream &operator<<(std::ostream &os,
                         const SparseLinearAlgebra::PLU_Factorization<T> &D)
{
  using namespace std;
  os << "{P=" << D.P() << "," << std::endl
     << " L=" << D.L() << "," << std::endl
     << " U=" << D.U() << "}";

  return os;
}

template<typename T>
std::ostream &operator<<(std::ostream &os,
                         const DenseLinearAlgebra::PLU_Factorization<T> &D)
{
  using namespace std;
  os << "{P=" << D._P << "," << std::endl << " LU=" << D._LU << "}";

  return os;
}
}
#endif // LINEAR_ALGEBRA_H
