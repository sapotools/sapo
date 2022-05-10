/**
 * @file LinearAlgebra.h
 * @author Alberto Casagrande <acasagrande@units.it>
 * @brief Contains linear algebra classes and functions code
 * @version 0.1
 * @date 2021-11-20
 * 
 * @copyright Copyright (c) 2021-2022
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

namespace LinearAlgebra {

/**
 * @brief An alias for vector type
 * 
 * @tparam T is the scalar value type
 */
template<typename T>
using Vector = std::vector<T>;

/**
 * @brief Approximate all the values in a vector up to given decimal digit
 *
 * @tparam T is a numerical type
 * @param v is a the vector
 * @param decimal is the number of decimal number in the result
 * @return the given vector approximated up to `decimal` decimal digits
 */
template<typename T>
inline Vector<T> approx(const Vector<T> &v,
                        const unsigned short decimal
                        = std::numeric_limits<T>::digits10);

/**
 * @brief Approximate all the values in a vector up to given decimal digit
 *
 * This function can be mainly used to speed up the computation.
 *
 * @tparam T is a numerical type
 * @param v is a the vector
 * @param decimal is the number of decimal number in the result
 * @return the given vector approximated up to `decimal` decimal digits
 */
template<typename T>
inline Vector<T> approx(Vector<T> &&v, const unsigned short decimal
                                        = std::numeric_limits<T>::digits10);

/**
 * @brief Approximate all the values in a vector up to given decimal digit
 *
 * This function can be mainly used to speed up the computation.
 *
 * @param v is a the vector
 * @param decimal is the number of decimal number in the result
 * @return the given vector approximated up to `decimal` decimal digits
 */
template<>
inline Vector<double> approx<double>(Vector<double> &&v,
                                     const unsigned short decimal)
{
  if (decimal == std::numeric_limits<double>::digits10) {
    return std::move(v);
  }

  const double mult = std::pow(10, decimal);
  for (auto it = std::begin(v); it != std::end(v); ++it) {
    *it = floor(*it * mult) / mult;
  }

  return std::move(v);
}

/**
 * @brief Approximate all the values in a vector up to given decimal digit
 *
 * This function can be mainly used to speed up the computation.
 *
 * @param v is a the vector
 * @param decimal is the number of decimal number in the result
 * @return the given vector approximated up to `decimal` decimal digits
 */
template<>
inline Vector<double> approx<double>(const Vector<double> &v,
                                     const unsigned short decimal)
{
  return approx<double>(Vector<double>(v), decimal);
}

/**
 * @brief Compute the Euclidean norm of a vector
 *
 * @tparam T is a numeric type supporting the `sqrt` function
 * @param v is a vector
 * @return the Euclidean norm of the given vector
 */
template<typename T>
T norm_2(const Vector<T> &v)
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
 * @tparam T is a numeric type supporting the `abs` function
 * @param v is a vector
 * @return the Manhattan norm of the given vector
 */
template<typename T>
T norm_1(const Vector<T> &v)
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
 * @tparam T is a numeric type supporting the `abs` function
 * @param v is a vector
 * @return the maximum norm of the given vector
 */
template<typename T>
T norm_infinity(const Vector<T> &v)
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
 * @param[in] orig the vector of values to be complementated
 * @return A vector containing the complementaries of the parameter values
 */
template<typename T>
Vector<T> operator-(const Vector<T> &orig)
{
  Vector<T> res = orig;

  transform(res.begin(), res.end(), res.begin(), std::negate<T>());

  return res;
}

/**
 * @brief Check whether two vectors are the same
 *
 * @tparam T is the scalar value type
 * @param a is the first vector to be compared
 * @param b is the second vector to be compared
 * @return `true` if and only if the two vectors are the same
 */
template<typename T>
bool operator==(const Vector<T> &a, const Vector<T> &b)
{
  if (a.size() != b.size()) {
    throw std::domain_error("The two vectors have different dimensions");
  }

  auto b_it = std::cbegin(b);
  for (auto a_it = std::cbegin(a);
       a_it != std::cend(a); ++a_it, ++b_it) {
    if (*a_it != *b_it) {
      return false;
    }
  }

  return true;
}

/**
 * @brief Check whether two vectors differ
 *
 * @tparam T is the scalar value type
 * @param a is the first vector to be compared
 * @param b is the second vector to be compared
 * @return `true` if and only if the two vectors differ
 */
template<typename T>
bool operator!=(const Vector<T> &a, const Vector<T> &b)
{
  return !(a == b);
}

/**
 * Compute the complementary of a vector of values.
 *
 * @param[in] v the vector of values to be complementated
 * @return A vector containing the complementaries of the parameter values
 */
template<typename T>
Vector<T> operator-(Vector<T> &&v)
{
  transform(v.begin(), v.end(), v.begin(), std::negate<T>());

  return std::move(v);
}

/**
 * @brief Compute the element-wise sum of two vectors
 *
 * @tparam T any numeric type
 * @param a is the first vector to be added
 * @param b is the second vector to be added
 * @return the element-wise sum of the two parameters
 */
template<typename T>
Vector<T> operator+(const Vector<T> &a, const Vector<T> &b)
{
  if (a.size() != b.size()) {
    throw std::domain_error("The two vectors have different dimensions");
  }

  Vector<T> res(a.size());

  for (unsigned int i = 0; i < res.size(); ++i) {
    res[i] = a[i] + b[i];
  }

  return res;
}

/**
 * @brief Compute the element-wise sum of two vectors
 *
 * @tparam T any numeric type
 * @param a is the first vector to be added
 * @param b is the second vector to be added
 * @return the element-wise sum of the two parameters
 */
template<typename T>
Vector<T> operator+(Vector<T> &&a, const Vector<T> &b)
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
 * @brief Compute the element-wise sum of two vectors
 *
 * @tparam T any numeric type
 * @param a is the first vector to be added
 * @param b is the second vector to be added
 * @return the element-wise sum of the two parameters
 */
template<typename T>
inline Vector<T> operator+(const Vector<T> &a, Vector<T> &&b)
{
  return operator+(b, a);
}

/**
 * @brief Compute the element-wise difference of two vectors
 *
 * @tparam T any numeric type
 * @param a is the vector from which the second parameter must be subtracted
 * @param b is the vector that must be subtracted
 * @return the element-wise difference of the second parameter from the first
 *         one
 */
template<typename T>
Vector<T> operator-(const Vector<T> &a, const Vector<T> &b)
{
  if (a.size() != b.size()) {
    throw std::domain_error("The two vectors have different dimensions");
  }

  Vector<T> res(a.size());

  for (unsigned int i = 0; i < res.size(); ++i) {
    res[i] = a[i] - b[i];
  }

  return res;
}

/**
 * @brief Compute the element-wise difference of two vectors
 *
 * @tparam T any numeric type
 * @param a is the vector from which the second parameter must be subtracted
 * @param b is the vector that must be subtracted
 * @return the element-wise difference of the second parameter from the first
 *         one
 */
template<typename T>
Vector<T> operator-(Vector<T> &&a, const Vector<T> &b)
{
  if (a.size() != b.size()) {
    throw std::domain_error("The two vectors have different dimensions");
  }

  for (unsigned int i = 0; i < a.size(); ++i) {
    a[i] -= b[i];
  }

  return std::move(a);
}

/**
 * @brief Compute the element-wise difference of two vectors
 *
 * @tparam T any numeric type
 * @param a is the vector from which the second parameter must be subtracted
 * @param b is the vector that must be subtracted
 * @return the element-wise difference of the second parameter from the 
 *         first one
 */
template<typename T>
Vector<T> operator-(const Vector<T> &a, Vector<T> &&b)
{
  if (a.size() != b.size()) {
    throw std::domain_error("The two vectors have different dimensions");
  }

  for (unsigned int i = 0; i < a.size(); ++i) {
    b[i] -= a[i];
    b[i] *= -1;
  }

  return std::move(b);
}

/**
 * @brief Compute the element-wise scalar product
 *
 * @tparam T any numeric type
 * @param s is a scalar value
 * @param v is a vector
 * @return the element-wise scalar product \f$s * \f$
 */
template<typename T>
Vector<T> operator*(const T s, const Vector<T> &v)
{
  Vector<T> res(v.size());

  for (unsigned int i = 0; i < res.size(); ++i) {
    res[i] = s * v[i];
  }

  return res;
}

/**
 * @brief Compute the element-wise scalar product
 *
 * @tparam T any numeric type
 * @param v is a vector
 * @param s is a scalar value
 * @return the element-wise scalar product \f$v * \f$
 */
template<typename T>
Vector<T> operator*(const Vector<T> &v, const T s)
{
  return operator*(s, v);
}

/**
 * @brief Compute the vector product
 *
 * @tparam T any numeric type
 * @param v1 is a vector
 * @param v2 is a vector
 * @return the vector product \f$v1 \cdot v\f$
 */
template<typename T>
T operator*(const Vector<T> &v1, const Vector<T> &v2)
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
 * @brief Compute the Hadamard product of two vectors
 *
 * This method computes the Hadamard product of two vectors
 * \f$v_1 \circ v_2\f$, i.e., it returns a vector whose elements
 * are obtained by element-wise multiplying of the
 * parameters.
 * \f$(v_1 \circ v_2)[i] = v_1[i]*v_2[i]\f$
 *
 * @tparam T any numeric type
 * @param v1 is a vector
 * @param v2 is a vector
 * @return the vector product \f$v1 \circ v2\f$
 */
template<typename T>
Vector<T> H_prod(const Vector<T> &v1, const Vector<T> &v2)
{
  if (v1.size() != v2.size()) {
    throw std::domain_error("The two vectors have not the same dimensions.");
  }

  Vector<T> res(v1.size());
  for (unsigned int i = 0; i < v1.size(); ++i) {
    res[i] = v1[i] * v2[i];
  }

  return res;
}

/**
 * @brief Compute the element-wise scalar product
 *
 * @tparam T any numeric type
 * @param s is a scalar value
 * @param v is a vector
 * @return the element-wise scalar product \f$s * \f$
 */
template<typename T>
Vector<T> operator*(const T s, Vector<T> &&v)
{
  for (auto v_it = std::begin(v); v_it != std::end(v); ++v_it) {
    *v_it *= s;
  }

  return v;
}

/**
 * @brief Compute the element-wise scalar product
 *
 * @tparam T any numeric type
 * @param v is a vector
 * @param s is a scalar value
 * @return the element-wise scalar product \f$v * \f$
 */
template<typename T>
inline Vector<T> operator*(Vector<T> &&v, const T s)
{
  return operator*(s, v);
}

/**
 * @brief Compute the element-wise scalar division
 *
 * @tparam T any numeric type
 * @param v is a vector
 * @param s is a scalar value
 * @return the element-wise scalar division \f$v/\f$
 */
template<typename T>
Vector<T> operator/(const Vector<T> &v, const T s)
{
  if (s==0) {
    throw std::domain_error("Division by 0");
  }

  Vector<T> res(v.size());

  for (unsigned int i = 0; i < res.size(); ++i) {
    res[i] = v[i] / s;
  }

  return res;
}

/**
 * @brief Compute the element-wise scalar division
 *
 * @tparam T any numeric type
 * @param v is a vector
 * @param s is a scalar value
 * @return the element-wise scalar division \f$v/\f$
 */
template<typename T>
Vector<T> operator/(Vector<T> &&v, const T s)
{
  if (s==0) {
    throw std::domain_error("Division by 0");
  }

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
inline T angle(const Vector<T> &v1, const Vector<T> &v2)
{
  return acos(v1 * v2 / (norm_2(v1) * norm_2(v2)));
}

/**
 * @brief A class to represent pemutations
 *
 * Permutation are square matrices whose row and column
 * sums are one. However, storing a permutation for
 * n-dim vectors by using this natural representation
 * requires \f$n^2\f$ integer numbers and its
 * application costs \f$O(n^2)\f$.
 * This class stores permutation is a more efficient
 * way representing only swapped element positions in
 * a map from integer to integer. The application cost
 * belongs to \f$O(m*log m)\f$ where \f$m\f$ is the
 * number of swapped positions.
 */
class Permutation : public std::map<unsigned int, unsigned int>
{
public:
  /**
   * @brief Create a new Permutation object
   */
  Permutation(): std::map<unsigned int, unsigned int>() {}

  /**
   * @brief Copy constructor for Permutation object
   *
   * @param orig
   */
  Permutation(const Permutation &orig): 
    std::map<unsigned int, unsigned int>(orig) {}

  /**
   * @brief Swaps two positions
   *
   * This method swaps two positions in the permutation.
   *
   * @param i first position to be swapped
   * @param j second position to be swapped
   * @return The updated permutation
   */
  Permutation &swap(const unsigned int i, const unsigned int j)
  {
    // for both the positions
    for (const unsigned int &pos: {i, j}) {
      // if the position is not stored among the swaps yet
      if (this->find(pos) == this->end()) {

        // store it
        this->operator[](pos) = pos;
      }
    }

    // swap the two positions
    std::swap(this->at(i), this->at(j));

    return *this;
  }

  /**
   * @brief Apply the permutation to a vector
   *
   * This method does not work in-place and
   * creates a swapped-version of the input
   * vector to be returned.
   *
   * @tparam T is the scalar value type
   * @param v is the vector to be swapped
   * @return the swapped vector
   */
  template<typename T>
  Vector<T> operator()(const Vector<T> &v) const
  {
    Vector<T> pv = v;

    for (auto it = this->begin(); it != this->end(); ++it) {
      if (it->first >= v.size() || it->second >= v.size()) {
        throw std::domain_error(
            "Permutation and vector dimensions are not compatible.");
      }
      pv[it->first] = v[it->second];
    }

    return pv;
  }

  /**
   * @brief Assignment operator
   *
   * @param P is the template permutation
   * @return a reference to the updated object
   */
  Permutation &operator=(const Permutation &P)
  {
    static_cast<std::map<unsigned int, unsigned int>>(*this) = 
      static_cast<std::map<unsigned int, unsigned int>>(P);

    return *this;
  }

  /**
   * @brief Print a permutation in a stream
   * 
   * @param out is the output stream
   * @param P is the permutation to be printed
   * @return a reference to the output stream
   */
  friend inline std::ostream &operator<<(std::ostream &out,
                                         const Permutation &P);
};

/**
 * @brief Establish whether two vectors are linearly dependent
 *
 * @tparam T is the type of the scalar values
 * @param v1 is the first vector to be tested
 * @param v2 is the second vector to be tested
 * @return `true` if and only if the two vectors are linearly dependent
 */
template<typename T>
bool are_linearly_dependent(const Vector<T> &v1, const Vector<T> &v2)
{
  if (v1.size() != v2.size()) {
    throw std::domain_error("The two vectors have not the same dimensions.");
  }

  if (v1.size() == 0) {
    throw std::domain_error("The two vectors must have size greater than 0.");
  }

  // initialize the candidate index to be the first one
  // different from zero in v1 or v2
  unsigned int non_zero_idx = v1.size();
  for (unsigned int i = 0; i < v1.size(); ++i) {

    // if non zero values has not been discovered yet and
    // at least one of the two elements in the two vectors
    // differs from zero
    if (non_zero_idx == v1.size()) {
      if (v1[i] != 0 || v2[i] != 0) {
        if (v1[i] == 0 || v2[i] == 0) {
          return false;
        }

        // set the first non-zero index
        non_zero_idx = i;
      }
    } else {
      if ((v1[i] == 0 && v2[i] != 0)||(v1[i] != 0 && v2[i] == 0)||
          (v1[i] * v2[non_zero_idx] != v2[i] * v1[non_zero_idx])) {

        // return that their are not linearly dependent
        return false;
      }
    }
  }

  if (non_zero_idx==v1.size()) {
    return false;
  }

  return true;
}

/**
 * @brief Compute the ratio between two linearly dependent vectors
 *
 * @tparam T is the type of the scalar values
 * @param v1 is the first vector
 * @param v2 is the second vector
 * @return The ratio between two linearly dependent vectors
 */
template<typename T>
T operator/(const Vector<T> &v1, const Vector<T> &v2)
{
  if (v1.size() != v2.size()) {
    throw std::domain_error("The two vectors have not the same dimensions.");
  }

  // initialize the candidate index to be the first one
  // different from zero in v1 or v2
  unsigned int non_zero_idx = v1.size();
  for (unsigned int i = 0; i < v1.size(); ++i) {

    // if non zero values has not been discovered yet and
    // at least one of the two elements in the two vectors
    // differs from zero
    if (non_zero_idx == v1.size()) {
      if (v1[i] != 0 || v2[i] != 0) {
        if (v1[i] == 0 || v2[i] == 0) {
          throw std::domain_error("The two vectors are not linealy dependent.");
        }

        // set the first non-zero index
       non_zero_idx = i;
      }
    } else {

      if ((v1[i] == 0 && v2[i] != 0)||(v1[i] != 0 && v2[i] == 0)||
          (v1[i] * v2[non_zero_idx] != v2[i] * v1[non_zero_idx])) {
        throw std::domain_error("The two vectors are not linealy dependent.");
      }
    }
  }

  if (non_zero_idx==v1.size()) {
    throw std::domain_error("The two vectors are not linealy dependent.");
  }

  return v1[non_zero_idx] / v2[non_zero_idx];
}

/*!
 *  \addtogroup Dense
 *  @{
 */

//! Linear algebra for dense matrices
namespace Dense
{

/**
 * @brief An alias for dense matrices
 * 
 * @tparam T is the scalar value type
 */
template<typename T>
using Matrix = std::vector<Vector<T>>;

/**
 * @brief Transpose a matrix
 * 
 * @tparam T is the scalar value type
 * @param A is the matrix to be transposed
 * @return the transposed matrix
 */
template<typename T>
Matrix<T> transpose(const Matrix<T>& A)
{
  size_t num_of_rows = A.size();
  size_t num_of_cols = (num_of_rows == 0 ? 0: A[0].size());
  Matrix<T> TA(num_of_cols, Vector<T>(num_of_rows));

  for (size_t i=0; i<num_of_rows; ++i) {
    for (size_t j=0; j<num_of_cols; ++j) {
      TA[j][i] = A[i][j];
    }
  }

  return TA;
}

/**
 * @brief Compute the row-column matrix-vector multiplication
 *
 * @tparam T is a numeric type
 * @param A is a dense matrix
 * @param v is a vector
 * @return the row-column matrix-vector multiplication \f$A \cdot \f$
 */
template<typename T>
Vector<T> operator*(const Matrix<T> &A, const Vector<T> &v)
{
  if (((A.size() != 0 && A.front().size() == 0) || (A.size() == 0))
      && v.size() == 0) {
    return Vector<T>();
  }

  if (A.size() == 0 || A.front().size() != v.size()) {
    throw std::domain_error("The matrix and the vector are not compatible for the "
                            "matrix-vector multiplication.");
  }

  Vector<T> res(A.size(), 0);

  auto res_it = std::begin(res);
  for (auto A_row_it = std::begin(A);A_row_it != std::end(A); ++A_row_it, ++res_it) {
    auto A_elem_it = std::begin(*A_row_it);
    for (auto v_it = std::begin(v); v_it != std::end(v); ++v_it, ++A_elem_it) {
      *res_it += (*A_elem_it) * (*v_it);
    }
  }

  return res;
}

/**
 * @brief Compute the row-column matrix-matrix multiplication
 *
 * @tparam T is a numeric type
 * @param A is a dense matrix
 * @param B is a dense matrix
 * @return the row-column matrix-matrix multiplication \f$A \cdot \f$
 */
template<typename T>
Matrix<T> operator*(const Matrix<T> &A, const Matrix<T> &B)
{
  if (((A.size() != 0 && A.front().size() == 0) || (A.size() == 0))
      && B.size() == 0) {
    return Matrix<T>();
  }

  if (A.front().size() != B.size()) {
    throw std::domain_error("The two matrices are not compatible for the "
                            "matrix-matrix multiplication.");
  }

  Matrix<T> res(A.size(), Vector<T>(B.front().size(), 0));

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
 * @brief A LUP factorization for dense matrices
 *
 * This class represents a Permutation, Lower-triangular,
 * Upper-triangular factorization for dense matrices. It
 * also offers a method to solve linear systems
 * represented by dense matrices. The permutation is not
 * represented by a matrix, but instead by using a
 * vector of indexes swaps.
 *
 * @tparam T any numeric type
 */
template<typename T>
class LUP_Factorization
{
  Permutation _P; //!< the permutation
  Matrix<T> _LU;  //!< the L and U matrices

  /**
   * @brief Solve the linear system \f$L\cdot x = \f$
   *
   * @param[in] b is the known term of the linear system
   * @return the vector \f$x = L^{-1} \cdot \f$
   */
  Vector<T> solveL(const Vector<T> &b) const
  {
    Vector<T> x(b.size());

    for (unsigned int j = 0; j < _LU.size(); ++j) {
      T value = b[j];
      const Vector<T> &row = _LU[j];
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
   * @brief Solve the linear system \f$U\cdot x = \f$
   *
   * @param[in] b is the known term of the linear system
   * @return the vector \f$x = U^{-1} \cdot \f$
   */
  Vector<T> solveU(const Vector<T> &b) const
  {
    Vector<T> x(b.size());

    for (unsigned int rj = 0; rj < _LU.size(); ++rj) {
      const unsigned int j = _LU.size() - rj - 1;
      T value = b[j];
      const Vector<T> &row = _LU[j];

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
   * @brief Create a new LUP_Factorization object
   *
   */
  LUP_Factorization(): _P(), _LU() {}

  /**
   * @brief Create a new LUP_Factorization object
   *
   * @param orig is the model for the new object
   */
  LUP_Factorization(const LUP_Factorization &orig): _P(orig._P), _LU(orig._LU)
  {
  }

  /**
   * @brief Create a new LUP_Factorization object
   *
   * This constructor performs the LUP factorization of
   * a dense matrix.
   *
   * @param M is the matrix whose LUP factorization must be computed
   */
  LUP_Factorization(const Matrix<T> &M): _P(), _LU(M)
  {
    if (M.size() == 0) {
      throw std::domain_error("The parameter is an empty matrix.");
    }

    for (unsigned int j = 0; j < _LU.size(); ++j) {
      if (j == M.front().size()) {
        return;
      }

      // find the non-null value below A[j-1][j] whose absolute value
      // is the nearest to 1
      unsigned int k = _LU.size();
      T delta = 1;
      unsigned int l=j;
      while (l<_LU.size() && delta > 0) {
        if (_LU[l][j] != 0) {
          T new_delta = std::abs(1-std::abs(_LU[l][j]));
          if (k==_LU.size()) {
            delta = new_delta;
            k = l;
          } else {
            if (new_delta < delta) {
              delta = new_delta;
              k = l;
            }
          }
        }
        ++l;
      }

      // if it does not exist, skip to the next row/column
      if (k < _LU.size()) {
        // otherwise, swap rows
        if (j != k) {
          std::swap(_LU[j], _LU[k]);
          _P.swap(j, k);
        }

        // nullify all the values below A[j][j]
        const Vector<T> &row_j = _LU[j];
        const T &v = row_j[j];
        k = j + 1;
        while (k < _LU.size()) {
          Vector<T> &row_k = _LU[k];
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
   * @brief Compute the solution of \f$(L\cdot U)\cdot x=P^{-1}\cdot b\f$
   *
   * @param b is the known term of the linear system
   * @return the vector \f$x = (U^{-1} \cdot L^{-1} \cdot P^{-1}) \cdot b\f$
   */
  Vector<T> solve(const Vector<T> &b) const
  {
    if (b.size() != _LU.size()) {
      throw std::domain_error("Wrong dimensions");
    }

    if (_LU.size() != _LU.front().size()) {
      throw std::domain_error("The factorization is not square and the linear "
                              "system cannot be solved.");
    }

    Vector<T> Pb = _P(b);

    // Forward substitution for L
    Pb = solveL(Pb);

    // Backward substitution for U
    return solveU(Pb);
  }

  /**
   * @brief Copy a LUP_Factorization object
   *
   * @param orig is the model that must be copied
   * @return a reference to the updated object
   */
  LUP_Factorization<T> &operator=(const LUP_Factorization<T> &orig)
  {
    _P = orig._P;
    _LU = orig._LU;

    return *this;
  }

  /**
   * @brief Copy a LUP_Factorization object
   *
   * @param orig is the model that must be copied
   * @return a reference to the updated object
   */
  LUP_Factorization<T> &operator=(LUP_Factorization<T> &&orig)
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
   * @return the factorization lower-triangular matrix
   */
  Matrix<T> L() const
  {
    Matrix<T> res(_LU.size(), Vector<T>(_LU.size(), 0));

    for (unsigned int row_idx = 0; row_idx < _LU.size(); ++row_idx) {
      const Vector<T> &row = _LU[row_idx];
      Vector<T> &res_row = res[row_idx];
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
   * @return the factorization upper-triangular matrix
   */
  Matrix<T> U() const
  {
    Matrix<T> res(_LU.size(), Vector<T>(_LU.front().size(), 0));

    for (unsigned int row_idx = 0; row_idx < _LU.size(); ++row_idx) {
      const Vector<T> &row = _LU[row_idx];
      Vector<T> &res_row = res[row_idx];
      for (unsigned int col_idx = row_idx; col_idx < row.size(); ++col_idx) {
        res_row[col_idx] = row[col_idx];
      }
    }

    return res;
  }

  /**
   * @brief Print a dense LUP factorization in a stream
   * 
   * @tparam E is the type of the scalar values
   * @param out is the output stream
   * @param F is the LUP factorization to be printed
   * @return a reference to the output stream
   */
  template<typename E>
  friend std::ostream &std::operator<<(std::ostream &out,
                                       const LUP_Factorization<E> &F);
  /**
   * @brief Compute the rank of a dense matrix
   *
   * @tparam E is the scalar value type
   * @param A is the matrix whose rank must be computed
   * @return The rank of the matrix `A`
   */
  template<typename E>
  friend unsigned int rank(const Matrix<E> &A);
};

/**
 * @brief Compute the rank of a dense matrix
 *
 * @tparam T is the scalar value type
 * @param A is the matrix whose rank must be computed
 * @return The rank of the matrix `A`
 */
template<typename T>
unsigned int rank(const Matrix<T> &A)
{
  if (A.size() == 0) {
    return 0;
  }

  LUP_Factorization<T> fact(A);

  unsigned int rank = 0;
  const unsigned max_diag = std::min(A.size(), A[0].size());
  for (unsigned int i = 0; i < max_diag; ++i) {
    if (fact._LU[i][i] != 0) {
      rank++;
    }
  }

  return rank;
}

/**
 * @brief Compute the determinant of a matrix
 * 
 * The determinant is evaluated by using three steps: first
 * of all, the LUP decomposition of the matrix is computed,
 * then the product of the elements in \f$U\f$ main diagonal
 * is computed, and, finally, this product is multiplied by 
 * \f$-1\f$ if the number of swaps in \f$P\f$ are odd.
 * 
 * @tparam T is the scalar domain type
 * @param A is the matrix whose determinant must be computed
 * @return the determinant of the matrix `A`
 */
template<typename T>
T determinant(const Matrix<T> &A)
{
  if (A.size()==0 || A.front().size()==0) {
    throw std::domain_error("Determinant of a 0x0-matrix not supported");
  }
  
  // Compute the factorization
  LUP_Factorization<T> fact(A);

  // multiply the elements in U main diagonal
  T det = 1;
  unsigned int elem_in_diag = std::min(A.size(), 
                                       A.front().size());
  for (unsigned int i=0; i<elem_in_diag; ++i)
  {
    det *= fact.U()[i][i];

    if (det == 0) {
      return 0;
    }
  }

  // Evaluate the number of row swaps and determine
  // the {1, -1} multiplier 
  auto P = fact.P();
  std::vector<bool> swaps(A.size(), true);
  for (auto P_it = P.begin(); P_it != P.end(); ++P_it) {
    unsigned int j=P_it->first;
    while (swaps[j]) {
      swaps[j] = false;
      j = P[j];
      if (j != P_it->first) {
        det = -1*det;
      }
    }
  }

  return det;
}

/**
 * @brief Compute the inverse matrix
 * 
 * @tparam T is the type of the scalar values
 * @param A is the matrix to be inverted
 * @return a matrix \f$A^{-1}\f$ such that \f$A\cdot A^{-1} = I\f$
 */
template<typename T>
Matrix<T> inverse(Matrix<T> &A)
{
  if (A.size()==0 || A.front().size()==0) {
    throw std::domain_error("Inverse of a 0x0-matrix not supported");
  }
  
  Matrix<T> Ainv;

  // Compute the factorization
  LUP_Factorization<T> fact(A);
  std::vector<T> vi(A.size(),0);

  for (unsigned int i=0; i<A.size(); ++i) {
    vi[i] = 1;
    Ainv.push_back(fact.solve(vi));
    vi[i] = 0;
  }

  return transpose(Ainv);
}

}
/*! @} End of DenseLinearAlgebra group */

/**
 * @brief Compute the maximum norm of a dense matrix
 *
 * @tparam T is a numeric type supporting the `abs` function
 * @param A is a dense matrix
 * @return the maximum norm of the given dense matrix
 */
template<typename T>
T norm_infinity(const Dense::Matrix<T> &A)
{
  T norm = 0;

  for (auto row_it = std::begin(A); row_it != std::end(A); ++row_it) {
    norm = std::max(norm, norm_infinity(*row_it));
  }

  return norm;
}

/**
 * @brief Element-wise scalar multiplication for matrices
 * 
 * @tparam T is the type of scalar values
 * @param A is a matrix
 * @param scalar is a scalar value
 * @return the matrix \f$A*\textrm{scalar}\f$
 */
template<typename T>
Dense::Matrix<T> operator*(Dense::Matrix<T> &&A, const T scalar)
{
  for (auto row_it = std::begin(A); row_it != std::end(A); ++row_it) {
    for (auto elem_it = std::begin(*row_it); elem_it != std::end(*row_it); ++elem_it) {
      *elem_it *= scalar;
    }
  }

  return A;
}

/**
 * @brief Element-wise scalar multiplication for matrices
 * 
 * @tparam T is the type of scalar values
 * @param A is a matrix
 * @param scalar is a scalar value
 * @return the matrix \f$A*\textrm{scalar}\f$
 */
template<typename T>
Dense::Matrix<T> operator*(const Dense::Matrix<T>& A, const T scalar)
{
  return operator*(Dense::Matrix<T>(A), scalar); 
}

/**
 * @brief Element-wise scalar multiplication for matrices
 * 
 * @tparam T is the type of scalar values
 * @param scalar is a scalar value
 * @param A is a matrix
 * @return the matrix \f$\textrm{scalar}*A\f$
 */
template<typename T>
Dense::Matrix<T> operator*(const T scalar, const Dense::Matrix<T>& A)
{
  return A*scalar;
}

/**
 * @brief Element-wise scalar multiplication for matrices
 * 
 * @tparam T is the type of scalar values
 * @param scalar is a scalar value
 * @param A is a matrix
 * @return the matrix \f$\textrm{scalar}*A\f$
 */
template<typename T>
Dense::Matrix<T> operator*(const T scalar, Dense::Matrix<T>&& A)
{
  return operator*(std::move(A),scalar);
}

/**
 * @brief Element-wise scalar division for matrices
 *
 * @tparam T is the type of scalar values
 * @param A is a matrix
 * @param scalar is a scalar value
 * @return the matrix \f$A/\textrm{scalar}\f$
 */
template<typename T>
Dense::Matrix<T> operator/(Dense::Matrix<T> &&A, const T scalar)
{
  if (scalar==0) {
    throw std::domain_error("Division by 0");
  }

  for (auto row_it = std::begin(A); row_it != std::end(A); ++row_it) {
    for (auto elem_it = std::begin(*row_it); elem_it != std::end(*row_it); ++elem_it) {
      *elem_it /= scalar;
    }
  }

  return A;
}

/**
 * @brief Element-wise scalar division for matrices
 * 
 * @tparam T is the type of scalar values
 * @param A is a matrix
 * @param scalar is a scalar value
 * @return the matrix \f$A/\textrm{scalar}\f$
 */
template<typename T>
Dense::Matrix<T> operator/(const Dense::Matrix<T>& A, const T scalar)
{
  return operator/(Dense::Matrix<T>(A), scalar); ; 
}

/*!
 *  \addtogroup Sparse
 *  @{
 */

//! Linear algebra for sparse matrices
namespace Sparse
{

/**
 * @brief A class to represent sparse matrices
 *
 * The matrix is internally represented as a map 
 * index-row type and each row is a index-element 
 * type.
 * 
 * @tparam T any numeric type
 */
template<typename T>
class Matrix
{
public:
  typedef std::map<unsigned int, T> RowType; //!< the index-row map

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
     * @param value is value that must be assign to the referenced element
     * @return a reference to the assigned matrix element
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
     * @brief The assignment method
     *
     * This method assign a value to the referenced matrix element. If the
     * element already exists, it reassigns it with a new value. If the
     * matrix element does not exists yet, it creates the referenced
     * element (and its row if necessary) and assigns it with `value`.
     *
     * @param value is value that must be assign to the referenced element
     * @return a reference to the assigned matrix element
     */
    T &operator+=(const T &value)
    {
      // TODO: delete the row if it only contains zeros
      auto row_it = A._matrix.find(row_idx);

      if (row_it == std::end(A._matrix)) {
        A._matrix[row_idx] = Matrix<T>::RowType();

        row_it = A._matrix.find(row_idx);
      }

      row_it->second[col_idx] = row_it->second[col_idx] + value;

      return row_it->second[col_idx];
    }


    /**
     * @brief The assignment method
     *
     * This method assign a value to the referenced matrix element. If the
     * element already exists, it reassigns it with a new value. If the
     * matrix element does not exists yet, it creates the referenced
     * element (and its row if necessary) and assigns it with `value`.
     *
     * @param value is value that must be assign to the referenced element
     * @return a reference to the assigned matrix element
     */
    T &operator-=(const T &value)
    {
      // TODO: delete the row if it only contains zeros
      auto row_it = A._matrix.find(row_idx);

      if (row_it == std::end(A._matrix)) {
        A._matrix[row_idx] = Matrix<T>::RowType();

        row_it = A._matrix.find(row_idx);
      }

      row_it->second[col_idx] = row_it->second[col_idx] - value;

      return row_it->second[col_idx];
    }

    /**
     * @brief The typecast operator for the type T
     *
     * @return The value of the matrix element
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
     * @param orig is the model for the new object
     */
    _row_ref_type(const _row_ref_type &orig): A(orig.A), row_idx(orig.row_idx)
    {
    }

    /**
     * @brief Return a reference of the element in the `col_idx`-th position of
     * the row
     *
     * @param col_idx is the index of the aimed value
     * @return a reference of the element in position `col_idx`
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
     * @param orig is the model for the new object
     */
    _const_row_ref_type(const _const_row_ref_type &orig):
        A(orig.A), row_idx(orig.row_idx)
    {
    }

    /**
     * @brief Return the value in the `col_idx`-th position of the row
     *
     * @param col_idx is the index of the aimed value
     * @return a copy of the value in position `col_idx`
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
   * @return the RowType representing `row`
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
   * @param num_of_rows is the number of rows in the new matrix
   * @param num_of_cols is the number of columns in the new matrix
   */
  Matrix(const unsigned int num_of_rows, const unsigned int num_of_cols):
      _matrix(), _num_of_rows(num_of_rows), _num_of_cols(num_of_cols)
  {
  }

  /**
   * @brief Create a new Matrix object
   *
   * @param orig is the model for the new matrix
   * @param up_to_row is the number of rows to be considered in `orig`
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
   * @param orig is the model for the new matrix
   */
  Matrix(const std::vector<std::vector<T>> &orig): Matrix(orig, orig.size()) {}

  /**
   * @brief Clone Matrix object
   *
   * @param orig is the model for the new object
   */
  Matrix(const Matrix<T> &orig):
      _matrix(orig._matrix), _num_of_rows(orig._num_of_rows),
      _num_of_cols(orig._num_of_cols)
  {
  }

  /**
   * @brief Add a row to the matrix
   *
   * @param row_idx is the index of the new row
   * @param row is the new row
   * @return a reference to the updated object
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
   * @brief Add a row to the matrix
   *
   * @param row_idx is the index of the new row
   * @param row is the new row
   * @return a reference to the updated object
   */
  inline Matrix<T> &add_row(unsigned row_idx, const std::vector<T> &row)
  {
    return add_row(row_idx, get_RowType(row));
  }

  /**
   * @brief Add a row as last row of the matrix
   *
   * @param row is the new row
   * @return a reference to the updated object
   */
  inline Matrix<T> &add_row(const RowType &row)
  {
    return add_row(this->num_of_rows(), row);
  }

  /**
   * @brief Add a row as last row of the matrix
   *
   * @param row is the new row
   * @return a reference to the updated object
   */
  inline Matrix<T> &add_row(const std::vector<T> &row)
  {
    return add_row(this->num_of_rows(), get_RowType(row));
  }

  /**
   * @brief Get a constant reference to one of the matrix rows
   *
   * @param row_idx is the index of the aimed row
   * @return the constant reference of the `row_idx`-th row
   */
  typename Matrix<T>::_const_row_ref_type
  operator[](const unsigned int row_idx) const
  {
    return _const_row_ref_type(*this, row_idx);
  }

  /**
   * @brief Get a reference to one of the matrix rows
   *
   * @param row_idx is the index of the aimed row
   * @return the reference of the `row_idx`-th row
   */
  typename Matrix<T>::_row_ref_type operator[](const unsigned int row_idx)
  {
    return _row_ref_type(*this, row_idx);
  }

  /**
   * @brief Return the number of columns in the matrix
   *
   * @return the number of columns in the matrix
   */
  size_t num_of_cols() const
  {
    return _num_of_cols;
  }

  /**
   * @brief Return the number of rows in the matrix
   *
   * @return the number of rows in the matrix
   */
  size_t num_of_rows() const
  {
    return _num_of_rows;
  }

  /**
   * @brief Compare two matrices
   * 
   * @param A is a matrix
   * @return `true` if and only if this object and the parameter 
   *         have the same sizes and contains the same values
   */
  bool operator==(const Matrix<T>& A) const
  {
    if (A.num_of_rows()!=this->num_of_rows() || 
        A.num_of_cols()!=this->num_of_cols()) {
      return false;
    }

    auto A_row_it = std::begin(A._matrix);
    auto row_it = std::begin(this->_matrix);

    for (; A_row_it != std::end(A._matrix); ++A_row_it, ++row_it) {
      auto A_elem_it = std::begin(A_row_it->second);
      auto elem_it = std::begin(row_it->second);
      for (; A_elem_it != std::end(A_row_it->second); ++A_elem_it, ++elem_it) {
        if (A_elem_it->second != elem_it->second) {
          return false;
        }
      } 
    }

    return true;
  }

  /**
   * @brief Matrix-vector product
   *
   * @param v is the vector that must the multiplied to the matrix
   * @return the product between the current matrix and `v`
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
   * @param orig is the matrix that must be copied
   * @return a reference to the updated object
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
   * @param orig is the matrix that must be copied
   * @return a reference to the updated object
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
   * @return The row-column multiplication between this object and `A`
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

  /**
   * @brief Element-wise scalar sum for matrices
   * 
   * @tparam E is the type of scalar values
   * @param A is a sparse matrix
   * @param B is a sparse matrix
   * @return the sparse matrix \f$A+B\f$
   */
  template<typename E>
  friend Matrix<E> operator+(Matrix<E>&& A, const Matrix<E>& B)
  {
    for (auto row_it = std::begin(B._matrix); row_it != std::end(B._matrix);
          ++row_it) {
      for (auto elem_it = std::begin(row_it->second);
            elem_it != std::end(row_it->second); ++elem_it) {
        A[row_it->first][elem_it->first] += elem_it->second;
      }
    }

    return std::move(A);
  }

  /**
   * @brief Element-wise scalar subtraction for matrices
   * 
   * @tparam E is the type of scalar values
   * @param A is a sparse matrix
   * @param B is a sparse matrix
   * @return the sparse matrix \f$A-B\f$
   */
  template<typename E>
  friend Matrix<E> operator-(Matrix<E>&& A, const Matrix<E>& B)
  {
    for (auto row_it = std::begin(B._matrix); row_it != std::end(B._matrix);
          ++row_it) {
      for (auto elem_it = std::begin(row_it->second);
            elem_it != std::end(row_it->second); ++elem_it) {
        A[row_it->first][elem_it->first] -= elem_it->second;
      }
    }

    return std::move(A);
  }

  /**
   * @brief Element-wise scalar subtraction for matrices
   * 
   * @tparam E is the type of scalar values
   * @param A is a sparse matrix
   * @param B is a sparse matrix
   * @return the sparse matrix \f$A-B\f$
   */
  template<typename E>
  friend Matrix<E> operator-(const Matrix<E>& A, Matrix<E>&& B)
  {
    for (auto row_it = std::begin(B._matrix); row_it != std::end(B._matrix);
          ++row_it) {
      for (auto elem_it = std::begin(row_it->second);
            elem_it != std::end(row_it->second); ++elem_it) {
        elem_it->second *= -1;
      }
    }

    return std::move(B)+A;
  }

  /**
   * @brief Element-wise scalar multiplication for matrices
   * 
   * @tparam E is the type of scalar values
   * @param A is a matrix
   * @param scalar is a scalar value
   * @return the matrix \f$A*\textrm{scalar}\f$
   */
  template<typename E>
  friend Matrix<E> operator*(Matrix<E>&& A, const T scalar)
  {
    for (auto row_it = std::begin(A._matrix); row_it != std::end(A._matrix);
          ++row_it) {
      for (auto elem_it = std::begin(row_it->second);
            elem_it != std::end(row_it->second); ++elem_it) {
        elem_it->second *= scalar;
      }
    }

    return std::move(A);
  }

  /**
   * @brief Element-wise scalar division for matrices
   * 
   * @tparam E is the type of scalar values
   * @param A is a matrix
   * @param scalar is a scalar value
   * @return the matrix \f$A/\textrm{scalar}\f$
   */
  template<typename E>
  friend Matrix<E> operator/(Matrix<E>&& A, const E scalar)
  {
    if (scalar==0) {
      throw std::domain_error("Division by 0");
    }

    for (auto row_it = std::begin(A._matrix); row_it != std::end(A._matrix);
          ++row_it) {
      for (auto elem_it = std::begin(row_it->second);
            elem_it != std::end(row_it->second); ++elem_it) {
        elem_it->second /= scalar;
      }
    }

    return std::move(A);
  }

  template<typename E>
  friend class LUP_Factorization;

  /**
   * @brief Get the transpose matrix
   * 
   * @tparam E is the scalar value type
   * @param A is a sparse matrix
   * @return the transpose matrix \f$A^T\f$
   */
  template<typename E>
  friend Matrix<E> transpose(const Matrix<E>& A);

  /**
   * @brief Print a sparse matrix in a stream
   * 
   * @tparam E is the type of the scalar values
   * @param out is the output stream
   * @param A is the matrix to be printed
   * @return a reference to the output stream
   */
  template<typename E>
  friend std::ostream &std::operator<<(std::ostream &out, 
                                       const Matrix<E> &A);
};

/**
 * @brief Transpose a matrix
 * 
 * @tparam T is the scalar value type
 * @param A is the matrix to be transposed
 * @return the transposed matrix
 */
template<typename T>
Matrix<T> transpose(const Matrix<T>& A)
{ 
  Matrix<T> TA(A.num_of_cols(), A.num_of_rows());

  for (auto row_it = std::cbegin(A._matrix); row_it != std::cend(A._matrix);
         ++row_it) {
    for (auto elem_it = std::cbegin(row_it->second);
               elem_it != std::cend(row_it->second); ++elem_it) {
      TA[elem_it->first][row_it->first] = elem_it->second;
    }
  }

  return TA;
}

/**
 * @brief A LUP factorization for sparse matrices
 *
 * This class represents a Permutation, Lower-triangular,
 * Upper-triangular factorization for sparse matrices. It
 * also offers a method to solve linear systems
 * represented by sparse matrices. The permutation is not
 * represented by a matrix, but instead by using a
 * vector of indexes swaps.
 *
 * @tparam T any numeric type
 */
template<typename T>
class LUP_Factorization
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
   * @param M is a sparse matrix
   * @return the vector of non-zero rows below the diagonal
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
   * @brief Swap two matrix rows
   *
   * This method swaps the `row_idx`-th row and the row
   * below the diagonal whose value in column `row_idx` is not
   * zero and is the nearest to \f$1\f$.
   *
   * @param non_zero_below_diag is the vector such that
   *     `non_zero_below_diag[i]` is the set of row indexes `j`
   *     greater than `i` and having non-zero value in column `i`.
   * @param row_idx is the index of the row that must be swaped down
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

    // get the index of the row having a non-empty value in this
    // column which is the nearest to \f$1\f$.
    auto l = non_zero_below_diag[row_idx].begin();
    auto new_row_idx = *l;
    T delta = std::abs(1-std::abs(_U._matrix[new_row_idx][row_idx]));

    for (++l; l != non_zero_below_diag[row_idx].end(); ++l) {
      T new_delta = std::abs(1-std::abs(_U._matrix[*l][row_idx]));
      if (new_delta < delta) {
        new_row_idx = *l;
        delta = new_delta;
      }
    }
  
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
   * @brief Solve the linear system \f$L\cdot x = \f$
   *
   * @param[in] b is the known term of the linear system
   * @return the vector \f$x = L^{-1} \cdot \f$
   */
  std::vector<T> solveL(const std::vector<T> &b) const
  {
    std::vector<T> x(b.size());

    unsigned int i=0;
    for (auto row_it = std::begin(_L._matrix); row_it != std::end(_L._matrix);
         ++row_it, ++i) {
     if (std::begin(row_it->second)==std::end(row_it->second) ||   
            row_it->first != i) {
          throw std::domain_error("The linear system is underdetermined.");
      }
      T value = b[row_it->first];
      for (auto elem_it = ++std::rbegin(row_it->second);
           elem_it != std::rend(row_it->second); ++elem_it) {
        value -= elem_it->second * x[elem_it->first];
      }
      const T &diag_value = std::rbegin(row_it->second)->second;
      x[row_it->first] = value / diag_value;
    }

    if (i != _L.num_of_cols()) {
      throw std::domain_error("The linear system is underdetermined.");
    }

    return x;
  }

  /**
   * @brief Solve the linear system \f$U\cdot x = \f$
   *
   * @param[in] b is the known term of the linear system
   * @return the vector \f$x = U^{-1} \cdot \f$
   */
  std::vector<T> solveU(const std::vector<T> &b) const
  {
    std::vector<T> x(b.size());

    unsigned int i=_U.num_of_rows();
    for (auto row_it = std::rbegin(_U._matrix);
         row_it != std::rend(_U._matrix); ++row_it, --i) {

      if (std::begin(row_it->second)==std::end(row_it->second) ||   
            row_it->first != i-1) {
          throw std::domain_error("The linear system is underdetermined.");
      }
      T value = b[row_it->first];
      for (auto elem_it = ++std::begin(row_it->second);
           elem_it != std::end(row_it->second); ++elem_it) {
        value -= elem_it->second * x[elem_it->first];
      }
      const T &diag_value = std::begin(row_it->second)->second;
      x[row_it->first] = value / diag_value;
    }

    if (i != 0) {
      throw std::domain_error("The linear system is underdetermined.");
    }

    return x;
  }

public:
  /**
   * @brief Create a new LUP_Factorization object
   */
  LUP_Factorization(): _P(), _L(), _U() {}

  /**
   * @brief Copy construct a new LUP_Factorization object
   *
   * @param orig is the model for the new object
   */
  LUP_Factorization(const LUP_Factorization &orig):
      _P(orig._P), _L(orig._L), _U(orig._U)
  {
  }

  /**
   * @brief Create a new LUP_Factorization object
   *
   * This constructor performs the LUP factorization of
   * a sparse matrix.
   *
   * @param M is the matrix whose LUP factorization must be computed
   */
  LUP_Factorization(const Matrix<T> &M):
      _P(), _L(M.num_of_rows(), M.num_of_rows()), _U(M)
  {
    if (M.num_of_rows() == 0 || M.num_of_cols() == 0) {
      throw std::domain_error("The parameter is an empty matrix.");
    }

    std::vector<std::set<unsigned int>> non_zero_below_diag
        = get_non_zero_below_diag(M);

    for (auto row_it = std::begin(_U._matrix); row_it != std::end(_U._matrix);
         ++row_it) {
      const unsigned int &row_idx = row_it->first;

      if (row_idx == M.num_of_cols()) {
        return;
      }


      typename Matrix<T>::RowType::const_iterator diag_elem;
      if (!non_zero_below_diag[row_idx].empty()) {

        // swap the current row with the row below whose new 
        // diagonal element is the nearest to 1 different from 0
        swap_with_the_leastest_non_zero_row_in_column(non_zero_below_diag,
                                                      row_idx);

        diag_elem = _U._matrix[row_it->first].find(row_idx);
      } else {

        // get the diagonal element on the row
        diag_elem = row_it->second.find(row_idx);
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
              if (elem_it->first < *nz_row_it) {
                non_zero_below_diag[elem_it->first].insert(*nz_row_it);
              }
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

    for (unsigned int row_idx = 0; row_idx < _L.num_of_rows(); ++row_idx) {
      // set the L diagonal to 1
      _L._matrix[row_idx][row_idx] = 1;
    }
  }

  /**
   * @brief Compute the solution of \f$(L\cdot U)\cdot x=P^{-1}\cdot b\f$
   *
   * @param b is the known term of the linear system
   * @return the vector \f$x = (U^{-1} \cdot L^{-1} \cdot P^{-1}) \cdot b\f$
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
   * @brief Copy a LUP_Factorization object
   *
   * @param orig is the model that must be copied
   * @return a reference to the updated object
   */
  LUP_Factorization<T> &operator=(const LUP_Factorization<T> &orig)
  {
    _P = orig._P;
    _L = orig._L;
    _U = orig._U;

    return *this;
  }

  /**
   * @brief Copy a LUP_Factorization object
   *
   * @param orig is the model that must be copied
   * @return a reference to the updated object
   */
  LUP_Factorization<T> &operator=(LUP_Factorization<T> &&orig)
  {
    std::swap(_P, orig._P);
    std::swap(_L, orig._L);
    std::swap(_U, orig._U);

    return *this;
  }
};

/**
 * @brief Compute the rank of a sparse matrix
 *
 * @tparam T is the type of the matrix elements
 * @param A is the matrix whose rank must be computed
 * @return The rank of the matrix `A`
 */
template<typename T>
unsigned int rank(const Matrix<T> &A)
{
  const Matrix<T> U = LUP_Factorization<T>(A);

  unsigned int rank = 0;

  const unsigned max_diag = std::min(A.num_of_cols(), A.num_of_rows());
  for (unsigned int i = 0; i < max_diag; ++i) {
    if (U[i][i] != 0) {
      rank++;
    }
  }

  return rank;
}

/**
 * @brief Compute the determinant of a matrix
 * 
 * The determinant is evaluated by using three steps: first
 * of all, the LUP decomposition of the matrix is computed,
 * then the product of the elements in \f$U\f$ main diagonal
 * is computed, and, finally, this product is multiplied by 
 * \f$-1\f$ if the number of swaps in \f$P\f$ are odd.
 * 
 * @tparam T is the scalar domain type
 * @param A is the matrix whose determinant must be computed
 * @return the determinant of the matrix `A`
 */
template<typename T>
T determinant(const Matrix<T> &A)
{
  if (A.num_of_rows()==0 || A.num_of_cols()==0) {
    throw std::domain_error("Determinant of a 0x0-matrix not supported");
  }
  
  // Compute the factorization
  LUP_Factorization<T> fact(A);

  T det = 1;
  // multiply the elements in U main diagonal
  unsigned int elem_in_diag = std::min(A.num_of_rows(), 
                                       A.num_of_cols());
  for (unsigned int i=0; i<elem_in_diag; ++i)
  {
    det *= fact.U()[i][i];

    if (det == 0) {
      return 0;
    }
  }

  // Evaluate the number of row swaps and determine
  // the {1, -1} multiplier
  auto P = fact.P();
  std::vector<bool> swaps(A.num_of_rows(), true);
  for (auto P_it = P.begin(); P_it != P.end(); ++P_it) {
    unsigned int j=P_it->first;
    while (swaps[j]) {
      swaps[j] = false;
      j = P[j];
      if (j != P_it->first) {
        det = -1*det;
      }
    }
  }

  return det;
}

/**
 * @brief Compute the inverse matrix
 * 
 * @tparam T is the type of the scalar values
 * @param A is the matrix to be inverted
 * @return a matrix \f$A^{-1}\f$ such that \f$A\cdot A^{-1} = I\f$
 */
template<typename T>
Matrix<T> inverse(Matrix<T> &A)
{
  if (A.num_of_rows()==0 || A.num_of_cols()==0) {
    throw std::domain_error("Inverse of a 0x0-matrix not supported");
  }
  
  Matrix<T> Ainv;

  // Compute the factorization
  LUP_Factorization<T> fact(A);
  std::vector<T> vi(A.num_of_rows(),0);

  for (unsigned int i=0; i<A.num_of_rows(); ++i) {
    vi[i] = 1;
    Ainv.add_row(fact.solve(vi));
    vi[i] = 0;
  }

  return transpose(Ainv);
}

}
/*! @} End of SparseLinearAlgebra group */

/**
 * @brief Element-wise scalar sum for sparse matrices
 * 
 * @tparam T is the type of scalar values
 * @param A is a sparse matrix
 * @param B is a sparse matrix
 * @return the sparse matrix \f$A+B\f$
 */
template<typename T>
Sparse::Matrix<T> operator+(const Sparse::Matrix<T>& A, 
                            const Sparse::Matrix<T>&& B)
{
  return std::move(B)+A;
}

/**
 * @brief Element-wise scalar sum for sparse matrices
 * 
 * @tparam T is the type of scalar values
 * @param A is a sparse matrix
 * @param B is a sparse matrix
 * @return the sparse matrix \f$A+B\f$
 */
template<typename T>
Sparse::Matrix<T> operator+(const Sparse::Matrix<T>& A, 
                                         const Sparse::Matrix<T>& B)
{
  Sparse::Matrix<T> C(B);

  return std::move(C)+A;
}

/**
 * @brief Element-wise scalar subtraction for sparse matrices
 * 
 * @tparam T is the type of scalar values
 * @param A is a sparse matrix
 * @param B is a sparse matrix
 * @return the sparse matrix \f$A-B\f$
 */
template<typename T>
Sparse::Matrix<T> operator-(const Sparse::Matrix<T>& A, 
                                         const Sparse::Matrix<T>& B)
{
  Sparse::Matrix<T> C(A);

  return std::move(C)-B;
}

/**
 * @brief Element-wise scalar multiplication for matrices
 * 
 * @tparam T is the type of scalar values
 * @param A is a matrix
 * @param scalar is a scalar value
 * @return the matrix \f$A*\textrm{scalar}\f$
 */
template<typename T>
Sparse::Matrix<T> operator*(const Sparse::Matrix<T>& A, const T scalar)
{
  return operator*(Sparse::Matrix<T>(A), scalar); 
}

/**
 * @brief Element-wise scalar multiplication for matrices
 * 
 * @tparam T is the type of scalar values
 * @param scalar is a scalar value
 * @param A is a matrix
 * @return the matrix \f$\textrm{scalar}*A\f$
 */
template<typename T>
Sparse::Matrix<T> operator*(const T scalar, const Sparse::Matrix<T>& A)
{
  return A*scalar;
}

/**
 * @brief Element-wise scalar multiplication for matrices
 * 
 * @tparam T is the type of scalar values
 * @param scalar is a scalar value
 * @param A is a matrix
 * @return the matrix \f$\textrm{scalar}*A\f$
 */
template<typename T>
Sparse::Matrix<T> operator*(const T scalar, Sparse::Matrix<T>&& A)
{
  return std::move(operator*(std::move(A),scalar));
}

/**
 * @brief Element-wise scalar division for matrices
 * 
 * @tparam T is the type of scalar values
 * @param A is a matrix
 * @param scalar is a scalar value
 * @return the matrix \f$A/\textrm{scalar}\f$
 */
template<typename T>
Sparse::Matrix<T> operator/(const Sparse::Matrix<T>& A, const T scalar)
{
  return operator/(Sparse::Matrix<T>(A), scalar); 
}
}

#endif // LINEAR_ALGEBRA_H
