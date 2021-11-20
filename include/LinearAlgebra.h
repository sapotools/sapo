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
std::vector<T> &operator-(std::vector<T> &&v)
{
  transform(v.begin(), v.end(), v.begin(), std::negate<T>());

  return v;
}

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
  std::vector<unsigned int> _P;    //!< the permutation
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
      for (unsigned int i = j + 1; i < row.size(); ++i) {
        value -= row[i] * x[i];
      }
      x[j] = value / row[j];
    }

    return x;
  }

public:
  /**
   * @brief Construct a new PLU_Factorization object
   *
   */
  PLU_Factorization(): _P(), _LU() {}

  /**
   * @brief Construct a new PLU_Factorization object
   *
   * @param orig is the model for the new object.
   */
  PLU_Factorization(const PLU_Factorization &orig): _P(orig._P), _LU(orig._LU)
  {
  }

  /**
   * @brief Construct a new PLU_Factorization object
   *
   * This constructor performs the P-LU factorization of
   * a dense matrix.
   *
   * @param M is the matrix whose P-LU factorization must be computed.
   */
  PLU_Factorization(const std::vector<std::vector<T>> &M): _P(M.size()), _LU(M)
  {
    if (M.size() == 0) {
      throw std::domain_error("The parameter is an empty matrix.");
    }

    if (M.size() != M.front().size()) {
      throw std::domain_error("The parameter is not a square matrix.");
    }

    // fill _P with the identity permutation
    std::iota(std::begin(_P), std::end(_P), 0);

    for (unsigned int j = 0; j < _LU.size(); ++j) {

      // find the first non-null value below A[j-1][j]
      unsigned int k = j;
      while (k < _LU.size() && _LU[k][j] == 0) {
        k++;
      }

      // if it does not exist, return an empty solution
      if (k >= _LU.size()) {
        throw std::domain_error("The matrix is singular.");
      }

      // otherwise, swap rows
      if (j != k) {
        std::swap(_LU[j], _LU[k]);
        std::swap(_P[j], _P[k]);
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

  /**
   * @brief Compute the solution of the system $(L \cot U) \cdot x = P^{-1}
   * \cdot b$
   *
   * @param is the known term of the linear system.
   * @return the vector $x = (U^{-1} \cdot L^{-1} \cdot P^{-1}) \cdot b$.
   */
  std::vector<T> solve(const std::vector<T> &b) const
  {
    if (b.size() != _P.size()) {
      throw std::domain_error("Wrong dimensions");
    }

    std::vector<T> Pb(b.size());

    // Apply the permutation to b
    for (unsigned int i = 0; i < b.size(); ++i) {
      Pb[i] = b[_P[i]];
    }

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
  /**
   * @brief A proxy for rows of sparse matrices
   */
  class _const_RowType
  {
    const std::map<unsigned int, T> *_row; //!< a pointer to the actual row
    size_t _num_of_cols;                   //!< the row size

    /**
     * @brief Construct a new const RowType object
     *
     * @param num_of_cols is the row size
     */
    _const_RowType(const size_t num_of_cols):
        _row(NULL), _num_of_cols(num_of_cols)
    {
    }

    /**
     * @brief Construct a new const RowType object
     *
     * @param row is a pointer to the actual row
     * @param num_of_cols is the row size
     */
    _const_RowType(const std::map<unsigned int, T> *row,
                   const size_t num_of_cols):
        _row(row),
        _num_of_cols(num_of_cols)
    {
    }

  public:
    /**
     * @brief Construct a new const RowType object
     *
     * @param orig is the model for the new object.
     */
    _const_RowType(const _const_RowType &orig):
        _row(orig._row), _num_of_cols(orig._num_of_cols)
    {
    }

    /**
     * @brief Return the value in the `i`-th position of the row
     *
     * @param i is the index of the aimed value.
     * @return a copy of the value in position `i`.
     */
    T operator[](const unsigned int i) const
    {
      if (i >= _num_of_cols) {
        throw std::domain_error("Out-of_bound: the column does not exist.");
      }

      if (_row == NULL) {
        return 0;
      }

      auto it = _row->find(i);
      if (it == std::end(*_row)) {
        return 0;
      }

      return it->second;
    }

    friend Matrix;
  };

  typedef std::map<unsigned int, T> RowType;

  std::map<unsigned int, RowType> _matrix; //!< indexed list of rows
  size_t _num_of_rows;                     //!< number of rows
  size_t _num_of_cols;                     //!< number of columns

public:
  /**
   * @brief Construct a new Matrix object
   *
   */
  Matrix(): _matrix(), _num_of_rows(0), _num_of_cols(0) {}

  /**
   * @brief Construct a new Matrix object
   *
   * @param num_of_rows is the number of rows in the new matrix.
   * @param num_of_cols is the number of columns in the new matrix.
   */
  Matrix(const unsigned int num_of_rows, const unsigned int num_of_cols):
      _matrix(), _num_of_rows(num_of_rows), _num_of_cols(num_of_cols)
  {
  }

  /**
   * @brief Construct a new Matrix object
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
   * @brief Construct a new Matrix object
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

    auto key_compare = [](const std::pair<unsigned int, T> &a,
                          const std::pair<unsigned int, T> &b) -> bool {
      return a.first < b.first;
    };
    auto greatest
        = std::max_element(std::begin(row), std::end(row), key_compare);

    if (_num_of_cols < greatest->first) {
      _num_of_cols = greatest->first;
    }

    _matrix[row_idx] = row;

    return *this;
  }

  /**
   * @brief Get a proxy of one of the matrix rows
   *
   * @param j is the index of the aimed row.
   * @return the proxy of the `j`-th row.
   */
  typename Matrix<T>::_const_RowType operator[](const unsigned int j) const
  {
    if (j >= _num_of_rows) {
      throw std::domain_error("Out-of_bound: the row does not exist");
    }

    auto it = _matrix.find(j);
    if (it == std::end(_matrix)) {
      return _const_RowType(_num_of_cols);
    }

    return _const_RowType(&it->second, _num_of_cols);
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
  std::vector<unsigned int> _P; //!< the permutation
  Matrix<T> _L;                 //!< the lower-triangular sparse matrix
  Matrix<T> _U;                 //!< the upper-triangular sparse matrix

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
    std::swap(_P[row_idx], _P[new_row_idx]);
    std::swap(_U._matrix[row_idx], _U._matrix[new_row_idx]);
    std::swap(_L._matrix[row_idx], _L._matrix[new_row_idx]);

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
      T value = b[row_it->first];
      for (auto elem_it = ++std::rbegin(row_it->second);
           elem_it != std::rend(row_it->second); ++elem_it) {
        value -= elem_it->second * x[elem_it->first];
      }
      x[row_it->first] = value / (std::rbegin(row_it->second)->second);
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
      T value = b[row_it->first];
      for (auto elem_it = ++std::begin(row_it->second);
           elem_it != std::end(row_it->second); ++elem_it) {
        value -= elem_it->second * x[elem_it->first];
      }
      x[row_it->first] = value / (std::begin(row_it->second)->second);
    }

    return x;
  }

public:
  /**
   * @brief Construct a new PLU_Factorization object
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
   * @brief Construct a new PLU_Factorization object
   *
   * This constructor performs the P-LU factorization of
   * a sparse matrix.
   *
   * @param M is the matrix whose P-LU factorization must be computed.
   */
  PLU_Factorization(const Matrix<T> &M):
      _P(M.num_of_rows()), _L(M.num_of_rows(), M.num_of_cols()), _U(M)
  {
    if (M.num_of_rows() == 0 || M.num_of_cols() == 0) {
      throw std::domain_error("The parameter is not an empty matrix.");
    }

    if (M.num_of_cols() != M.num_of_rows()) {
      throw std::domain_error("The parameter is not a square matrix.");
    }

    // fill _P with the identity permutation
    std::iota(std::begin(_P), std::end(_P), 0);

    std::vector<std::set<unsigned int>> non_zero_below_diag
        = get_non_zero_below_diag(M);

    for (auto row_it = std::begin(_U._matrix); row_it != std::end(_U._matrix);
         ++row_it) {
      const unsigned int &row_idx = row_it->first;

      // get the diagonal element on the row
      auto diag_elem = row_it->second.find(row_idx);

      if (diag_elem == std::end(row_it->second)) {

        // the diagonal element is empty
        swap_with_the_leastest_non_zero_row_in_column(non_zero_below_diag,
                                                      row_idx);

        diag_elem = _U._matrix[row_it->first].find(row_idx);
      }

      // the diagonal element is non-empty
      auto &U_row = _U._matrix[row_it->first];

      // set the L diagonal to 1
      _L._matrix[row_idx][row_idx] = 1;

      // for any non-zero row in this column
      for (auto nz_row_it = std::cbegin(non_zero_below_diag[row_idx]);
           nz_row_it != std::cend(non_zero_below_diag[row_idx]); ++nz_row_it) {

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

  /**
   * @brief Compute the solution of the system $(L \cot U) \cdot x = P^{-1}
   * \cdot b$
   *
   * @param is the known term of the linear system.
   * @return the vector $x = (U^{-1} \cdot L^{-1} \cdot P^{-1}) \cdot b$.
   */
  std::vector<T> solve(const std::vector<T> &v) const
  {
    if (v.size() != _P.size()) {
      throw std::domain_error("Wrong dimensions");
    }

    std::vector<T> b(v.size());

    // Apply the permutation to b
    for (unsigned int i = 0; i < v.size(); ++i) {
      b[i] = v[_P[i]];
    }

    // Forward substitution for L
    b = solveL(b);

    // Backward substitution for U
    return solveU(b);
  }

  /**
   * @brief Get the permutation
   *
   * @return the factorization permutation
   */
  const std::vector<unsigned int> &P() const
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
  os << "[ #rows: " << M.num_of_rows() << " #cols: " << M.num_of_cols();
  for (auto row_it = std::begin(M._matrix); row_it != std::end(M._matrix);
       ++row_it) {
    os << std::endl << " row " << row_it->first << ": [";
    for (auto elem_it = std::begin(row_it->second);
         elem_it != std::end(row_it->second); ++elem_it) {
      if (elem_it != std::begin(row_it->second)) {
        os << ", ";
      }
      os << "(" << elem_it->first << "," << elem_it->second << ")";
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

#endif // LINEAR_ALGEBRA_H
