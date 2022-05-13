/**
 * @file LinearAlgebraIO.h
 * @author Alberto Casagrande <acasagrande@units.it>
 * @brief Contains linear algebra print functions
 * @version 0.1
 * @date 2021-11-26
 * 
 * @copyright Copyright (c) 2021-2022
 */

#ifndef LINEAR_ALGEBRA_IO_H
#define LINEAR_ALGEBRA_IO_H

#include <iostream>

#include "LinearAlgebra.h"

namespace std
{

/**
 * @brief Print a vector in a stream
 * 
 * @tparam T is the type of the scalar values
 * @param out is the output stream
 * @param v is the vector to be printed
 * @return a reference to the output stream
 */
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

/**
 * @brief Print a set in a stream
 * 
 * @tparam T is the type of the set values
 * @param out is the output stream
 * @param s is the set to be printed
 * @return a reference to the output stream
 */
template<typename T>
std::ostream &operator<<(std::ostream &out, const std::set<T> &s)
{
  out << "{";
  for (auto el_it = std::begin(s); el_it != std::end(s); ++el_it) {
    if (el_it != std::begin(s)) {
      out << ",";
    }
    out << *el_it;
  }
  out << "}";

  return out;
}

/**
 * @brief Print a dense matrix in a stream
 * 
 * @tparam T is the type of the scalar values
 * @param out is the output stream
 * @param A is the matrix to be printed
 * @return a reference to the output stream
 */
template<typename T>
std::ostream &operator<<(std::ostream &out,
                         const LinearAlgebra::Dense::Matrix<T> &A)
{
  for (auto row_it = std::begin(A); row_it != std::end(A); ++row_it) {
    if (row_it != std::begin(A)) {
      out << std::endl;
    }

    out << *row_it;
  }

  return out;
}

/**
 * @brief Print a sparse matrix in a stream
 * 
 * @tparam T is the type of the scalar values
 * @param out is the output stream
 * @param A is the matrix to be printed
 * @return a reference to the output stream
 */
template<typename T>
std::ostream &operator<<(std::ostream &out,
                         const LinearAlgebra::Sparse::Matrix<T> &A)
{
  out << "[ #rows: " << A.num_of_rows() << " #cols: " << A.num_of_cols() << " ";
  for (auto row_it = std::begin(A._matrix); row_it != std::end(A._matrix);
       ++row_it) {
    out << std::endl << " row " << row_it->first << ": [";
    for (auto elem_it = std::begin(row_it->second);
         elem_it != std::end(row_it->second); ++elem_it) {
      if (elem_it != std::begin(row_it->second)) {
        out << ", ";
      }
      out << "col " << elem_it->first << ": " << elem_it->second;
    }
    out << "]";
  }
  out << "]";

  return out;
}

/**
 * @brief Print a sparse LUP factorization in a stream
 * 
 * @tparam T is the type of the scalar values
 * @param out is the output stream
 * @param F is the LUP factorization to be printed
 * @return a reference to the output stream
 */
template<typename T>
std::ostream &operator<<(std::ostream &out,
                         const LinearAlgebra::Sparse::LUP_Factorization<T> &F)
{
  using namespace std;
  out << "{P=" << F.P() << "," << std::endl
      << " L=" << F.L() << "," << std::endl
      << " U=" << F.U() << "}";

  return out;
}

/**
 * @brief Print a dense LUP factorization in a stream
 * 
 * @tparam T is the type of the scalar values
 * @param out is the output stream
 * @param F is the LUP factorization to be printed
 * @return a reference to the output stream
 */
template<typename T>
std::ostream &operator<<(std::ostream &out,
                         const LinearAlgebra::Dense::LUP_Factorization<T> &F)
{
  using namespace std;
  out << "{P=" << F._P << "," << std::endl << " LU=" << F._LU << "}";

  return out;
}

}
#endif // LINEAR_ALGEBRA_H
