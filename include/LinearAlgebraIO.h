/**
 * @file LinearAlgebra.h
 * Contains linear algebra classes and functions code
 *
 * @author Alberto Casagrande <acasagrande@units.it>
 * @version 0.1
 */

#ifndef LINEAR_ALGEBRA_IO_H
#define LINEAR_ALGEBRA_IO_H

#include <iostream>

#include "LinearAlgebra.h"

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
                         const SparseLinearAlgebra::LUP_Factorization<T> &D)
{
  using namespace std;
  os << "{P=" << D.P() << "," << std::endl
     << " L=" << D.L() << "," << std::endl
     << " U=" << D.U() << "}";

  return os;
}

template<typename T>
std::ostream &operator<<(std::ostream &os,
                         const DenseLinearAlgebra::LUP_Factorization<T> &D)
{
  using namespace std;
  os << "{P=" << D._P << "," << std::endl << " LU=" << D._LU << "}";

  return os;
}
}
#endif // LINEAR_ALGEBRA_H
