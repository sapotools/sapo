/**
 * @file Parallelotope.h
 * @author Tommaso Dreossi (tommasodreossi@berkeley.edu)
 * @author Alberto Casagrande <acasagrande@units.it>
 * @brief Represent and manipulate a parallelotope
 * @version 0.1
 * @date 2015-10-14
 *
 * @copyright Copyright (c) 2015-2022
 */

#ifndef PARALLELOTOPE_H_
#define PARALLELOTOPE_H_

#include <vector>

#include "Polytope.h"

/**
 * @brief A class to represent Parallelotope.
 *
 * Parallelotopes are internally represented as a base vertex,
 * a set of versors, and the lengths of the corresponding
 * tensors.
 *
 * @tparam T is the numeric type used to represent parallelotopes
 * @tparam APPROX_TYPE is the type of approximation used in place of T
 *      during the computation
 */
template<typename T = double, typename APPROX_TYPE = T>
class Parallelotope;

/**
 * @brief Swap two parallelotope
 *
 * This method swap two parallelotope objects.
 *
 * @tparam T is the numeric type used to represent parallelotopes
 * @tparam APPROX_TYPE is the type of approximation used in place of T
 *      during the computation
 * @param A is a parallelotope
 * @param B is a parallelotope
 */
template<typename T, typename APPROX_TYPE>
void swap(Parallelotope<T, APPROX_TYPE> &A, Parallelotope<T, APPROX_TYPE> &B);

template<typename T, typename APPROX_TYPE>
class Parallelotope
{
private:
  std::vector<LinearAlgebra::Vector<T>> _generators; //!< generators
  LinearAlgebra::Vector<T> _base_vertex;             //!< base vertex
  std::vector<T> _lengths; //!< tensors original lengths

public:
  /**
   * Constructor
   *
   * @param[in] directions is the vector of the parallelotope directions
   * @param[in] lower_bound is the lower bound offsets of the parallelotope
   * @param[in] upper_bound is the upper bound offsets of the parallelotope
   */
  Parallelotope(const std::vector<LinearAlgebra::Vector<T>> &directions,
                const LinearAlgebra::Vector<T> &lower_bound,
                const LinearAlgebra::Vector<T> &upper_bound);

  /**
   * Get the parallelotope dimension
   *
   * @returns parallelotope dimension
   */
  size_t dim() const;

  /**
   * Build a polytope representing the parallelotope.
   *
   * @returns a polytope representing the parallelotope
   */
  operator Polytope<T, APPROX_TYPE>() const;

  /**
   * @brief Get the base vertex
   *
   * @return a reference to the base vertex
   */
  const LinearAlgebra::Vector<T> &base_vertex() const;

  /**
   * @brief Get the parallelotope's center
   *
   * @return the parallelotope's center
   */
  LinearAlgebra::Vector<T> center() const;

  /**
   * @brief Get the generator lengths
   *
   * @return a reference to the generator lengths
   */
  const std::vector<T> &lengths() const;

  /**
   * @brief Get the parallelotope generators
   *
   * @return a reference to the parallelotope generators
   */
  const std::vector<LinearAlgebra::Vector<T>> &generators() const;

  friend void swap<T, APPROX_TYPE>(Parallelotope<T, APPROX_TYPE> &A,
                                   Parallelotope<T, APPROX_TYPE> &B);
};

/**
 * @brief Test whether two sets are the same
 *
 * @tparam T is the numeric type used to represent parallelotopes
 * @tparam APPROX_TYPE is the type of approximation used in place of T
 *      during the computation
 * @param A is a polytope
 * @param B is a parallelotope
 * @return `true` when the method can establish that `A` and `B`
 *         represent the same set. `false` when it can
 *         establish `A` and `B` differ. `uncertain` in the
 *         remaining cases
 */
template<typename T, typename APPROX_TYPE>
TriBool operator==(const Polytope<T, APPROX_TYPE> &A,
                   const Parallelotope<T, APPROX_TYPE> &B);

/**
 * @brief Test whether two sets are the same
 *
 * @tparam T is the numeric type used to represent parallelotopes
 * @tparam APPROX_TYPE is the type of approximation used in place of T
 *      during the computation
 * @param A is a parallelotope
 * @param B is a polytope
 * @return `true` when the method can establish that `A` and `B`
 *         represent the same set. `false` when it can
 *         establish `A` and `B` differ. `uncertain` in the
 *         remaining cases
 */
template<typename T, typename APPROX_TYPE>
TriBool operator==(const Parallelotope<T, APPROX_TYPE> &A,
                   const Polytope<T, APPROX_TYPE> &B);

/**
 * @brief Test whether two parallelotopes are the same
 *
 * @tparam T is the numeric type used to represent parallelotopes
 * @tparam APPROX_TYPE is the type of approximation used in place of T
 *      during the computation
 * @param A is a parallelotope
 * @param B is a parallelotope
 * @return `true` when the method can establish that `A` and `B`
 *         represent the same set. `false` when it can
 *         establish `A` and `B` differ. `uncertain` in the
 *         remaining cases
 */
template<typename T, typename APPROX_TYPE>
TriBool operator==(const Parallelotope<T, APPROX_TYPE> &A,
                   const Parallelotope<T, APPROX_TYPE> &B);

/**
 * @brief Test whether two sets differ
 *
 * @tparam T is the numeric type used to represent parallelotopes
 * @tparam APPROX_TYPE is the type of approximation used in place of T
 *      during the computation
 * @param A is a polytope
 * @param B is a parallelotope
 * @return `true` when the method can establish that `A` and `B`
 *         differ. `false` when it can establish `A` and `B`
 *         represent the same set. `uncertain` in the remaining
 *         cases
 */
template<typename T, typename APPROX_TYPE>
TriBool operator!=(const Polytope<T, APPROX_TYPE> &A,
                   const Parallelotope<T, APPROX_TYPE> &B);

/**
 * @brief Test whether two sets differ
 *
 * @tparam T is the numeric type used to represent parallelotopes
 * @tparam APPROX_TYPE is the type of approximation used in place of T
 *      during the computation
 * @param A is a parallelotope
 * @param B is a polytope
 * @return `true` when the method can establish that `A` and `B`
 *         differ. `false` when it can establish `A` and `B`
 *         represent the same set. `uncertain` in the remaining
 *         cases
 */
template<typename T, typename APPROX_TYPE>
TriBool operator!=(const Parallelotope<T, APPROX_TYPE> &A,
                   const Polytope<T, APPROX_TYPE> &B);

/**
 * @brief Test whether two parallelotopes differ
 *
 * @tparam T is the numeric type used to represent parallelotopes
 * @tparam APPROX_TYPE is the type of approximation used in place of T
 *      during the computation
 * @param A is a parallelotope
 * @param B is a parallelotope
 * @return `true` when the method can establish that `A` and `B`
 *         differ. `false` when it can establish `A` and `B`
 *         represent the same set. `uncertain` in the remaining
 *         cases
 */
template<typename T, typename APPROX_TYPE>
TriBool operator!=(const Parallelotope<T, APPROX_TYPE> &A,
                   const Parallelotope<T, APPROX_TYPE> &B);

template<typename T, typename APPROX_TYPE>
Parallelotope<T, APPROX_TYPE>::Parallelotope(
    const std::vector<LinearAlgebra::Vector<T>> &directions,
    const LinearAlgebra::Vector<T> &lower_bound,
    const LinearAlgebra::Vector<T> &upper_bound)
{
  using namespace LinearAlgebra;

  if (lower_bound.size() != upper_bound.size()) {
    SAPO_ERROR("the lower and upper bounds vectors must "
               "have the same number of dimensions",
               std::domain_error);
  }

  if (lower_bound.size() != directions.size()) {
    SAPO_ERROR("the lower bounds vector and the direction vector must "
               "have the same number of dimensions",
               std::domain_error);
  }

  Dense::LUP_Factorization<double> factorization(directions);

  try {
    // store the base vertex
    _base_vertex = factorization.solve(lower_bound);

  } catch (std::domain_error &) { // if a domain_error is raised, then the
                                  // template is singular
    SAPO_ERROR("singular parallelotope", std::runtime_error);
  }

  // Compute the versors
  std::vector<T> offset(upper_bound.size(), 0);
  for (unsigned int k = 0; k < upper_bound.size(); k++) {
    const T delta_k = upper_bound[k] - lower_bound[k];
    offset[k] = 1;

    const std::vector<T> tensor = factorization.solve(offset);
    offset[k] = 0;

    // compute its length and store it
    const T length = norm_2(tensor);
    if (delta_k != 0) {
      _lengths.push_back(delta_k * length);
    } else {
      _lengths.push_back(0);
    }

    // compute and store the corresponding versor
    _generators.push_back(tensor / length);
  }
}

template<typename T, typename APPROX_TYPE>
inline size_t Parallelotope<T, APPROX_TYPE>::dim() const
{
  return (_generators.size() == 0 ? 0 : _generators.front().size());
}

template<typename T, typename APPROX_TYPE>
Parallelotope<T, APPROX_TYPE>::operator Polytope<T, APPROX_TYPE>() const
{
  using namespace LinearAlgebra;

  std::vector<Vector<T>> A(2 * dim());
  Vector<T> b(2 * dim());
  for (unsigned int i = 0; i < dim(); ++i) {
    const T base_bound = _generators[i] * _base_vertex;
    const T vertex_bound = base_bound + _lengths[i];

    A[i] = Vector<T>(_generators[i]);
    A[i + dim()] = -_generators[i];
    b[i] = std::max(base_bound, vertex_bound);
    b[i + dim()] = -std::min(base_bound, vertex_bound);
  }

  return Polytope(A, b);
}

template<typename T, typename APPROX_TYPE>
inline const LinearAlgebra::Vector<T> &
Parallelotope<T, APPROX_TYPE>::base_vertex() const
{
  return this->_base_vertex;
}

template<typename T, typename APPROX_TYPE>
LinearAlgebra::Vector<T> Parallelotope<T, APPROX_TYPE>::center() const
{
  using namespace LinearAlgebra;

  Vector<T> c{base_vertex()};

  for (size_t i = 0; i < dim(); ++i) {
    c = c + _lengths[i] / 2 * _generators[i];
  }

  return c;
}

template<typename T, typename APPROX_TYPE>
inline const std::vector<T> &Parallelotope<T, APPROX_TYPE>::lengths() const
{
  return this->_lengths;
}

template<typename T, typename APPROX_TYPE>
inline const std::vector<LinearAlgebra::Vector<T>> &
Parallelotope<T, APPROX_TYPE>::generators() const
{
  return this->_generators;
}

template<typename T, typename APPROX_TYPE>
void swap(Parallelotope<T, APPROX_TYPE> &A, Parallelotope<T, APPROX_TYPE> &B)
{
  std::swap(A._generators, B._generators);
  std::swap(A._base_vertex, B._base_vertex);
  std::swap(A._lengths, B._lengths);
}

template<typename T, typename APPROX_TYPE>
inline TriBool operator==(const Polytope<T, APPROX_TYPE> &A,
                          const Parallelotope<T, APPROX_TYPE> &B)
{
  return static_cast<Polytope<T, APPROX_TYPE>>(B) == A;
}

template<typename T, typename APPROX_TYPE>
inline TriBool operator==(const Parallelotope<T, APPROX_TYPE> &A,
                          const Polytope<T, APPROX_TYPE> &B)
{
  return B == A;
}

template<typename T, typename APPROX_TYPE>
inline TriBool operator==(const Parallelotope<T, APPROX_TYPE> &A,
                          const Parallelotope<T, APPROX_TYPE> &B)
{
  return static_cast<Polytope<T, APPROX_TYPE>>(A)
         == static_cast<Polytope<T, APPROX_TYPE>>(B);
}

template<typename T, typename APPROX_TYPE>
inline TriBool operator!=(const Polytope<T, APPROX_TYPE> &A,
                          const Parallelotope<T, APPROX_TYPE> &B)
{
  return !(A == B);
}

template<typename T, typename APPROX_TYPE>
inline TriBool operator!=(const Parallelotope<T, APPROX_TYPE> &A,
                          const Polytope<T, APPROX_TYPE> &B)
{
  return !(A == B);
}

template<typename T, typename APPROX_TYPE>
inline TriBool operator!=(const Parallelotope<T, APPROX_TYPE> &A,
                          const Parallelotope<T, APPROX_TYPE> &B)
{
  return !(A == B);
}
#endif /* PARALLELOTOPE_H_ */
