/**
 * @file Polytope.h
 * @author <acasagrande@units.it>
 * @brief  Represent and manipulate polytopes
 * @version 0.1
 * @date 2021-11-16
 *
 * @copyright Copyright (c) 2021-2022
 */

#ifndef POLYTOPE_H_
#define POLYTOPE_H_

#include "LinearSystem.h"

/**
 * @brief A polytope representation class
 *
 * Polytopes are bounded convex sets representable
 * by a linear system \f[A \cdot x \leq b\f].
 * This class extends the `LinearSystem` class and
 * provides some set-based methods and functions
 * such as `is_empty` or `intersect`.
 *
 * @tparam T is the numeric type used to represent polytopes
 * @tparam APPROX_TYPE is the type of approximation used in place of T
 *      during the computation
 */
template<typename T = double, typename APPROX_TYPE = T>
class Polytope : public LinearSystem<T, APPROX_TYPE>
{

private:
  /// @private
  std::list<Polytope<T, APPROX_TYPE>>
  split(const std::vector<bool> &bvect_base, const unsigned int cidx,
        std::list<Polytope<T, APPROX_TYPE>> &tmp_covering,
        std::vector<std::vector<T>> &A, std::vector<T> &b,
        const unsigned int num_of_splits) const;

public:
  /**
   * Constructor that instantiates a unbounded polytope
   */
  Polytope();

  /**
   * Copy constructor
   *
   * @tparam APPROX2 is the approximation type of the
   *         original polytope
   * @param[in] orig the original polytope
   */
  template<typename APPROX2 = APPROX_TYPE>
  Polytope(const Polytope<T, APPROX2> &orig);

  /**
   * Move constructor
   *
   * @tparam APPROX2 is the approximation type of the
   *         original polytope
   * @param[in] orig the original polytope
   */
  template<typename APPROX2 = APPROX_TYPE>
  Polytope(Polytope<T, APPROX2> &&orig);

  /**
   * Constructor
   *
   * @param[in] A template matrix
   * @param[in] b offset vector
   * @todo Test polytope closeness
   */
  Polytope(const std::vector<std::vector<T>> &A, const std::vector<T> &b);

  /**
   * Move Constructor
   *
   * @param[in] A template matrix
   * @param[in] b offset vector
   */
  Polytope(std::vector<std::vector<T>> &&A, std::vector<T> &&b);

  /**
   * Constructor from a set of symbolic expressions
   *
   * @param[in] vars list of variables appearing in the constraints
   * @param[in] constraints symbolic constraints
   * @todo Test polytope closeness
   */
  Polytope(const std::vector<SymbolicAlgebra::Symbol<T>> &vars,
           const std::vector<SymbolicAlgebra::Expression<T>> &constraints);

  /**
   * @brief Copy operator
   *
   * @tparam APPROX2 is the approximation type of the
   *         original polytope
   * @param orig is the original Polytope
   * @return a reference to the updated object
   */
  template<typename APPROX2>
  Polytope<T, APPROX_TYPE> &operator=(const Polytope<T, APPROX2> &orig);

  /**
   * @brief Copy operator
   *
   * @tparam APPROX2 is the approximation type of the
   *         original polytope
   * @param orig is the original Polytope
   * @return a reference to the updated object
   */
  template<typename APPROX2>
  Polytope<T, APPROX_TYPE> &operator=(Polytope<T, APPROX2> &&orig);

  /**
   * Establish whether a polytope is empty
   *
   * @param[in] strict_inequality specifies whether the polytope is
   *         defined by a strict inequality (i.e., Ax < b).
   * @return `true` when it can establish that the polytope is
   *         empty. `false` when it can establish that the polytope
   *         is not empty. `uncertain` in the other cases
   */
  TriBool is_empty(const bool strict_inequality = false) const;

  /**
   * Establish whether a polytope interior is empty
   *
   * @return `true` when it can establish that the polytope interior
   *         is empty. `false` when it can establish that the polytope
   *         interior is not empty. `uncertain` in the other cases
   */
  TriBool is_interior_empty() const;

  /**
   * Check whether one polytope is a subset of another polytope
   *
   * This method establishes whether the current polytope is a
   * subset of another polytope.
   *
   * @param[in] P is the polytope that are compared to the current
   *             object
   * @return `true` if the method can establish that this polytope
   *         is a subset of `P`. `false` when it can establish that
   *         this polytope is not a subset of `P`. `uncertain` in
   *         the remaining cases
   */
  TriBool is_subset_of(const Polytope<T, APPROX_TYPE> &P) const;

  /**
   * Check whether one polytope is a superset of another polytope
   *
   * This method establishes whether the current polytope is a
   * superset of another polytope.
   *
   * @param[in] P is the polytope that are compared to the current
   *             object
   * @return `true` if the method can establish that this polytope
   *         is a superset of `P`. `false` when it can establish that
   *         this polytope is not a superset of `P`. `uncertain` in
   *         the remaining cases
   */
  TriBool includes(const Polytope<T, APPROX_TYPE> &P) const;

  /**
   * @brief Expand the polytope
   *
   * This method expands the polytope so that each of its boundaries
   * is moved by a value `delta`.
   *
   * @param delta is the aimed expansion
   * @return a reference to the updated polytope
   */
  Polytope<T, APPROX_TYPE> &expand_by(const T delta);

  /**
   *  Split a polytope in a list of polytopes.
   *
   *  This method splits a polytope in a list of polytopes such
   *  that their set union equals the original polytope.
   *
   *  @return A list of polytopes such that their union equals
   *      the original polytope.
   */
  std::list<Polytope<T, APPROX_TYPE>> split() const;

  /**
   *  Split a polytope in a list of polytopes.
   *
   *  This method splits a polytope in a list of polytopes such
   *  that their set union equals the original polytope. Each
   *  polytope in the original list is split in
   *  \f$2^\textrm{num_of_splits}\f$ polytopes at most.
   *
   * @param num_of_splits is the number of splits to performed.
   * @return A list of polytopes such that their union equals
   *      the original polytope.
   */
  std::list<Polytope<T, APPROX_TYPE>>
  split(const unsigned int num_of_splits) const;

  /**
   * Update a polytope by intersecting it with another one.
   *
   * This method works in-place and changes the calling object.
   *
   * @param[in] P is a polytope.
   * @return a reference to the updated object.
   */
  Polytope<T, APPROX_TYPE> &intersect_with(const Polytope<T, APPROX_TYPE> &P);

  /**
   * Determine the volume of the bounding box of the polytope
   *
   * @return volume of the bounding box
   */
  T bounding_box_volume() const;

  template<typename T2, typename APPROX2>
  friend void swap(Polytope<T2, APPROX2> &P1, Polytope<T2, APPROX2> &P2);

  template<typename T2, typename APPROX2>
  friend Polytope<T2, APPROX2> intersect(const Polytope<T2, APPROX2> &P1,
                                         const Polytope<T2, APPROX2> &P2);

  template<typename T2, typename APPROX2>
  friend Polytope<T2, APPROX2>
  over_approximate_union(const Polytope<T2, APPROX2> &P1,
                         const Polytope<T2, APPROX2> &P2);
};

/**
 * @brief Test whether two polytopes are disjoint
 *
 * @tparam T is the numeric type used to represent polytopes
 * @tparam APPROX_TYPE is the type of approximation used in place of T
 *      during the computation
 * @param P1 is a polytope
 * @param P2 is a polytope
 * @return `true` if and only if `P1` and `P2` are disjoint
 * @return `true` if the method can establish that `P1` and `P2`
 *         are disjoint. `false` when it can establish that
 *         `P1` and `P2` are not disjoint. `uncertain` in
 *         the remaining cases
 */
template<typename T, typename APPROX_TYPE>
TriBool are_disjoint(const Polytope<T, APPROX_TYPE> &P1,
                     const Polytope<T, APPROX_TYPE> &P2);

/**
 * @brief Get the expansion of a linear set
 *
 * This method expands a linear set so that each of its boundaries
 * is moved by a value `epsilon`.
 *
 * @tparam BASIC_SET_TYPE is the basic set type
 * @tparam T is the numeric type
 * @param S is the set to be expanded
 * @param epsilon is the aimed expansion
 * @return an expanded version of `S`
 */
template<class BASIC_SET_TYPE, typename T>
BASIC_SET_TYPE expand(const BASIC_SET_TYPE &S, const T epsilon);

/**
 * Compute and over-approximation of the union of two polytopes
 *
 * @tparam T is the numeric type used to represent polytopes
 * @tparam APPROX_TYPE is the type of approximation used in place of T
 *      during the computation
 * @param[in] P1 is a polytope
 * @param[in] P2 is a polytope
 * @return an over-approximation of the union of `P1` and `P2`
 */
template<typename T, typename APPROX_TYPE>
Polytope<T, APPROX_TYPE>
over_approximate_union(const Polytope<T, APPROX_TYPE> &P1,
                       const Polytope<T, APPROX_TYPE> &P2);

template<typename T, typename APPROX_TYPE>
Polytope<T, APPROX_TYPE>::Polytope(): LinearSystem<T, APPROX_TYPE>()
{
}

template<typename T, typename APPROX_TYPE>
template<typename APPROX2>
Polytope<T, APPROX_TYPE>::Polytope(const Polytope<T, APPROX2> &orig):
    LinearSystem<T, APPROX_TYPE>(orig)
{
}

template<typename T, typename APPROX_TYPE>
template<typename APPROX2>
Polytope<T, APPROX_TYPE>::Polytope(Polytope<T, APPROX2> &&orig):
    LinearSystem<T, APPROX_TYPE>(std::move(orig))
{
}

template<typename T, typename APPROX_TYPE>
Polytope<T, APPROX_TYPE>::Polytope(const std::vector<std::vector<T>> &A,
                                   const std::vector<T> &b):
    LinearSystem<T, APPROX_TYPE>(A, b)
{
}

template<typename T, typename APPROX_TYPE>
Polytope<T, APPROX_TYPE>::Polytope(std::vector<std::vector<T>> &&A,
                                   std::vector<T> &&b):
    LinearSystem<T, APPROX_TYPE>(std::move(A), std::move(b))
{
}

template<typename T, typename APPROX_TYPE>
Polytope<T, APPROX_TYPE>::Polytope(
    const std::vector<SymbolicAlgebra::Symbol<T>> &vars,
    const std::vector<SymbolicAlgebra::Expression<T>> &constraints):
    LinearSystem<T, APPROX_TYPE>(vars, constraints)
{
}

template<typename T, typename APPROX_TYPE>
template<typename APPROX2>
Polytope<T, APPROX_TYPE> &
Polytope<T, APPROX_TYPE>::operator=(const Polytope<T, APPROX2> &orig)
{
  static_cast<LinearSystem<T, APPROX_TYPE> *>(this)->operator=(orig);

  return *this;
}

template<typename T, typename APPROX_TYPE>
template<typename APPROX2>
Polytope<T, APPROX_TYPE> &
Polytope<T, APPROX_TYPE>::operator=(Polytope<T, APPROX2> &&orig)
{
  static_cast<LinearSystem<T, APPROX_TYPE> *>(this)->operator=(
      std::move(orig));

  return *this;
}

template<typename T, typename APPROX_TYPE>
inline TriBool
Polytope<T, APPROX_TYPE>::is_empty(const bool strict_inequality) const
{
  return !this->has_solutions(strict_inequality);
}

template<typename T, typename APPROX_TYPE>
inline TriBool Polytope<T, APPROX_TYPE>::is_interior_empty() const
{
  return this->is_empty(true);
}

template<typename T, typename APPROX_TYPE>
inline TriBool
Polytope<T, APPROX_TYPE>::is_subset_of(const Polytope<T, APPROX_TYPE> &P) const
{
  return this->satisfies(P);
}

template<typename T, typename APPROX_TYPE>
inline TriBool
Polytope<T, APPROX_TYPE>::includes(const Polytope<T, APPROX_TYPE> &P) const
{
  return P.satisfies(*this);
}

template<typename T, typename APPROX_TYPE>
Polytope<T, APPROX_TYPE> &Polytope<T, APPROX_TYPE>::expand_by(const T delta)
{
  for (size_t i = 0; i < this->_b.size(); ++i) {
    this->_b[i] += delta * LinearAlgebra::norm_2(this->_A[i]);
  }

  return *this;
}

template<typename T, typename APPROX_TYPE>
inline std::list<Polytope<T, APPROX_TYPE>>
Polytope<T, APPROX_TYPE>::split() const
{
  return this->split(this->_A.size());
}

template<typename T, typename APPROX_TYPE>
std::list<Polytope<T, APPROX_TYPE>> Polytope<T, APPROX_TYPE>::split(
    const std::vector<bool> &bvect_base, const unsigned int cidx,
    std::list<Polytope<T, APPROX_TYPE>> &tmp_covering,
    std::vector<std::vector<T>> &A, std::vector<T> &b,
    const unsigned int num_of_splits) const
{
  using namespace LinearAlgebra;

  if (this->_A.size() == cidx) {
    Polytope<T, APPROX_TYPE> ls(A, b);
    ls.simplify();

    tmp_covering.push_back(ls);
    return tmp_covering;
  }

  if (bvect_base[cidx]) {
    A.push_back(this->_A[cidx]);
    b.push_back(this->_b[cidx]);

    try {
      const T min_value = this->minimize(this->_A[cidx]).objective_value();

      A.push_back(-this->_A[cidx]);
      if (num_of_splits > 0) {
        const T avg_value = (this->_b[cidx] + min_value) / 2;
        b.push_back(-avg_value);

        split(bvect_base, cidx + 1, tmp_covering, A, b, num_of_splits - 1);

        b[b.size() - 1] = -min_value;
        b[b.size() - 2] = avg_value;

        split(bvect_base, cidx + 1, tmp_covering, A, b, num_of_splits - 1);
      } else {
        b.push_back(-min_value);
        split(bvect_base, cidx + 1, tmp_covering, A, b, num_of_splits);
      }

      A.pop_back();
      b.pop_back();
    } catch (std::logic_error &e) {
      std::cerr << "This polytope is open." << std::endl;

      split(bvect_base, cidx + 1, tmp_covering, A, b, num_of_splits);
    }

    A.pop_back();
    b.pop_back();
  } else {
    split(bvect_base, cidx + 1, tmp_covering, A, b, num_of_splits);
  }

  return tmp_covering;
}

/**
 * @brief Get a set of linear independent rows in a matrix
 *
 * @tparam T is s numeric type
 * @param A is the matrix whose linear independent rows must be found
 * @return a vector containing the indices of linear independent
 *         rows in the matrix `A`
 */
template<typename T>
std::vector<size_t>
get_a_polytope_base(const LinearAlgebra::Dense::Matrix<T> &A)
{
  using namespace LinearAlgebra;

  if (A.size() == 0) {
    return std::vector<size_t>();
  }

  std::vector<size_t> base{0};

  size_t row_idx = 1;
  while (row_idx < A.size()) {
    bool indep_from_base = true;
    auto b_it = std::begin(base);
    while (indep_from_base && b_it != std::end(base)) {
      indep_from_base = !are_linearly_dependent(A[row_idx], A[*b_it]);
      ++b_it;
    }

    if (indep_from_base) {
      base.push_back(row_idx);
    }
    ++row_idx;
  }

  return base;
}

/**
 * @brief Get the bit-vector of a set of linear independent rows in a matrix
 *
 * @tparam T is s numeric type
 * @param A is the matrix whose linear independent rows must be found
 * @return a Boolean vector such that the indices of all the `true`
 *         values correspond to a linearly independent set in `A`
 */
template<typename T>
std::vector<bool>
get_a_polytope_base_bit_vector(const LinearAlgebra::Dense::Matrix<T> &A)
{
  std::vector<size_t> base = get_a_polytope_base(A);

  std::vector<bool> bvect_base(A.size(), false);

  for (auto it = std::begin(base); it != std::end(base); ++it) {
    bvect_base[*it] = true;
  }

  return bvect_base;
}

template<typename T, typename APPROX_TYPE>
std::list<Polytope<T, APPROX_TYPE>>
Polytope<T, APPROX_TYPE>::split(const unsigned int num_of_splits) const
{
  std::list<Polytope<T, APPROX_TYPE>> result;

  std::vector<bool> bvect_base = get_a_polytope_base_bit_vector(this->_A);

  std::vector<std::vector<T>> A;
  std::vector<T> b;

  split(bvect_base, 0, result, A, b, num_of_splits);

  return result;
}

template<typename T, typename APPROX_TYPE>
Polytope<T, APPROX_TYPE> &
Polytope<T, APPROX_TYPE>::intersect_with(const Polytope<T, APPROX_TYPE> &P)
{
  if (this->dim() != P.dim()) {
    SAPO_ERROR("the intersecting parallelotope differ "
               "in space dimensions",
               std::domain_error);
  }

  for (unsigned int i = 0; i < P.size(); i++) {
    if (!is_true(this->satisfies(P._A[i], P._b[i]))) {
      (this->_A).push_back(P._A[i]);
      (this->_b).push_back(P._b[i]);
    }
  }

  return *this;
}

template<typename T, typename APPROX_TYPE>
T Polytope<T, APPROX_TYPE>::bounding_box_volume() const
{

  std::vector<T> base(this->dim(), 0);
  T vol{1};

  for (unsigned int i = 0; i < this->dim(); i++) {
    base[i] = 1;
    const T b_plus = maximize(base).objective_value();
    const T b_minus = minimize(base).objective_value();
    vol = vol * (b_plus + b_minus);
    base[i] = 0;
  }

  return vol;
}

template<typename T, typename APPROX_TYPE>
inline void swap(Polytope<T, APPROX_TYPE> &P1, Polytope<T, APPROX_TYPE> &P2)
{
  swap(*(static_cast<LinearSystem<T, APPROX_TYPE> *>(&P1)),
       *(static_cast<LinearSystem<T, APPROX_TYPE> *>(&P2)));
}

template<typename T, typename APPROX_TYPE>
inline TriBool are_disjoint(const Polytope<T, APPROX_TYPE> &P1,
                            const Polytope<T, APPROX_TYPE> &P2)
{
  return have_disjoint_solutions(P1, P2);
}

template<class BASIC_SET_TYPE, typename T>
BASIC_SET_TYPE expand(const BASIC_SET_TYPE &S, const T epsilon)
{
  BASIC_SET_TYPE expS(S);

  expS.expand_by(epsilon);

  return expS;
}

template<typename T, typename APPROX_TYPE>
Polytope<T, APPROX_TYPE> intersect(const Polytope<T, APPROX_TYPE> &P1,
                                   const Polytope<T, APPROX_TYPE> &P2)
{
  Polytope<T, APPROX_TYPE> result(P1._A, P1._b);

  result.intersect_with(P2);

  return result;
}

template<typename T, typename APPROX_TYPE>
Polytope<T, APPROX_TYPE>
over_approximate_union(const Polytope<T, APPROX_TYPE> &P1,
                       const Polytope<T, APPROX_TYPE> &P2)
{
  using namespace LinearAlgebra;
  using namespace LinearAlgebra::Dense;

  if (P1.dim() != P2.dim()) {
    SAPO_ERROR("the two polytopes differ in space dimensions",
               std::domain_error);
  }

  if (is_true(P1.is_empty())) {
    return P2;
  }

  if (is_true(P2.is_empty())) {
    return P1;
  }

  Matrix<T> A;
  Vector<T> b;

  A.reserve(P1.size() + P2.size());
  b.reserve(P1.size() + P2.size());

  for (size_t i = 0; i < P1.size(); ++i) {
    A.push_back(P1._A[i]);
    b.push_back(std::max(P2.maximize(A.back()).objective_value(), P1._b[i]));
  }

  for (size_t i = 0; i < P2.size(); ++i) {
    A.push_back(P2._A[i]);
    b.push_back(std::max(P1.maximize(A.back()).objective_value(), P2._b[i]));
  }

  Polytope<T, APPROX_TYPE> res(std::move(A), std::move(b));

  res.simplify();

  return res;
}

#endif /* POLYTOPE_H_ */
