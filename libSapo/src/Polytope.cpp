/**
 * @file Polytope.cpp
 * @author <acasagrande@units.it>
 * @brief  Represent and manipulate polytopes (reached states, parameters,
 * etc.)
 * @version 0.1
 * @date 2021-11-16
 *
 * @copyright Copyright (c) 2021-2022
 */

#include "Polytope.h"
#include "PolytopesUnion.h"

#include "LinearAlgebra.h"
#include "SetsUnion.h"

#include "ErrorHandling.h"

/**
 * @brief Get a set of linear independent rows in a matrix
 *
 * @param A is the matrix whose linear independent rows must be found
 * @return a vector containing the indices of linear independent
 *         rows in the matrix `A`
 */
std::vector<unsigned int>
get_a_polytope_base(const LinearAlgebra::Dense::Matrix<double> &A)
{
  using namespace LinearAlgebra;

  if (A.size() == 0) {
    return std::vector<unsigned int>();
  }

  std::vector<unsigned int> base{0};

  unsigned int row_idx = 1;
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
 * @param A is the matrix whose linear independent rows must be found
 * @return a Boolean vector such that the indices of all the `true`
 *         values correspond to a linearly independent set in `A`
 */
std::vector<bool>
get_a_polytope_base_bit_vector(const LinearAlgebra::Dense::Matrix<double> &A)
{
  std::vector<unsigned int> base = get_a_polytope_base(A);

  std::vector<bool> bvect_base(A.size(), false);

  for (auto it = std::begin(base); it != std::end(base); ++it) {
    bvect_base[*it] = true;
  }

  return bvect_base;
}

std::list<Polytope> Polytope::split(const std::vector<bool> &bvect_base,
                                    const unsigned int cidx,
                                    std::list<Polytope> &tmp_covering,
                                    std::vector<std::vector<double>> &A,
                                    std::vector<double> &b,
                                    const unsigned int num_of_splits) const
{
  using namespace LinearAlgebra;

  if (this->_A.size() == cidx) {
    Polytope ls(A, b);
    ls.simplify();

    tmp_covering.push_back(ls);
    return tmp_covering;
  }

  if (bvect_base[cidx]) {
    A.push_back(this->_A[cidx]);
    b.push_back(this->_b[cidx]);

    try {
      const double min_value = minimize(this->_A[cidx]).optimum();

      A.push_back(-this->_A[cidx]);
      if (num_of_splits > 0) {
        const double avg_value = (this->_b[cidx] + min_value) / 2;
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

std::list<Polytope> Polytope::split(const unsigned int num_of_splits) const
{
  std::list<Polytope> result;

  std::vector<bool> bvect_base = get_a_polytope_base_bit_vector(this->_A);

  std::vector<std::vector<double>> A;
  std::vector<double> b;

  split(bvect_base, 0, result, A, b, num_of_splits);

  return result;
}

/**
 * Update a polytope by intersecting it with another one.
 *
 * This method works in-place and changes the calling object.
 *
 * @param[in] P is a polytope.
 * @return a reference to the updated object.
 */
Polytope &Polytope::intersect_with(const Polytope &P)
{
  if (this->dim() != P.dim()) {
    SAPO_ERROR("the intersecting parallelotope differ "
               "in space dimensions",
               std::domain_error);
  }

  for (unsigned int i = 0; i < P.size(); i++) {
    if (!this->satisfies(P._A[i], P._b[i])) {
      (this->_A).push_back(P._A[i]);
      (this->_b).push_back(P._b[i]);
    }
  }

  return *this;
}

Polytope &Polytope::expand_by(const double delta)
{
  for (auto b_it = std::begin(_b); b_it != std::end(_b); ++b_it) {
    *b_it += delta;
  }

  return *this;
}

Polytope intersect(const Polytope &P1, const Polytope &P2)
{
  Polytope result(P1._A, P1._b);

  result.intersect_with(P2);

  return result;
}

/**
 * Determine the volume of the bounding box of the polytope
 *
 * @return volume of the bounding box
 */
double Polytope::bounding_box_volume() const
{

  std::vector<double> zeros(this->dim(), 0);
  double vol = 1;

  for (unsigned int i = 0; i < this->dim(); i++) {
    std::vector<double> facet = zeros;
    facet[i] = 1;
    const double b_plus = maximize(facet).optimum();
    facet[i] = -1;
    const double b_minus = minimize(facet).optimum();
    vol = vol * (b_plus + b_minus);
  }

  return vol;
}

Polytope over_approximate_union(const Polytope &P1, const Polytope &P2)
{
  using namespace LinearAlgebra;
  using namespace LinearAlgebra::Dense;

  if (P1.dim() != P2.dim()) {
    SAPO_ERROR("the two polytopes differ in space dimensions",
               std::domain_error);
  }

  if (P1.is_empty()) {
    return P2;
  }

  if (P2.is_empty()) {
    return P1;
  }

  Matrix<double> A;
  Vector<double> b;

  A.reserve(P1.size() + P2.size());
  b.reserve(P1.size() + P2.size());

  for (unsigned int i = 0; i < P1.size(); ++i) {
    A.push_back(P1._A[i]);
    b.push_back(std::max(P2.maximize(A.back()).optimum(), P1._b[i]));
  }

  for (unsigned int i = 0; i < P2.size(); ++i) {
    A.push_back(P2._A[i]);
    b.push_back(std::max(P1.maximize(A.back()).optimum(), P2._b[i]));
  }

  Polytope res(std::move(A), std::move(b));

  res.simplify();

  return res;
}

/**
 * @brief Subtract two polytopes and close the result
 *
 * @param[in] P1 is a polytope
 * @param[in] P2 is a polytope
 * @return a union of polytopes obtained by closing the set
 *         \f$P1\setminus P2\f$
 */
SetsUnion<Polytope> subtract_and_close(const Polytope &P1, const Polytope &P2)
{
  using namespace LinearAlgebra;

  SetsUnion<Polytope> su;

  if (P2.includes(P1)) {
    return su;
  }

  if (are_disjoint(P1, P2)) {
    su.add(P1);

    return su;
  }

  for (unsigned int i = 0; i < P2.size(); ++i) {
    auto min = P1.minimize(P2.A(i)).optimum();

    if (min < P2.b(1)) {
      Polytope new_P1 = P1;
      new_P1._A.push_back(-P2.A(i));
      new_P1._b.push_back(-P2.b(1));

      su.add(std::move(new_P1));
    }
  }

  return su;
}