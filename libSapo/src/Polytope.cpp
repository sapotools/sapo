/**
 * @file Polytope.cpp
 * Represent and manipulate polytopes (reached states, parameters, etc.)
 *
 * @author Alberto Casagrande <acasagrande@units.it>
 * @version 0.1
 */

#include "Polytope.h"
#include "PolytopesUnion.h"

#include "LinearAlgebra.h"

using namespace std;

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
  for (unsigned int i = 0; i < P.size(); i++) {
    if (!this->satisfies(P._A[i], P._b[i])) {
      (this->_A).push_back(P._A[i]);
      (this->_b).push_back(P._b[i]);
    }
  }

  return *this;
}

Polytope &Polytope::expand_by(const double epsilon)
{
  for (auto b_it = std::begin(_b); b_it != std::end(_b); ++b_it) {
    *b_it += epsilon;
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

  vector<double> zeros(this->dim(), 0);
  double vol = 1;

  for (unsigned int i = 0; i < this->dim(); i++) {
    vector<double> facet = zeros;
    facet[i] = 1;
    const double b_plus = maximize(facet).optimum();
    facet[i] = -1;
    const double b_minus = minimize(facet).optimum();
    vol = vol * (b_plus + b_minus);
  }

  return vol;
}
