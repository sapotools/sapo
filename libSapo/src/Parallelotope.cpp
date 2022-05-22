/**
 * @file Parallelotope.cpp
 * Represent and manipulate a parallelotope
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#include "Parallelotope.h"

#include "LinearAlgebra.h"

#include <cmath>
#include <sstream>

/**
 * Constructor
 *
 * @param[in] directions is the vector of the parallelotope directions
 * @param[in] lower_bound is the lower bound offsets of the parallelotope
 * @param[in] upper_bound is the upper bound offsets of the parallelotope
 */
Parallelotope::Parallelotope(
    const std::vector<LinearAlgebra::Vector<double>> &directions,
    const LinearAlgebra::Vector<double> &lower_bound,
    const LinearAlgebra::Vector<double> &upper_bound)
{
  using namespace LinearAlgebra;

  if (lower_bound.size() != upper_bound.size()) {
    std::ostringstream oss;

    oss << "The lower and upper bounds vectors must "
        << "have the same dimension: they are " << lower_bound.size()
        << " and " << upper_bound.size() << ", respectively.";
    throw std::domain_error(oss.str());
  }

  if (lower_bound.size() != directions.size()) {
    std::ostringstream oss;

    oss << "The lower bounds vector and the direction vector must "
        << "have the same dimension: they are " << lower_bound.size()
        << " and " << directions.size() << ", respectively.";
    throw std::domain_error(oss.str());
  }

  const Sparse::Matrix<double> dmatrix(directions);
  Sparse::LUP_Factorization<double> factorization(dmatrix);

  try {
    // store the base vertex
    _base_vertex = factorization.solve(lower_bound);

  } catch (std::domain_error &) { // if a domain_error is raised, then the
                                  // template is singular
    std::domain_error("Parallelotope::Parallelotope: singular parallelotope");
  }

  // Compute the versors
  std::vector<double> offset(upper_bound.size(), 0);
  for (unsigned int k = 0; k < upper_bound.size(); k++) {
    const double delta_k = upper_bound[k] - lower_bound[k];
    if (delta_k != 0) {
      offset[k] = delta_k;
    } else {
      offset[k] = -1;
    }

    const std::vector<double> tensor = factorization.solve(offset);
    offset[k] = 0;

    // compute its length and store it
    const double length = norm_2(tensor);
    if (delta_k != 0) {
      _lengths.push_back(length);
    } else {
      _lengths.push_back(0);
    }

    // compute and store the corresponding versor
    _generators.push_back(tensor / length);
  }
}

/**
 * Build a polytope representing the parallelotope.
 *
 * @returns a polytope representing the parallelotope
 */
Parallelotope::operator Polytope() const
{
  using namespace LinearAlgebra;

  std::vector<Vector<double>> A(2 * dim());
  Vector<double> b(2 * dim());
  for (unsigned int i = 0; i < dim(); ++i) {
    const double base_bound = _generators[i] * _base_vertex;
    const double vertex_bound = base_bound + _lengths[i];

    A[i] = Vector<double>(_generators[i]);
    A[i + dim()] = -_generators[i];
    b[i] = std::max(base_bound, vertex_bound);
    b[i + dim()] = -std::min(base_bound, vertex_bound);
  }

  return Polytope(A, b);
}

void swap(Parallelotope &A, Parallelotope &B)
{
  std::swap(A._generators, B._generators);
  std::swap(A._base_vertex, B._base_vertex);
  std::swap(A._lengths, B._lengths);
}
