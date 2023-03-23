/**
 * @file Parallelotope.cpp
 * @author Tommaso Dreossi (tommasodreossi@berkeley.edu)
 * @author Alberto Casagrande <acasagrande@units.it>
 * @brief Represent and manipulate a parallelotope
 * @version 0.1
 * @date 2015-10-14
 *
 * @copyright Copyright (c) 2015-2022
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
  std::vector<double> offset(upper_bound.size(), 0);
  for (unsigned int k = 0; k < upper_bound.size(); k++) {
    const double delta_k = upper_bound[k] - lower_bound[k];
    offset[k] = 1;

    const std::vector<double> tensor = factorization.solve(offset);
    offset[k] = 0;

    // compute its length and store it
    const double length = norm_2(tensor);
    if (delta_k != 0) {
      _lengths.push_back(delta_k * length);
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
Parallelotope::operator Polytope<double>() const
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

/**
 * @brief Get the parallelotope's center
 *
 * @return the parallelotope's center
 */
LinearAlgebra::Vector<double> Parallelotope::center() const
{
  using namespace LinearAlgebra;

  Vector<double> c{base_vertex()};

  for (size_t i = 0; i < dim(); ++i) {
    c = c + _lengths[i] / 2 * _generators[i];
  }

  return c;
}

void swap(Parallelotope &A, Parallelotope &B)
{
  std::swap(A._generators, B._generators);
  std::swap(A._base_vertex, B._base_vertex);
  std::swap(A._lengths, B._lengths);
}
