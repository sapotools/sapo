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

    // The approximation below improves performances
    // TODO: set the approximation dynamically, possibly
    //       by using SIL input.
    _generators.push_back(approx(tensor / length, 11));
  }
}

/**
 * Determine the equation of the hyperplane passing through some linearly
 * independent points, i.e., res[0]*x_0 + .. + res[n-1]*x_{n-1} + res[n] = 0
 *
 * @param[in] pts interpolation points
 * @returns interpolating function coefficients
 */
std::vector<double>
hyperplane_through_points(const std::list<std::vector<double>> &pts)
{
  using namespace std;
  using namespace LinearAlgebra;

  if (pts.size() == 0 || pts.size() != pts.front().size()) {
    std::domain_error("hyperplane_through_points: pts must contain "
                      "n non-linearly dependent n-dimensional points");
  }

  // build the matrix A=[pts[1]-pts[0],...,pts[n-2]-pts[0]]^T
  Sparse::Matrix<double> A;
  for (auto pts_it = ++std::begin(pts); pts_it != std::end(pts); ++pts_it) {
    A.add_row(pts.front() - *pts_it);
  }

  // factorize the A
  Sparse::LUP_Factorization<double> fact(A);

  // A is underdetermined because it is a (n-1)x n matrix
  // find the underdetermined dimension udim
  unsigned int udim = 0;
  const Sparse::Matrix<double> &U = fact.U();
  while (udim < U.num_of_rows() && U[udim][udim] != 0) {
    ++udim;
  }

  // add to A the versor on udim
  Sparse::Matrix<double>::RowType row;
  row[udim] = 1;
  A.add_row(row);

  // factorize the new matrix
  fact = LinearAlgebra::Sparse::LUP_Factorization<double>(A);

  std::vector<double> lambda;
  try {
    // try to solve the new factorization: if it is not possible
    // the points did not define an n-dimensional hyperplane.
    lambda = fact.solve(std::vector<double>(A.num_of_rows(), 0));

  } catch (std::domain_error &) {
    std::domain_error("hyperplane_through_points: the points "
                      "are linearly dependent.");
  }

  // add the offset to the result
  const double offset = -(lambda * pts.front());
  lambda.push_back(offset);

  return lambda;
}

// TODO: Fix the method below
/**
 * Build a polytope representing the parallelotope.
 *
 * @returns a polytope representing the parallelotope
 */
Parallelotope::operator Polytope() const
{
  using namespace std;
  using namespace LinearAlgebra;

  const unsigned int dim = this->dim();

  std::vector<std::vector<double>> hps; // hyperplane equations

  // first set of hyperplanes
  for (unsigned int i = 0; i < dim; i++) {

    std::list<std::vector<double>> pts;
    pts.push_back(_base_vertex); // add base vertex

    for (unsigned int j = 0; j < dim; j++) { // for all the generators u
      if (i != j) {
        if (_lengths[j] != 0) {
          pts.push_back(pts.front() + _lengths[j] * this->_generators[j]);
        } else {
          // these are fake points to define one hypeplane in
          // degenerous parallelotypes, for instance, when one
          // of the variables admits one single value.
          pts.push_back(pts.front() + this->_generators[j]);
        }
      }
    }
    hps.push_back(hyperplane_through_points(pts));
  }

  // second set of hyperplanes
  for (unsigned int i = 0; i < dim; i++) {

    std::list<std::vector<double>> pts;

    // add base vertex
    pts.push_back(_base_vertex + _lengths[i] * _generators[i]);

    for (unsigned int j = 0; j < dim; j++) {
      // for all the generators u
      if (i != j) {
        if (_lengths[j] != 0) {
          pts.push_back(pts.front() + _lengths[j] * this->_generators[j]);
        } else {
          // these are fake points to define one hypeplane in
          // degenerous parallelotypes, for instance, when one
          // of the variables admits one single value.
          pts.push_back(pts.front() + this->_generators[j]);
        }
      }
    }

    hps.push_back(hyperplane_through_points(pts));
  }

  std::vector<double> d(dim * 2, 0);
  std::vector<std::vector<double>> Lambda(d.size(),
                                          std::vector<double>(dim, 0));

  for (unsigned int i = 0; i < dim; i++) {

    const std::vector<double> &hps_i = hps[i];
    const std::vector<double> &hps_i_plus = hps[i + dim];

    std::vector<double> &Lambda_i = Lambda[i];
    std::vector<double> &Lambda_i_plus = Lambda[i + dim];

    // find hyperplane with smaller direction
    double b1 = -hps_i[hps[i].size() - 1];
    double b2 = -hps_i_plus[hps[i].size() - 1];

    double sign = (b1 > b2 ? 1 : -1);

    d[i] = sign * b1;
    d[i + dim] = -sign * b2;
    for (unsigned int j = 0; j < dim; j++) {
      Lambda_i[j] = sign * hps_i[j];
      Lambda_i_plus[j] = -sign * hps_i_plus[j];
    }
  }

  return Polytope(Lambda, d);
}

void swap(Parallelotope &A, Parallelotope &B)
{
  std::swap(A._generators, B._generators);
  std::swap(A._base_vertex, B._base_vertex);
  std::swap(A._lengths, B._lengths);
}
