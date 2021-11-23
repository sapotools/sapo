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

std::list<std::vector<double>> get_parallelotope_vertices(
    const std::vector<std::vector<double>> &template_matrix,
    const std::vector<double> &lower_bound,
    const std::vector<double> &upper_bound)
{
  const unsigned int dim = upper_bound.size();

  std::list<std::vector<double>> vertices;

  std::vector<double> offset = upper_bound;

  SparseLinearAlgebra::Matrix<double> tmatrix(template_matrix);
  SparseLinearAlgebra::PLU_Factorization<double> factorization(tmatrix);

  try {
    // store the base vertex
    vertices.push_back(factorization.solve(offset));

  } catch (std::domain_error &) { // if a domain_error is raised, then the
                                  // template is singular
    std::cerr << "singular parallelotope" << std::endl << std::endl;
    exit(EXIT_FAILURE);
  }

  // Compute the vertices v
  for (unsigned int k = 0; k < dim; k++) {
    double tmp = offset[k];
    offset[k] = -lower_bound[k];
    vertices.push_back(factorization.solve(offset));
    offset[k] = tmp;
  }

  return vertices;
}

/**
 * Constructor
 *
 * @param[in] template_matrix is the template matrix of the parallelotope
 * @param[in] lower_bound is the lower bound offsets of the parallelotope
 * @param[in] upper_bound is the upper bound offsets of the parallelotope
 */
Parallelotope::Parallelotope(const Matrix &template_matrix,
                             const Vector &lower_bound,
                             const Vector &upper_bound)
{
  using namespace std;

  // convert the linear system to vectors
  std::list<Vector> vertices
      = get_parallelotope_vertices(template_matrix, lower_bound, upper_bound);

  this->_base_vertex = vertices.front();

  // Compute the generators, their lengths, and the versors
  for (auto v_it = ++std::begin(vertices); v_it != std::end(vertices);
       ++v_it) {

    // add a new generator
    const Vector gen = *v_it - _base_vertex;

    // compute its length and store it
    const double length = norm_2(gen);
    _lengths.push_back(length);

    // compute and store the corresponding versor

    // The approximation below improves performances
    // TODO: set the approximation dynamically, possibly
    //       by using SIL input.
    if (length == 0) {
      // TODO: check when a versor can be null
      _versors.push_back(approx(gen, 11));
    } else {
      _versors.push_back(approx(gen / length, 11));
    }
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

  if (pts.size() == 0 || pts.size() != pts.front().size()) {
    std::cerr << "hyperplane_through_points: pts must contain "
              << "n non-linearly dependent n-dimensional points with n!=0"
              << std::endl;
    exit(EXIT_FAILURE);
  }

  // build the matrix A=[pts[1]-pts[0],...,pts[n-2]-pts[0]]^T
  SparseLinearAlgebra::Matrix<double> A;
  for (auto pts_it = ++std::begin(pts); pts_it != std::end(pts); ++pts_it) {
    A.add_row(pts.front() - *pts_it);
  }

  // factorize the A
  SparseLinearAlgebra::PLU_Factorization<double> fact(A);

  // A is underdetermined because it is a (n-1)x n matrix
  // find the underdetermined dimension udim
  unsigned int udim = 0;
  const SparseLinearAlgebra::Matrix<double> &U = fact.U();
  while (udim < U.num_of_rows() && U[udim][udim] != 0) {
    ++udim;
  }

  // add to A the versor on udim
  SparseLinearAlgebra::Matrix<double>::RowType row;
  row[udim] = 1;
  A.add_row(row);

  // factorize the new matrix
  fact = SparseLinearAlgebra::PLU_Factorization<double>(A);

  std::vector<double> lambda;
  try {
    // try to solve the new factorization: if it is not possible
    // the points did not define an n-dimensional hyperplane.
    lambda = fact.solve(std::vector<double>(A.num_of_rows(), 0));

  } catch (std::domain_error &) {
    std::cerr << "hyperplane_through_points: the points "
              << "were not linearly independent." << std::endl;
    exit(EXIT_FAILURE);
  }

  // add the offset to the result
  const double offset = -(lambda * pts.front());
  lambda.push_back(offset);

  return lambda;
}

/**
 * Constructor generator to constraints representation
 *
 * @param[in] q numerical base vertex
 * @param[in] beta numerical generator lengths
 * @returns polytope representing the parallelotope
 */
Polytope Parallelotope::gen2const(const Vector &q, const Vector &beta) const
{
  using namespace std;

  const unsigned int dim = this->dim();

  if (q.size() != dim) {
    std::cerr << "Parallelotope::gen2const : q must have dimension " << dim
              << std::endl;
    exit(EXIT_FAILURE);
  }

  std::vector<std::vector<double>> hps; // hyperplane equations

  // first set of hyperplanes
  for (unsigned int i = 0; i < dim; i++) {

    std::list<std::vector<double>> pts;
    pts.push_back(q); // add base vertex

    for (unsigned int j = 0; j < dim; j++) { // for all the generators u
      if (i != j) {
        pts.push_back(q + (beta[j] * this->_versors[j]));
      }
    }
    hps.push_back(hyperplane_through_points(pts));
  }

  // second set of hyperplanes
  for (unsigned int i = 0; i < dim; i++) {

    std::list<std::vector<double>> pts;

    const std::vector<double> &versor_i = _versors[i];

    std::vector<double> qt; // traslate q
    for (unsigned int j = 0; j < dim; j++) {
      qt.push_back(q[j] + beta[i] * versor_i[j]);
    }
    pts.push_back(qt); // add base vertex

    for (unsigned int j = 0; j < dim; j++) {
      // for all the generators u
      if (i != j) {
        std::vector<double> p;
        for (unsigned int k = 0; k < dim; k++) {
          // coordinate of the point
          p.push_back((q[k] + this->_versors[j][k] * beta[j])
                      + beta[i] * versor_i[k]);
        }
        pts.push_back(p);
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
  std::swap(A._versors, B._versors);
  std::swap(A._base_vertex, B._base_vertex);
  std::swap(A._lengths, B._lengths);
}
