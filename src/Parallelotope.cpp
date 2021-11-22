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
 * Constructor that instantiates a parallelotope
 *
 * @param[in] vars vector with the list of variables used for the generator
 * functions
 * @param[in] u collection of generator versors
 */
Parallelotope::Parallelotope(const std::vector<GiNaC::lst> &vars,
                             const Matrix &u)
{
  using namespace std;
  using namespace GiNaC;

  if (vars.size() != 3) {
    cerr << "Parallelotope::Parallelotope : vars must contain 3 "
         << "collections of variable names (q,alpha,beta)" << endl;
    exit(EXIT_FAILURE);
  }

  this->vars.push_back(vars[0]);
  this->vars.push_back(vars[1]);
  this->vars.push_back(vars[2]);

  // get the dimension of the parallelotope
  this->dim = vars[0].nops();
  // and store its variable names
  for (int i = 0; i < 3; i++) {
    if (vars[i].nops() != this->dim) {
      std::cerr << "Parallelotope::Parallelotope : vars[" << i
                << "] must have " << this->dim << " variables" << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  // extract variable names
  lst q = this->vars[0];
  lst alpha = this->vars[1];
  lst beta = this->vars[2];

  // initialize generator function
  for (auto it = std::begin(q); it != std::end(q); ++it) {
    this->generator_function.append(*it);
  }

  // create the generation function accumulating the versor values
  for (unsigned int i = 0; i < this->dim; i++) {

    // check dimension of versors
    if (u[i].size() != this->dim) {
      std::cerr << "Parallelotope::Parallelotope : dim and ui "
                << "dimensions must agree" << std::endl;
      exit(EXIT_FAILURE);
    }

    // Normalize the versors and generatre the generator function
    // double norm = euclidNorm(u[i]);
    // std::vector< double > norm_versor;

    for (unsigned int j = 0; j < this->dim; j++) {
      // norm_versor.push_back(u[i][j]/norm);
      this->generator_function[j]
          = this->generator_function[j] + alpha[i] * beta[i] * u[i][j];
    }
    // store the i-th versor
    this->u.push_back(u[i]);
  }

  // Initialize the template matrix
  std::vector<double> base_vertex(this->dim, 0);
  Vector lengths(this->dim, 1);
}

/**
 * Constructor that instantiates a parallelotope from a polytope
 *
 * @param[in] vars vector with the list of variables used for the generator
 * functions
 * @param[in] constr polytope constituting the parallelotope
 */
Parallelotope::Parallelotope(const std::vector<GiNaC::lst> &vars,
                             const Polytope &constr):
    Parallelotope(vars, constr.getA(), constr.getb())
{
}

std::list<std::vector<double>> get_parallelotope_vertices(
    const std::vector<std::vector<double>> &template_matrix,
    const std::vector<double> &offset)
{
  const unsigned int dim = offset.size() / 2;

  std::list<std::vector<double>> vertices;

  std::vector<double> offset_cut;
  std::copy(offset.begin(), offset.begin() + dim,
            std::back_inserter(offset_cut));

  SparseLinearAlgebra::Matrix<double> tmatrix(template_matrix, dim);
  SparseLinearAlgebra::PLU_Factorization<double> factorization(tmatrix);

  try {
    // store the base vertex
    vertices.push_back(factorization.solve(offset_cut));

  } catch (std::domain_error &) { // if a domain_error is raised, then the
                                  // template is singular
    std::cerr << "singular parallelotope" << std::endl
              << Polytope(template_matrix, offset) << std::endl
              << std::endl;
    exit(EXIT_FAILURE);
  }

  // Compute the vertices v
  for (unsigned int k = 0; k < dim; k++) {
    double tmp = offset_cut[k];
    offset_cut[k] = -offset[k + dim];
    vertices.push_back(factorization.solve(offset_cut));
    offset_cut[k] = tmp;
  }

  return vertices;
}

// TODO: this method assumes that the template matrix has the form
//       [A, -A]^T where A is square. It would be safier to have as
//       parameters a square template matrix A, the lower offset
//       boundaries, and the upper offset boundaries.
/**
 * Constructor that instantiates a parallelotope from a linear system
 *
 * @param[in] vars vector with the list of variables used for the generator
 * functions
 * @param[in] template_matrix is the template matrix of the parallelotope
 * @param[in] offset is the offset of the parallelotope
 */
Parallelotope::Parallelotope(const std::vector<GiNaC::lst> &vars,
                             const Matrix &template_matrix,
                             const Vector &offset)
{
  using namespace std;
  using namespace GiNaC;

  if (vars.size() != 3) {
    cerr << "Parallelotope::Parallelotope : vars must contain "
         << "3 collections of variable names (q,alpha,beta)" << endl;
    exit(EXIT_FAILURE);
  }

  this->vars.push_back(vars[0]);
  this->vars.push_back(vars[1]);
  this->vars.push_back(vars[2]);

  // get the dimension of the parallelotope
  this->dim = vars[0].nops();
  // and store its variable names
  for (int i = 0; i < 3; i++) {
    if (vars[i].nops() != this->dim) {
      std::cerr << "Parallelotope::Parallelotope : vars[" << i
                << "] must have " << this->dim << " variables" << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  // initialize generator function
  for (auto it = std::begin(vars[0]); it != std::end(vars[0]); ++it) {
    this->generator_function.append(*it);
  }

  // convert the linear system to vectors
  std::list<Vector> vertices
      = get_parallelotope_vertices(template_matrix, offset);

  this->base_vertex = vertices.front();

  // Compute the generators, their lengths, and the versors
  for (auto v_it = ++std::begin(vertices); v_it != std::end(vertices);
       ++v_it) {

    // add a new generator
    const Vector gen = *v_it - base_vertex;

    // compute its length and store it
    const double length = norm_2(gen);
    lengths.push_back(length);

    // compute and store the corresponding versor

    // The approximation below improves performances
    // TODO: check whether the performance-approximation
    //       corelation is due to GiNaC.
    // TODO: set the approximation dynamically, possibly
    //       by using SIL input.
    if (length == 0) {
      // TODO: check when a versor can be null
      u.push_back(approx(gen, 11));
    } else {
      u.push_back(approx(gen / length, 11));
    }
  }

  // create the generation function accumulating the versor values
  const lst &alpha = this->vars[1];
  const lst &beta = this->vars[2];
  for (unsigned int i = 0; i < this->dim; i++) {

    // check dimension of versors
    if (u[i].size() != this->dim) {
      std::cerr << "Parallelotope::Parallelotope : dim and ui "
                << "dimensions must agree" << std::endl;
      exit(EXIT_FAILURE);
    }

    // Generatre the generator function
    // double norm = euclidNorm(u[i]);
    // std::vector< double > norm_versor;

    for (unsigned int j = 0; j < this->dim; j++) {
      // norm_versor.push_back(u[i][j]/norm);
      this->generator_function[j]
          = this->generator_function[j] + alpha[i] * beta[i] * this->u[i][j];
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

  if (q.size() != this->dim) {
    std::cerr << "Parallelotope::gen2const : q must have dimension "
              << this->dim << std::endl;
    exit(EXIT_FAILURE);
  }

  std::vector<std::vector<double>> hps; // hyperplane equations

  // first set of hyperplanes
  for (unsigned int i = 0; i < this->dim; i++) {

    std::list<std::vector<double>> pts;
    pts.push_back(q); // add base vertex

    for (unsigned int j = 0; j < this->dim; j++) { // for all the generators u
      if (i != j) {
        pts.push_back(q + (beta[j] * this->u[j]));
      }
    }
    hps.push_back(hyperplane_through_points(pts));
  }

  // second set of hyperplanes
  for (unsigned int i = 0; i < this->dim; i++) {

    std::list<std::vector<double>> pts;

    std::vector<double> qt; // traslate q
    for (unsigned int j = 0; j < this->dim; j++) {
      qt.push_back(q[j] + beta[i] * this->u[i][j]);
    }
    pts.push_back(qt); // add base vertex

    for (unsigned int j = 0; j < this->dim; j++) {
      // for all the generators u
      if (i != j) {
        std::vector<double> p;
        for (unsigned int k = 0; k < this->dim; k++) {
          // coordinate of the point
          p.push_back((q[k] + this->u[j][k] * beta[j])
                      + beta[i] * this->u[i][k]);
        }
        pts.push_back(p);
      }
    }

    hps.push_back(hyperplane_through_points(pts));
  }

  std::vector<std::vector<double>> Lambda;
  std::vector<double> d(this->dim * 2, 0);

  // initialize template Lambda
  std::vector<double> Lambda_i(this->dim, 0);
  for (unsigned int i = 0; i < (this->dim) * 2; i++) {
    Lambda.push_back(Lambda_i);
  }

  for (unsigned int i = 0; i < this->dim; i++) {

    // find hyperplane with smaller direction
    double b1 = -hps[i][hps[i].size() - 1];
    double b2 = -hps[i + this->dim][hps[i].size() - 1];

    if (b1 > b2) {
      d[i] = b1;
      d[i + this->dim] = -b2;
      for (unsigned int j = 0; j < this->dim; j++) {
        Lambda[i][j] = hps[i][j];
        Lambda[i + this->dim][j] = -hps[i + this->dim][j];
      }
    } else {
      d[i] = -b1;
      d[i + this->dim] = b2;
      for (unsigned int j = 0; j < this->dim; j++) {
        Lambda[i][j] = -hps[i][j];
        Lambda[i + this->dim][j] = hps[i + this->dim][j];
      }
    }
  }

  return Polytope(Lambda, d);
}

/**
 * Convert the constraint representation to the generator one
 *
 * @param[in] constr constraint representation of the parallelotope
 * @returns numerical values to plug in the generator function for the
 * generator representation
 */
poly_values const2gen(const Polytope &P)
{
  using namespace std;

  // convert the linear system to vectors
  auto vertices = get_parallelotope_vertices(P.getA(), P.getb());

  // Compute the generators and their lengths
  std::vector<double> lengths;
  for (auto v_it = ++std::begin(vertices); v_it != std::end(vertices);
       ++v_it) {

    // compute its length and store it
    lengths.push_back(norm_2(*v_it - vertices.front()));
  }

  // Return the conversion
  return poly_values{vertices.front(), lengths};
}

/**
 * Convert a symbolic list of numbers to a vector
 *
 * @param[in] symbolic list
 * @returns numeric list of numbers
 */
std::vector<double> lst2vec(const GiNaC::ex &list)
{
  using namespace GiNaC;

  std::vector<double> res;

  for (unsigned int i = 0; i < list.nops(); i++)
    res.push_back(ex_to<numeric>(evalf(list[i])).to_double());

  return res;
}

Parallelotope::~Parallelotope()
{
  // TODO Auto-generated destructor stub
}

void swap(Parallelotope &A, Parallelotope &B)
{
  std::swap(A.dim, B.dim);
  std::swap(A.vars, B.vars);
  std::swap(A.generator_function, B.generator_function);
  std::swap(A.u, B.u);
}
