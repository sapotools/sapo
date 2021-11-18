/**
 * @file Parallelotope.cpp
 * Represent and manipulate a parallelotope
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#include "Parallelotope.h"

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
    // vector< double > norm_versor;

    for (unsigned int j = 0; j < this->dim; j++) {
      // norm_versor.push_back(u[i][j]/norm);
      this->generator_function[j]
          = this->generator_function[j] + alpha[i] * beta[i] * u[i][j];
    }
    // store the i-th versor
    this->u.push_back(u[i]);
  }

  // Initialize the template matrix
  vector<double> base_vertex(this->dim, 0);
  Vector lenghts(this->dim, 1);
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
  vector<Vector> vertices;

  // find base vertex
  // build the linear system
  ex q = this->vars[0];
  lst LS;

  //	for (int i=0; i<Lambda.size(); i++) {
  //		for (int j=0; j<Lambda[i].size(); j++) {
  //			cout<<Lambda[i][j]<<" ";
  //		}
  //		cout<<"\n";
  //	}

  for (unsigned int i = 0; i < this->dim; i++) {
    ex eq = 0;
    for (unsigned int j = 0; j < template_matrix[i].size(); j++) {
      eq = eq + template_matrix[i][j] * q[j];
    }
    eq = eq == offset[i];
    LS.append(eq);
  }

  //	cout<<LS;

  ex solLS = lsolve(LS, q);
  if (solLS.nops() == 0) { // the template is singular
    std::cerr << "singular parallelotope" << std::endl
              << Polytope(template_matrix, offset) << std::endl
              << LS << std::endl;
    exit(EXIT_FAILURE);
  }

  vertices.push_back(lst2vec((ex)q.subs(solLS))); // store the base_vertex

  // Compute the vertices v
  for (unsigned int k = 0; k < this->dim; k++) {
    ex a = this->vars[1];
    lst LS;

    for (unsigned int i = 0; i < this->dim; i++) {
      ex eq = 0;
      for (unsigned int j = 0; j < template_matrix[i].size(); j++) {
        eq = eq + template_matrix[i][j] * a[j];
      }
      if (i != k) {
        eq = eq == offset[i];
      } else {
        eq = eq == -offset[i + this->dim];
      }
      LS.append(eq);
    }
    ex solLS = lsolve(LS, a);
    vertices.push_back(lst2vec(a.subs(solLS))); // store the i-th vertex
  }

  this->base_vertex = vertices[0];

  // Compute the generators
  Matrix g(this->dim, Vector(this->dim));
  const Vector &first_vertex = vertices[0];
  for (unsigned int i = 0; i < this->dim; i++) {
    Vector &g_i = g[i];
    const Vector &vertex = vertices[i + 1];
    for (unsigned int j = 0; j < this->dim; j++) {
      g_i[j] = vertex[j] - first_vertex[j];
    }
  }

  lenghts.resize(this->dim);

  // Compute the generators lengths
  for (unsigned int i = 0; i < this->dim; i++) {
    lenghts[i] = euclidNorm(g[i]);
    // cout<<lengths[i]<<" ";
  }
  //	cout<<"\n";

  // cout<<"g\n";
  // Find the versors
  this->u = Matrix(this->dim, Vector(this->dim));
  for (unsigned int i = 0; i < this->dim; i++) {
    Vector &versor_i = this->u[i];
    Vector &g_i = g[i];
    for (unsigned int j = 0; j < this->dim; j++) {
      if (lenghts[i] != 0) {
        versor_i[j] = floor((g_i[j] / lenghts[i]) * 100000000000.0f)
                      / 100000000000.0f;
      } else {
        versor_i[j] = floor(g_i[j] * 100000000000.0f) / 100000000000.0f;
      }
    }
  }

  //	for (int i=0; i<this->u.size(); i++) {
  //		for (int j=0; j<this->u[i].size(); j++) {
  //			cout<<this->u[i][j]<<" ";
  //		}
  //		cout<<"\n";
  //	}

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
    // vector< double > norm_versor;

    for (unsigned int j = 0; j < this->dim; j++) {
      // norm_versor.push_back(u[i][j]/norm);
      this->generator_function[j]
          = this->generator_function[j] + alpha[i] * beta[i] * this->u[i][j];
    }
  }
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
  using namespace GiNaC;

  if (q.size() != this->dim) {
    std::cerr << "Parallelotope::gen2const : q must have dimension "
              << this->dim << std::endl;
    exit(EXIT_FAILURE);
  }

  vector<vector<double>> hps; // hyperplane equations

  // first set of hyperplanes
  for (unsigned int i = 0; i < this->dim; i++) {

    vector<vector<double>> pts;
    pts.push_back(q); // add base vertex

    for (unsigned int j = 0; j < this->dim; j++) { // for all the generators u
      if (i != j) {
        vector<double> p;
        for (unsigned int k = 0; k < this->dim;
             k++) { // coordinate of the point
          p.push_back(q[k] + this->u[j][k] * beta[j]);
        }
        pts.push_back(p);
      }
    }
    hps.push_back(hyperplaneThroughPts(pts));
  }

  // second set of hyperplanes
  for (unsigned int i = 0; i < this->dim; i++) {

    vector<vector<double>> pts;

    vector<double> qt; // traslate q
    for (unsigned int j = 0; j < this->dim; j++) {
      qt.push_back(q[j] + beta[i] * this->u[i][j]);
    }
    pts.push_back(qt); // add base vertex

    for (unsigned int j = 0; j < this->dim; j++) { // for all the generators u
      if (i != j) {
        vector<double> p;
        for (unsigned int k = 0; k < this->dim;
             k++) { // coordinate of the point
          p.push_back((q[k] + this->u[j][k] * beta[j])
                      + beta[i] * this->u[i][k]);
        }
        pts.push_back(p);
      }
    }

    hps.push_back(hyperplaneThroughPts(pts));
  }

  vector<vector<double>> Lambda;
  vector<double> d(this->dim * 2, 0);

  // initialize template Lambda
  vector<double> Lambda_i(this->dim, 0);
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
 * Determine the equation of the hyperplane passing through some points, i.e.,
 * res[0]*x_0 + .. + res[n]*x_n + res[n+1] = 0
 *
 * @param[in] pts interpolation points
 * @returns interpolating function coefficients
 */
std::vector<double> Parallelotope::hyperplaneThroughPts(
    const std::vector<std::vector<double>> &pts) const
{
  using namespace std;
  using namespace GiNaC;

  if (pts.size() != pts[0].size()) {
    std::cerr << "Parallelotope::hyperplaneThroughPts : pts must contain "
              << "n n-dimensional points" << std::endl;
    exit(EXIT_FAILURE);
  }

  Matrix A;
  // build the linear system Ax = 0
  for (unsigned int i = 1; i < pts.size(); i++) {
    A.push_back(Vector(this->dim, 0));
    Vector &last_row = A.back();
    for (unsigned int j = 0; j < this->dim; j++) {
      last_row[j] = pts[0][j] - pts[i][j];
    }
  }

  // Build the linear system to find the normal vector
  lst LS;
  lst a = this->vars[1];

  for (unsigned int i = 0; i < this->dim - 1; i++) {

    ex eq = 0;
    lst sub;

    for (unsigned int j = 0; j < this->dim; j++) {
      eq = eq + A[i][j] * a[j];
    }

    eq = eq == 0;
    LS.append(eq);
  }

  // Solve the linear system
  ex solLS = lsolve(LS, a);

  // Build the linear inequality
  lst x = this->vars[0];
  ex eqr = 0;
  vector<double> lambda(this->dim + 1, 0);

  ex sub;
  unsigned int sub_idx = 0;
  // search for the tautology
  for (unsigned int i = 0; i < solLS.nops(); i++) {
    if (solLS[i].is_equal(a[i] == a[i])) {
      sub = a[i] == 1;
      sub_idx = i;
    }
  }

  lambda[sub_idx] = 1;
  for (unsigned int i = 0; i < this->dim; i++) {
    if (i != sub_idx) {
      lambda[i]
          = ex_to<numeric>(evalf(a[i].subs(solLS[i].subs(sub)))).to_double();
    }
    eqr = eqr + a[i].subs(a[i] == lambda[i]) * pts[0][i];
  }

  lambda[this->dim] = -ex_to<numeric>(evalf(eqr)).to_double();

  return lambda;
}

/**
 * Convert the constraint representation to the generator one
 *
 * @param[in] constr constraint representation of the parallelotope
 * @returns numerical values to plug in the generator function for the
 * generator representation
 */
poly_values Parallelotope::const2gen(Polytope *constr) const
{
  using namespace std;
  using namespace GiNaC;

  const vector<vector<double>> &Lambda = constr->getA();
  const vector<double> &d = constr->getb();
  vector<vector<double>> vertices;

  // find base vertex
  // build the linear system
  ex q = this->vars[0];
  lst LS;

  for (unsigned int i = 0; i < this->dim; i++) {
    ex eq = 0;
    for (unsigned int j = 0; j < Lambda[i].size(); j++) {
      eq = eq + Lambda[i][j] * q[j];
    }
    eq = eq == d[i];
    LS.append(eq);
  }

  ex solLS = lsolve(LS, q);
  vertices.push_back(lst2vec(q.subs(solLS))); // store the base_vertex

  // Compute the vertices v
  for (unsigned int k = 0; k < this->dim; k++) {
    ex a = this->vars[1];
    lst LS;

    for (unsigned int i = 0; i < this->dim; i++) {
      ex eq = 0;
      for (unsigned int j = 0; j < Lambda[i].size(); j++) {
        eq = eq + Lambda[i][j] * a[j];
      }
      if (i != k) {
        eq = eq == d[i];
      } else {
        eq = eq == -d[i + this->dim];
      }
      LS.append(eq);
    }

    ex solLS = lsolve(LS, a);
    vertices.push_back(lst2vec(a.subs(solLS))); // store the i-th vertex
  }

  // Compute the generators
  vector<vector<double>> g;
  for (unsigned int i = 0; i < this->dim; i++) {
    vector<double> gi;
    for (unsigned int j = 0; j < this->dim; j++) {
      gi.push_back(vertices[i + 1][j] - vertices[0][j]);
    }
    g.push_back(gi);
  }

  // Compute the generators lengths
  vector<double> lengths;
  for (unsigned int i = 0; i < this->dim; i++) {
    lengths.push_back(euclidNorm(g[i]));
  }

  // Find the versors (optional since versors are specified by the user)
  vector<vector<double>> versors;
  for (unsigned int i = 0; i < this->dim; i++) {
    vector<double> versori;
    for (unsigned int j = 0; j < this->dim; j++) {
      versori.push_back(g[i][j] / lengths[i]);
    }
    versors.push_back(versori);
  }

  // Return the conversion
  poly_values result;
  result.base_vertex = vertices[0];
  result.lenghts = lengths;
  // result.versors = versors;

  return result;
}

/**
 * Convert a symbolic list of numbers to a vector
 *
 * @param[in] symbolic list
 * @returns numeric list of numbers
 */
std::vector<double> Parallelotope::lst2vec(const GiNaC::ex &list) const
{
  using namespace GiNaC;

  std::vector<double> res;

  for (unsigned int i = 0; i < list.nops(); i++)
    res.push_back(ex_to<numeric>(evalf(list[i])).to_double());

  return res;
}

// TODO: Turn the following method is a function.

/**
 * Compute the Euclidean norm of a vector
 *
 * @param[in] v vector to normalize
 * @returns norm of the given vector
 */
double Parallelotope::euclidNorm(const Vector &v) const
{

  double norm = 0;

  for (unsigned int i = 0; i < v.size(); i++) {
    norm = norm + (v[i] * v[i]);
  }

  return sqrt(norm);
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