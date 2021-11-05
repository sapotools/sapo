/**
 * @file LinearSystem.cpp
 * Represent and manipulate a linear system
 * It can be used to represent polytopes (reached states, parameters, etc.)
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#include <iostream>

#include <glpk.h>
#include <limits>

#include "LinearSystem.h"
#include "LinearSystemSet.h"

#define MAX_APPROX_ERROR 1e-8 // necessary for double comparison

using namespace std;
using namespace GiNaC;

std::ostream &operator<<(std::ostream &out, const LinearSystem &ls)
{
  for (unsigned int row_idx = 0; row_idx < ls.size(); row_idx++) {
    if (row_idx != 0) {
      out << std::endl;
    }

    for (unsigned int col_idx = 0; col_idx < ls.dim(); col_idx++) {
      out << ls.getA(row_idx, col_idx) << " ";
    }
    out << "<= " << ls.getb(row_idx);
  }

  return out;
}

/**
 * Optimize a linear system
 *
 * @param[in] A template matrix of the system to optimize
 * @param[in] b offset vector of the system to optimize
 * @param[in] obj_fun objective function
 * @param[in] min_max minimize of maximize Ax<=b (GLP_MIN=min, GLP_MAX=max)
 * @return optimum
 */
double solveLinearSystem(const vector<vector<double>> &A,
                         const vector<double> &b,
                         const vector<double> &obj_fun, const int min_max)
{
  unsigned int num_rows = A.size();
  unsigned int num_cols = obj_fun.size();
  unsigned int size_lp = num_rows * num_cols;

  int *ia, *ja;
  double *ar;

  ia = (int *)calloc(size_lp + 1, sizeof(int));
  ja = (int *)calloc(size_lp + 1, sizeof(int));
  ar = (double *)calloc(size_lp + 1, sizeof(double));

  glp_prob *lp;
  lp = glp_create_prob();
  glp_set_obj_dir(lp, min_max);

  // Turn off verbose mode
  glp_smcp lp_param;
  glp_init_smcp(&lp_param);
  lp_param.msg_lev = GLP_MSG_ERR;

  glp_add_rows(lp, num_rows);
  for (unsigned int i = 0; i < num_rows; i++) {
    glp_set_row_bnds(lp, i + 1, GLP_UP, 0.0, b[i]);
  }

  glp_add_cols(lp, num_cols);
  for (unsigned int i = 0; i < num_cols; i++) {
    glp_set_col_bnds(lp, i + 1, GLP_FR, 0.0, 0.0);
  }

  for (unsigned int i = 0; i < num_cols; i++) {
    glp_set_obj_coef(lp, i + 1, obj_fun[i]);
  }

  unsigned int k = 1;
  for (unsigned int i = 0; i < num_rows; i++) {
    for (unsigned int j = 0; j < num_cols; j++) {
      ia[k] = i + 1;
      ja[k] = j + 1;
      ar[k] = A[i][j]; /* a[i+1,j+1] = A[i][j] */
      k++;
    }
  }

  glp_load_matrix(lp, size_lp, ia, ja, ar);
  glp_exact(lp, &lp_param);
  //	glp_simplex(lp, &lp_param);

  double res;
  if (glp_get_status(lp) == GLP_UNBND) {
    if (min_max == GLP_MAX) {
      res = std::numeric_limits<double>::infinity();
    } else {
      res = -std::numeric_limits<double>::infinity();
    }
  } else {
    res = glp_get_obj_val(lp);
  }
  glp_delete_prob(lp);
  glp_free_env();
  free(ia);
  free(ja);
  free(ar);
  return res;
}

/**
 * Check if if a vector is null, i.e.,
 * it's a vector of zeros (used to detected useless constraints)
 *
 * @param[in] line vector to test
 * @return true is the vector is nulle
 */
bool zeroLine(const vector<double> &line)
{
  bool zeros = true;
  unsigned int i = 0;
  while (zeros && i < line.size()) {
    zeros = zeros && (abs(line[i]) < MAX_APPROX_ERROR);
    i++;
  }
  return zeros;
}

/**
 * Constructor that instantiates a linear system
 *
 * @param[in] A template matrix
 * @param[in] b offset vector
 */
LinearSystem::LinearSystem(const vector<vector<double>> &A,
                           const vector<double> &b)
{

  bool smart_insert = false;

  if (!smart_insert) {
    this->A = A;
    this->b = b;
  } else {
    for (unsigned int i = 0; i < A.size(); i++) {
      if (!this->isIn(A[i], b[i]) && (!zeroLine(A[i]))) {
        this->A.push_back(A[i]);
        this->b.push_back(b[i]);
      }
    }
  }
}

/**
 * Constructor that instantiates an empty linear system
 */
LinearSystem::LinearSystem(): A(), b() {}

/**
 * Copy constructor
 *
 * @param[in] orig the original linear system
 */
LinearSystem::LinearSystem(const LinearSystem &orig): A(orig.A), b(orig.b) {}

/**
 * Swap constructor
 *
 * @param[in] orig the original linear system
 */
LinearSystem::LinearSystem(LinearSystem &&orig)
{
  swap(*this, orig);
}

/**
 * Check if a constraint belongs to the linear system
 *
 * @param[in] Ai direction
 * @param[in] bi offset
 * @returns true is Ai x <= b is in the linear system
 */
bool LinearSystem::isIn(vector<double> Ai, const double bi) const
{
  Ai.push_back(bi);

  for (unsigned int i = 0; i < this->A.size(); i++) {
    vector<double> line = this->A[i];
    line.push_back(this->b[i]);
    bool is_in = true;
    for (unsigned int j = 0; j < Ai.size(); j++) {
      is_in = is_in && (abs(Ai[j] - line[j]) < MAX_APPROX_ERROR);
    }
    if (is_in)
      return true;
  }
  return false;
}

/**
 * Constructor that instantiates a linear system from a set of symbolic
 * expressions
 *
 * @param[in] vars list of variables appearing in the constraints
 * @param[in] constraints symbolic constraints
 */
LinearSystem::LinearSystem(const lst &vars, const lst &constraints)
{
  lst lconstraints = constraints;
  lconstraints.unique();

  for (lst::const_iterator c_it = begin(lconstraints);
       c_it != end(lconstraints); ++c_it) {
    vector<double> Ai;
    ex const_term(*c_it);

    for (lst::const_iterator v_it = begin(vars); v_it != end(vars); ++v_it) {
      // Extract the coefficient of the i-th variable (grade 1)
      double coeff = ex_to<numeric>(evalf(c_it->coeff(*v_it, 1))).to_double();
      Ai.push_back(coeff);

      // Project to obtain the constant term
      const_term = const_term.coeff(*v_it, 0);
    }

    double bi = ex_to<numeric>(evalf(const_term)).to_double();

    if (!this->isIn(Ai, -bi)) {
      this->A.push_back(Ai);
      this->b.push_back(-bi);
    }
  }
}

/**
 * Return the (i,j) element of the template matrix
 *
 * @param[in] i row index
 * @param[in] j column index
 * @return (i,j) element
 */
const double &LinearSystem::getA(unsigned int i, unsigned int j) const
{
  if (i < this->A.size() && j < this->A[j].size()) {
    return this->A[i][j];
  }
  std::cerr << "LinearSystem::getA : i and j must be within the LS->A size"
            << std::endl;
  exit(EXIT_FAILURE);
}

/**
 * Return the i-th element of the offset vector
 *
 * @param[in] i column index
 * @return i-th element
 */
const double &LinearSystem::getb(unsigned int i) const
{
  if (i < this->b.size()) {
    return this->b[i];
  }
  std::cerr << "LinearSystem::getb : i and j must be within the LS->b size"
            << std::endl;
  exit(EXIT_FAILURE);
}

/**
 * Determine whether this linear system is empty, i.e.,
 * the linear system has no solutions.
 *
 * Due to approximation errors, it may return false for some empty
 * systems too. However, when it returns true, the set is certainly empty.
 *
 * @param[in] strict_inequality specifies whether the linear system is a
 * 						strict inequality (i.e., Ax <
 * b).
 * @return a Boolean value. If the returned value is true, then the
 *       linear system is empty.
 */
bool LinearSystem::isEmpty(const bool strict_inequality) const
{
  if (this->size() == 0) {
    return false;
  }

  vector<vector<double>> extA(this->A);
  vector<double> obj_fun(this->A[0].size(), 0);
  obj_fun.push_back(1);

  // Add an extra variable to the linear system
  for (vector<vector<double>>::iterator row_it = begin(extA);
       row_it != end(extA); ++row_it) {
    row_it->push_back(-1);
  }

  const double z = solveLinearSystem(extA, this->b, obj_fun, GLP_MIN);

  return (z > MAX_APPROX_ERROR)
         || (strict_inequality && (z >= MAX_APPROX_ERROR));
}

/**
 * Check whether all the solutions of a linear system are also solutions for
 * another linear system.
 *
 * This method establishes whether all the solutions of a linear system
 * are are also solutions for another linear system. Due to approximation
 * errors, it may return false even if this is the case. However, whenever
 * it returns true, the all the solutions of the linear system are certainly
 * solutions for the linear system passed as parameter.
 *
 * @param[in] ls is the linear system whose set of solutions is compated this
 * 	that of this linear system.
 * @return a Boolean value. When some of the solutions of this linear system
 *     are not solutions for the parameter, the returned value is false. When
 *     the method returns true, all the solution of the object are also
 *     solutions for the parameter. There are cases in which the set of
 *     object solutions is a subset of the parameter solutions and, still,
 *     this method returns false.
 */
bool LinearSystem::satisfies(const LinearSystem &ls) const
{

  for (unsigned int i = 0; i < ls.size(); i++) {
    if (!this->satisfies(ls.A[i], ls.b[i])) {
      return false;
    }
  }

  return true;
}

template<typename T>
bool are_independent(const std::vector<T> &v1, const std::vector<T> &v2)
{
  if (v1.size() != v2.size()) {
    return true;
  }

  if (v1.size() == 0) {
    return false;
  }

  unsigned int fnz_v1(0);
  while (fnz_v1 < v1.size() && v1[fnz_v1] == 0)
    fnz_v1++;

  unsigned int fnz_v2(0);
  while (fnz_v2 < v2.size() && v2[fnz_v2] == 0)
    fnz_v2++;

  if (fnz_v1 != fnz_v2) {
    return true;
  }

  for (unsigned int i = 1; i < v1.size(); ++i) {
    if (v1[fnz_v1] * v2[i]
        != v2[fnz_v2] * v1[i]) { // this is to avoid numerical errors,
                                 // but it can produce overflows
      return true;
    }
  }

  return false;
}

std::vector<unsigned int>
get_a_linear_system_base(const std::vector<std::vector<double>> &A)
{
  if (A.size() == 0) {
    return std::vector<unsigned int>();
  }

  std::vector<unsigned int> base{0};

  unsigned int row_idx = 1;
  while (row_idx < A.size()) {
    bool indep_from_base = true;
    auto b_it = std::begin(base);
    while (indep_from_base && b_it != std::end(base)) {
      indep_from_base = are_independent(A[row_idx], A[*b_it]);
      ++b_it;
    }

    if (indep_from_base) {
      base.push_back(row_idx);
    }
    ++row_idx;
  }

  return base;
}

inline std::vector<bool>
get_a_ls_base_bit_vector(const std::vector<std::vector<double>> &A)
{
  std::vector<unsigned int> base = get_a_linear_system_base(A);

  std::vector<bool> bvect_base(A.size(), false);

  for (auto it = std::begin(base); it != std::end(base); ++it) {
    bvect_base[*it] = true;
  }

  return bvect_base;
}

std::list<LinearSystem> LinearSystem::get_a_finer_covering(
    const std::vector<bool> &bvect_base, const unsigned int cidx,
    std::list<LinearSystem> &tmp_covering, std::vector<std::vector<double>> &A,
    std::vector<double> &b) const
{
  if (this->A.size() == cidx) {
    LinearSystem ls(A, b);
    ls.simplify();

    tmp_covering.push_back(ls);
    return tmp_covering;
  }

  if (bvect_base[cidx]) {
    A.push_back(this->A[cidx]);
    b.push_back(this->b[cidx]);

    try {
      const double min_value = minLinearSystem(this->A[cidx]);
      const double avg_value = (this->b[cidx] + min_value) / 2;

      A.push_back(get_complementary(this->A[cidx]));
      b.push_back(-avg_value);

      get_a_finer_covering(bvect_base, cidx + 1, tmp_covering, A, b);

      b[b.size() - 1] = -min_value;
      b[b.size() - 2] = avg_value;

      get_a_finer_covering(bvect_base, cidx + 1, tmp_covering, A, b);

      A.pop_back();
      b.pop_back();

    } catch (std::logic_error &e) {
      std::cerr << "The linear system solutions are not a closed polyheadron."
                << std::endl;

      get_a_finer_covering(bvect_base, cidx + 1, tmp_covering, A, b);
    }

    A.pop_back();
    b.pop_back();
  } else {
    get_a_finer_covering(bvect_base, cidx + 1, tmp_covering, A, b);
  }

  return tmp_covering;
}

std::list<LinearSystem> LinearSystem::get_a_finer_covering() const
{
  std::list<LinearSystem> result;

  std::vector<bool> bvect_base = get_a_ls_base_bit_vector(this->A);

  std::vector<std::vector<double>> A;
  std::vector<double> b;

  get_a_finer_covering(bvect_base, 0, result, A, b);

  return result;
}

/**
 * Minimize the linear system
 *
 * @param[in] vars list of variables
 * @param[in] obj_fun objective function
 * @return minimum
 */
double LinearSystem::minLinearSystem(const lst &vars, const ex &obj_fun) const
{

  vector<double> obj_fun_coeffs;
  ex const_term(obj_fun);

  // Extract the coefficient of the i-th variable (grade 1)
  for (lst::const_iterator v_it = begin(vars); v_it != end(vars); ++v_it) {
    double coeff = ex_to<numeric>(evalf(obj_fun.coeff(*v_it, 1))).to_double();

    obj_fun_coeffs.push_back(coeff);
    const_term = const_term.coeff(*v_it, 0);
  }

  const double c = ex_to<numeric>(evalf(const_term)).to_double();
  const double min
      = solveLinearSystem(this->A, this->b, obj_fun_coeffs, GLP_MIN);

  return (min + c);
}

/**
 * Minimize the linear system
 *
 * @param[in] obj_fun objective function
 * @return minimum
 */
double
LinearSystem::minLinearSystem(const vector<double> &obj_fun_coeffs) const
{
  return solveLinearSystem(this->A, this->b, obj_fun_coeffs, GLP_MIN);
}

/**
 * Maximize the linear system
 *
 * @param[in] obj_fun objective function
 * @return maximum
 */
double
LinearSystem::maxLinearSystem(const vector<double> &obj_fun_coeffs) const
{
  return solveLinearSystem(this->A, this->b, obj_fun_coeffs, GLP_MAX);
}

/**
 * Maximize the linear system
 *
 * @param[in] vars list of variables
 * @param[in] obj_fun objective function
 * @return maximum
 */
double LinearSystem::maxLinearSystem(const lst &vars, const ex &obj_fun) const
{

  vector<double> obj_fun_coeffs;
  ex const_term(obj_fun);

  // Extract the coefficient of the i-th variable (grade 1)
  for (lst::const_iterator v_it = begin(vars); v_it != end(vars); ++v_it) {
    double coeff = ex_to<numeric>(evalf(obj_fun.coeff(*v_it, 1))).to_double();
    obj_fun_coeffs.push_back(coeff);
    const_term = const_term.coeff(*v_it, 0);
  }

  const double c = ex_to<numeric>(evalf(const_term)).to_double();
  const double max
      = solveLinearSystem(this->A, this->b, obj_fun_coeffs, GLP_MAX);

  return (max + c);
}

/**
 * Create a new linear system by joining the constraints of two linear system
 *
 * @param[in] ls a linear system
 * @return linear system obtained by joining the constraints of this object and
 *      those of of the parameter.
 */
LinearSystem &LinearSystem::intersectWith(const LinearSystem &ls)
{
  for (unsigned int i = 0; i < ls.size(); i++) {
    if (!this->satisfies(ls.A[i], ls.b[i])) { // check for duplicates
      (this->A).push_back(ls.A[i]);
      (this->b).push_back(ls.b[i]);
    }
  }

  return *this;
}

LinearSystem intersection(const LinearSystem &A, const LinearSystem &B)
{
  LinearSystem result(A.A, A.b);

  result.intersectWith(B);

  return result;
}

/**
 * Check whether all the solutions of a linear system satisfy a constraint.
 *
 * This method establishes whether all the solutions of a linear system
 * satisfy a constraint. Due to approximation errors, it may return
 * false even if this is the case. However, whenever it returns true, the
 * all the solutions of the linear system certainly satisfy the inequality.
 *
 * @param[in] i is the index of the constraint to be checked
 * @return a Boolean value. When some of the solutions of the linear system
 *     do not satisfy the inequality, the returned value is false. When the
 *     method returns true, the constraint is certainly satisfied by any of
 *     the solutions of the system. There are cases in which the constraint
 *     is satisfied by all the solutions and this method returns false.
 */
bool LinearSystem::satisfies(const std::vector<double> &Ai,
                             const double bi) const
{
  if (size() == 0)
    return false;

  if (isIn(Ai, bi)) {
    return true;
  }

  double max = this->maxLinearSystem(Ai);
  if (max + MAX_APPROX_ERROR
      <= bi) { /* This should be max <= bi,
                                  however, due to double approximation
                                  errors, testing whether the distance
                                  between max and bi is greater than a
                                  fixed positive approximation constant
                                  is more conservative */
    return true;
  }

  return false;

  auto max_coeff = max_element(std::begin(Ai), std::end(Ai));
  auto min_coeff = min_element(std::begin(Ai), std::end(Ai));
  return ((*max_coeff == *min_coeff) && (*min_coeff == 0) && (bi >= 0));
}

/**
 * Check whether one of the constraints in a linear system is redundant.
 *
 * This method establishes whether the i-th constraint of a linear
 * system is redundant. Due to approximation errors, it may return
 * false even if the constraint is redundant. However, whenever it
 * returns true, the constraint is certainly redundant.
 *
 * @param[in] i is the index of the constraint to be checked
 * @return a Boolean value. When the constraint is non-redundanct, the
 *     returned value is true. When it returns true, the constraint is
 *     certainly redundant. There are cases in which the constraint is
 *     redundant and this method returns false.
 */
bool LinearSystem::constraintIsRedundant(const unsigned int i) const
{
  LinearSystem tmp(*this);
  std::vector<double> Ai(dim(), 0);
  double bi(0);

  // replace the i-th constraint with the empty constraint
  std::swap(Ai, tmp.A[i]);
  std::swap(bi, tmp.b[i]);

  // check whether the i-th constraint is redundant
  bool result(tmp.satisfies(Ai, bi));

  // reinsert the i-th into the system
  std::swap(Ai, tmp.A[i]);
  std::swap(bi, tmp.b[i]);

  return result;
}

/**
 * Remove redundant constraints from a linear system.
 *
 * This method removes redundant constraints from the system.
 * The order of the non-redundant constraints can be shuffled after
 * the call.
 *
 * @return A reference to this object after removing all the
 *         redundant constraints.
 */
LinearSystem &LinearSystem::simplify()
{
  unsigned int i = 0, last_non_redundant = size() - 1;

  while (i < last_non_redundant) { // for every unchecked constraint

    // if it is redundant
    if (constraintIsRedundant(i)) {
      // swap it with the last non-reduntant constraint
      swap(A[i], A[last_non_redundant]);
      swap(b[i], b[last_non_redundant]);

      // decrease the number of the non-reduntant constraints
      last_non_redundant--;
    } else { // otherwise, i.e. if it is not redundant

      // increase the index of the next constraint to be checked
      i++;
    }
  }

  // if the last constraint to be checked is redundant
  if (constraintIsRedundant(last_non_redundant)) {
    // reduce the number of non-reduntant constraints
    last_non_redundant--;
  }

  // remove the redundant constraints that are at the end of the system
  A.resize(last_non_redundant + 1);
  b.resize(last_non_redundant + 1);

  return *this;
}

/**
 * Remove redundant constraints
 */
LinearSystem LinearSystem::get_simplified() const
{
  LinearSystem simplier(*this);

  return simplier.simplify();
}

/**
 * Determine the volume of the bounding box of the linear system
 *
 * @return volume of the bounding box
 */
double LinearSystem::volBoundingBox() const
{

  vector<double> zeros(this->dim(), 0);
  double vol = 1;

  for (unsigned int i = 0; i < this->dim(); i++) {
    vector<double> facet = zeros;
    facet[i] = 1;
    const double b_plus = solveLinearSystem(this->A, this->b, facet, GLP_MAX);
    facet[i] = -1;
    const double b_minus = solveLinearSystem(this->A, this->b, facet, GLP_MAX);
    vol = vol * (b_plus + b_minus);
  }

  return vol;
}

/**
 * Print the linear system in Matlab format (for plotregion script)
 *
 * @param[in] os is the output stream
 * @param[in] color color of the polytope to plot
 */
void LinearSystem::plotRegion(std::ostream &os, const char color) const
{

  if (this->dim() > 3) {
    std::cerr << "LinearSystem::plotRegion : maximum 3d sets are allowed"
              << std::endl;
    exit(EXIT_FAILURE);
  }

  os << "Ab = [" << std::endl;
  for (unsigned i = 0; i < A.size(); i++) {
    for (auto el = std::begin(A[i]); el != std::end(A[i]); ++el) {
      os << *el << " ";
    }
    os << " " << this->b[i] << ";" << std::endl;
  }
  os << "];" << std::endl;
  os << "plotregion(-Ab(:,1:" << this->A[0].size() << "),-Ab(:,"
     << this->A[0].size() + 1 << "),[],[],";

  if (color == ' ') {
    os << "color";
  } else {
    os << color;
  }
  os << ");" << std::endl;
}

/**
 * Print the 2d linear system in Matlab format (for plotregion script) over
 * time
 *
 * @param[in] os is the output stream
 * @param[in] t thickness of the set to plot
 */
void LinearSystem::plotRegionT(std::ostream &os, const double t) const
{
  if (this->dim() > 2) {
    std::cerr << "LinearSystem::plotRegionT : maximum 2d sets are allowed"
              << std::endl;
    exit(EXIT_FAILURE);
  }

  os << "Ab = [" << std::endl;
  os << " 1 ";
  for (unsigned int j = 0; j < this->A[0].size(); j++) {
    os << " 0 ";
  }
  os << t << ";" << std::endl;
  os << " -1 ";
  for (unsigned int j = 0; j < this->A[0].size(); j++) {
    os << " 0 ";
  }
  os << -t << ";" << std::endl;

  for (unsigned i = 0; i < A.size(); i++) {
    os << " 0 ";
    for (auto el = std::begin(A[i]); el != std::end(A[i]); ++el) {
      os << *el << " ";
    }
    os << this->b[i] << ";" << std::endl;
  }

  os << "];" << std::endl;
  os << "plotregion(-Ab(:,1:3),-Ab(:,4),[],[],color);" << std::endl;
}

/**
 * Print the specified projections in Matlab format (for plotregion script)
 * into a file
 *
 * @param[in] os is the output stream
 * @param[in] rows rows to be plot
 * @param[in] cols colors of the plots
 */
void LinearSystem::plotRegion(std::ostream &os, const vector<int> &rows,
                              const vector<int> &cols) const
{

  if (cols.size() > 3) {
    std::cerr << "LinearSystem::plotRegion : cols maximum 3d sets are allowed"
              << std::endl;
    exit(EXIT_FAILURE);
  }

  os << "Ab = [" << std::endl;
  for (auto r_it = std::begin(rows); r_it != std::end(rows); ++r_it) {
    for (auto c_it = std::begin(cols); c_it != std::end(cols); ++c_it) {
      os << this->A[*r_it][*c_it] << " ";
    }
    os << " " << this->b[*r_it] << ";" << std::endl;
  }
  os << "];" << std::endl;
  os << "plotregion(-Ab(:,1:" << cols.size() << "),-Ab(:," << cols.size() + 1
     << "),[],[],color);" << std::endl;
}
