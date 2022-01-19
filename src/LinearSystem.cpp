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
#include <cmath>

#include "LinearSystem.h"
#include "SymbolicAlgebra.h"

#define MAX_APPROX_ERROR 1e-8 // necessary for double comparison

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

JSON::ostream &operator<<(JSON::ostream &out, const LinearSystem &ls)
{
  out << "{\"A\":" << ls.getA() << ","
      << "\"b\":" << ls.getb() << "}";

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
double optimize(const std::vector<std::vector<double>> &A,
                const std::vector<double> &b,
                const std::vector<double> &obj_fun, const int min_max)
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
 * Optimize a linear system
 *
 * @param[in] obj_fun objective function
 * @param[in] min_max minimize of maximize Ax<=b (GLP_MIN=min, GLP_MAX=max)
 * @return optimum
 */
double LinearSystem::optimize(const std::vector<double> &obj_fun,
                              const int min_max) const
{
  return ::optimize(this->A, this->b, obj_fun, min_max);
}

/**
 * Check if if a vector is null, i.e.,
 * it's a vector of zeros (used to detected useless constraints)
 *
 * @param[in] line vector to test
 * @return true is the vector is nulle
 */
bool zeroLine(const std::vector<double> &line)
{
  bool zeros = true;
  unsigned int i = 0;
  while (zeros && i < line.size()) {
    zeros = zeros && (std::abs(line[i]) < MAX_APPROX_ERROR);
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
LinearSystem::LinearSystem(const std::vector<std::vector<double>> &A,
                           const std::vector<double> &b)
{

  bool smart_insert = false;

  if (!smart_insert) {
    this->A = A;
    this->b = b;
  } else {
    for (unsigned int i = 0; i < A.size(); i++) {
      if (!this->is_in(A[i], b[i]) && (!zeroLine(A[i]))) {
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
bool LinearSystem::is_in(std::vector<double> Ai, const double bi) const
{
  Ai.push_back(bi);

  for (unsigned int i = 0; i < this->A.size(); i++) {
    std::vector<double> line = this->A[i];
    line.push_back(this->b[i]);
    bool is_in = true;
    for (unsigned int j = 0; j < Ai.size(); j++) {
      is_in = is_in && (std::abs(Ai[j] - line[j]) < MAX_APPROX_ERROR);
    }
    if (is_in)
      return true;
  }
  return false;
}

// TODO: remove this method as it uses non-necessary symbolic expressions.
/**
 * Constructor that instantiates a linear system from a set of symbolic
 * expressions
 *
 * @param[in] vars list of variables appearing in the constraints
 * @param[in] constraints symbolic constraints
 */
LinearSystem::LinearSystem(
    const std::vector<SymbolicAlgebra::Symbol<>> &vars,
    const std::vector<SymbolicAlgebra::Expression<>> &constraints)
{
  using namespace SymbolicAlgebra;
  std::vector<Expression<>> lconstraints = constraints;
  // lconstraints.unique();  // remove multiple copies of the same expression

  for (auto c_it = begin(lconstraints); c_it != end(lconstraints); ++c_it) {
    std::vector<double> Ai;
    Expression<> const_term(*c_it);

    for (auto v_it = begin(vars); v_it != end(vars); ++v_it) {
      // Extract the coefficient of the i-th variable (grade 1)
      double coeff = (c_it->get_coeff(*v_it, 1)).evaluate<double>();
      Ai.push_back(coeff);

      // Project to obtain the constant term
      const_term = const_term.get_coeff(*v_it, 0);
    }

    double bi = const_term.evaluate<double>();

    if (!this->is_in(Ai, -bi)) {
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
 * Establish whether a linear system has solutions
 *
 * Due to approximation errors, it may return true for some systems
 * having no solution too. However, when it returns false, the linear
 * system certainly has no solution.
 *
 * @param[in] strict_inequality specifies whether the linear system is
 *         a strict inequality (i.e., Ax < b).
 * @return a Boolean value. If the returned value is false, then the
 *       linear system has no solution.
 */
bool LinearSystem::has_solutions(const bool strict_inequality) const
{
  if (this->size() == 0) {
    return true;
  }

  std::vector<std::vector<double>> extA(this->A);
  std::vector<double> obj_fun(this->A[0].size(), 0);
  obj_fun.push_back(1);

  // Add an extra variable to the linear system
  for (std::vector<std::vector<double>>::iterator row_it = begin(extA);
       row_it != end(extA); ++row_it) {
    row_it->push_back(-1);
  }

  const double z = ::optimize(extA, this->b, obj_fun, GLP_MIN);

  return (z <= MAX_APPROX_ERROR)
         || (!strict_inequality && (z < MAX_APPROX_ERROR));
}

/**
 * Minimize the linear system
 *
 * @param[in] vars list of variables
 * @param[in] obj_fun objective function
 * @return minimum
 */
double
LinearSystem::minimize(const std::vector<SymbolicAlgebra::Symbol<>> &vars,
                       const SymbolicAlgebra::Expression<> &obj_fun) const
{
  using namespace SymbolicAlgebra;

  std::vector<double> obj_fun_coeffs;
  Expression<> const_term(obj_fun);

  // Extract the coefficient of the i-th variable (grade 1)
  for (auto v_it = begin(vars); v_it != end(vars); ++v_it) {
    double coeff = (obj_fun.get_coeff(*v_it, 1)).evaluate<double>();

    obj_fun_coeffs.push_back(coeff);
    const_term = const_term.get_coeff(*v_it, 0);
  }

  const double c = const_term.evaluate<double>();
  const double min = optimize(obj_fun_coeffs, GLP_MIN);

  return (min + c);
}

/**
 * Minimize the linear system
 *
 * @param[in] obj_fun objective function
 * @return minimum
 */
double LinearSystem::minimize(const std::vector<double> &obj_fun_coeffs) const
{
  return optimize(obj_fun_coeffs, GLP_MIN);
}

/**
 * Maximize the linear system
 *
 * @param[in] obj_fun objective function
 * @return maximum
 */
double LinearSystem::maximize(const std::vector<double> &obj_fun_coeffs) const
{
  return optimize(obj_fun_coeffs, GLP_MAX);
}

/**
 * Maximize the linear system
 *
 * @param[in] vars list of variables
 * @param[in] obj_fun objective function
 * @return maximum
 */
double
LinearSystem::maximize(const std::vector<SymbolicAlgebra::Symbol<>> &vars,
                       const SymbolicAlgebra::Expression<> &obj_fun) const
{
  using namespace SymbolicAlgebra;

  std::vector<double> obj_fun_coeffs;
  Expression<> const_term(obj_fun);

  // Extract the coefficient of the i-th variable (grade 1)
  for (auto v_it = begin(vars); v_it != end(vars); ++v_it) {
    const double coeff = obj_fun.get_coeff(*v_it, 1).evaluate<double>();
    obj_fun_coeffs.push_back(coeff);
    const_term = const_term.get_coeff(*v_it, 0);
  }

  const double c = const_term.evaluate<double>();

  return maximize(obj_fun_coeffs) + c;
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

  if (is_in(Ai, bi)) {
    return true;
  }

  double max = this->maximize(Ai);
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
bool LinearSystem::constraint_is_redundant(const unsigned int i) const
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
  if (size() == 0) {
    return *this;
  }

  unsigned int i = 0, last_non_redundant = size() - 1;

  while (i < last_non_redundant) { // for every unchecked constraint

    // if it is redundant
    if (constraint_is_redundant(i)) {
      // swap it with the last non-reduntant constraint
      std::swap(A[i], A[last_non_redundant]);
      std::swap(b[i], b[last_non_redundant]);

      // decrease the number of the non-reduntant constraints
      last_non_redundant--;
    } else { // otherwise, i.e. if it is not redundant

      // increase the index of the next constraint to be checked
      i++;
    }
  }

  // if the last constraint to be checked is redundant
  if (constraint_is_redundant(last_non_redundant)) {
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
