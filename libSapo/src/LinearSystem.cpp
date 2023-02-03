/**
 * @file LinearSystem.cpp
 * @author Alberto Casagrande <acasagrande@units.it>
 * @brief Linear algebra classes and functions code
 * @version 0.1
 * @date 2021-11-20
 *
 * @copyright Copyright (c) 2021-2022
 */

#include <iostream>

#include <limits>
#include <cmath>
#include <sstream>

#include "LinearSystem.h"
#include "SymbolicAlgebra.h"
#include "LinearAlgebraIO.h"
#include "ErrorHandling.h"

/**
 * @brief Print a linear system in a stream
 *
 * @param out is the output stream
 * @param ls is the linear system to be print
 * @return a reference to the output stream
 */
std::ostream &operator<<(std::ostream &out, const LinearSystem &ls)
{
  for (unsigned int row_idx = 0; row_idx < ls.size(); row_idx++) {
    if (row_idx != 0) {
      out << std::endl;
    }

    for (unsigned int col_idx = 0; col_idx < ls.dim(); col_idx++) {
      out << ls.A(row_idx)[col_idx] << " ";
    }
    out << "<= " << ls.b(row_idx);
  }

  return out;
}

/// @private
template<typename T>
unsigned int find_in(const std::vector<LinearAlgebra::Vector<T>> &A,
                     const LinearAlgebra::Vector<T> &v)
{
  for (unsigned int i = 0; i < A.size(); ++i) {
    if (A[i] == v) {
      return i;
    }
  }

  return A.size();
}

OptimizationResult<double>
LinearSystem::optimize(const LinearAlgebra::Vector<double> &obj_fun,
                       const bool maximize) const
{
  SimplexMethodOptimizer optimizer;

  return optimizer(this->_A, this->_b, obj_fun,
                   (maximize ? MAXIMIZE : MINIMIZE));
}

LinearSystem &
LinearSystem::add_constraint(const LinearAlgebra::Vector<double> &v,
                             const double &b)
{
  if (dim() != v.size()) {
    SAPO_ERROR("the new constraint and the linear system "
               "differ in the number of variables",
               std::domain_error);
  }

  _A.push_back(v);
  _b.push_back(b);

  return *this;
}

LinearSystem &LinearSystem::operator=(const LinearSystem &orig)
{
  _A = orig._A;
  _b = orig._b;

  return *this;
}

bool have_disjoint_solutions(const LinearSystem &ls1, const LinearSystem &ls2)
{
  if (ls1.dim() != ls2.dim()) {
    SAPO_ERROR("the two linear systems must have the "
               "same number of variables",
               std::domain_error);
  }

  if (ls1.dim() == 0) {
    return false;
  }

  LinearAlgebra::Vector<double> obj(ls1.dim(), 1);

  auto A{ls1.A()};
  A.reserve(ls1.size() + ls2.size());
  std::copy(std::begin(ls2.A()), std::end(ls2.A()), std::back_inserter(A));

  auto b{ls1.b()};
  b.reserve(ls1.size() + ls2.size());
  std::copy(std::begin(ls2.b()), std::end(ls2.b()), std::back_inserter(b));

  SimplexMethodOptimizer optimizer;

  auto result = optimizer(A, b, obj);

  return result.status() == result.INFEASIBLE;
}

/**
 * Minimize the linear system
 *
 * @param[in] obj_fun objective function
 * @return minimum
 */
OptimizationResult<double>
LinearSystem::minimize(const LinearAlgebra::Vector<double> &obj_fun) const
{
  SimplexMethodOptimizer optimizer;

  return optimizer(this->_A, this->_b, obj_fun);
}

/**
 * Maximize the linear system
 *
 * @param[in] obj_fun objective function
 * @return maximum
 */
OptimizationResult<double>
LinearSystem::maximize(const LinearAlgebra::Vector<double> &obj_fun) const
{
  SimplexMethodOptimizer optimizer;

  return optimizer(this->_A, this->_b, obj_fun, MAXIMIZE);
}

/**
 * @brief Get the ordinal of a number
 *
 * @param number is the number whose ordinal must be returned
 * @return the ordinal of `number`
 */
std::string get_ordinal(unsigned int number)
{
  std::ostringstream oss;

  oss << number;

  number = number % 100;

  if (number / 10 == 1) {
    oss << "th";

    return oss.str();
  }

  switch (number % 10) {
  case 1:
    oss << "st";
    break;
  case 2:
    oss << "nd";
    break;
  case 3:
    oss << "rd";
    break;
  default:
    oss << "th";
  }

  return oss.str();
}

/**
 * Constructor
 *
 * @param[in] A is the matrix
 * @param[in] b is the offset vector
 */
LinearSystem::LinearSystem(const std::vector<LinearAlgebra::Vector<double>> &A,
                           const LinearAlgebra::Vector<double> &b)
{
  if (A.size() != b.size()) {
    SAPO_ERROR("the number of rows in A must equals the "
               "number of elements in b",
               std::domain_error);
  }

  for (unsigned int i = 0; i < A.size(); i++) {
    if (A[i].size() != A[0].size()) {
      std::ostringstream oss;
      oss << "the " << get_ordinal(i + 1) << " row " << A[i]
          << " and the 1st one " << A[0] << " differ in "
          << "size" << std::endl;

      SAPO_ERROR(oss.str(), std::domain_error);
    }
  }

  bool smart_insert = false;

  if (!smart_insert) {
    this->_A = A;
    this->_b = b;
  } else {
    for (unsigned int i = 0; i < A.size(); i++) {
      if ((LinearAlgebra::norm_infinity(A[i]) > 0)
          && !this->satisfies(A[i], b[i])) {
        this->_A.push_back(A[i]);
        this->_b.push_back(b[i]);
      }
    }
  }
}

/**
 * Move constructor
 *
 * @param[in] A is the matrix
 * @param[in] b is the offset vector
 */
LinearSystem::LinearSystem(std::vector<LinearAlgebra::Vector<double>> &&A,
                           LinearAlgebra::Vector<double> &&b):
    _A(std::move(A)),
    _b(std::move(b))
{
  if (_A.size() != _b.size()) {
    SAPO_ERROR("the number of rows in A must equals the "
               "number of elements in b",
               std::domain_error);
  }

  for (unsigned int i = 0; i < _A.size(); i++) {
    if (_A[i].size() != _A[0].size()) {
      std::ostringstream oss;
      oss << "The " << get_ordinal(i + 1) << " row " << A[i]
          << " and the 1st one " << A[0] << " differ in "
          << "size" << std::endl;

      SAPO_ERROR(oss.str(), std::domain_error);
    }
  }
}

/**
 * Constructor that instantiates an empty linear system
 */
LinearSystem::LinearSystem(): _A(), _b() {}

/**
 * Copy constructor
 *
 * @param[in] ls is the original linear system
 */
LinearSystem::LinearSystem(const LinearSystem &ls): _A(ls._A), _b(ls._b) {}

/**
 * Swap constructor
 *
 * @param[in] ls is the linear system
 */
LinearSystem::LinearSystem(LinearSystem &&ls)
{
  swap(*this, ls);
}

/**
 * @brief Check whether two linear equations are the same
 *
 * @param A1 non-constant coefficient of the first equation
 * @param b1 constant coefficient of the first equation
 * @param A2 non-constant coefficient of the second equation
 * @param b2 constant coefficient of the second equation
 * @return whether the two linear equations are the same
 */
bool same_constraint(const LinearAlgebra::Vector<double> &A1, const double &b1,
                     const LinearAlgebra::Vector<double> &A2, const double &b2)
{
  if (A1.size() != A2.size()) {
    return false;
  }

  if (b1 != b2) {
    return false;
  }

  for (unsigned int i = 0; i < A1.size(); ++i) {
    if (A1[i] != A2[i]) {
      return false;
    }
  }

  return true;
}

bool LinearSystem::contains(const LinearAlgebra::Vector<double> &Ai,
                            const double &bi) const
{
  for (unsigned int i = 0; i < this->_A.size(); i++) {
    if (same_constraint(this->_A[i], this->_b[i], Ai, bi)) {
      return true;
    }
  }
  return false;
}

LinearSystem::LinearSystem(
    const std::vector<SymbolicAlgebra::Symbol<>> &x,
    const std::vector<SymbolicAlgebra::Expression<>> &expressions)
{
  using namespace SymbolicAlgebra;
  std::vector<Expression<>> lexpressions = expressions;
  // lconstraints.unique();  // remove multiple copies of the same expression

  for (auto e_it = begin(lexpressions); e_it != end(lexpressions); ++e_it) {
    LinearAlgebra::Vector<double> Ai;
    Expression<> const_term(*e_it);

    for (auto x_it = begin(x); x_it != end(x); ++x_it) {
      if (e_it->degree(*x_it) > 1) {
        SAPO_ERROR("one of the expressions is non-linear", std::domain_error);
      }

      // Extract the coefficient of the i-th variable (grade 1)
      double coeff = (e_it->get_coeff(*x_it, 1)).evaluate();
      Ai.push_back(coeff);

      // Project to obtain the constant term
      const_term = const_term.get_coeff(*x_it, 0);
    }

    double bi = const_term.evaluate();

    if (!this->contains(Ai, -bi)) {
      this->_A.push_back(Ai);
      this->_b.push_back(-bi);
    }
  }
}

/**
 * Return the i-th row of the matrix
 *
 * @param[in] i row index
 * @return the i-th row
 */
const LinearAlgebra::Vector<double> &
LinearSystem::A(const unsigned int i) const
{
  if (i >= this->_A.size()) {
    SAPO_ERROR("the parameter must be a valid row index "
               "for the linear system matrix A",
               std::domain_error);
  }
  return this->_A[i];
}

/**
 * Return the i-th element of the offset vector
 *
 * @param[in] i column index
 * @return i-th element
 */
const double &LinearSystem::b(unsigned int i) const
{
  if (i >= this->_b.size()) {
    SAPO_ERROR("the parameter must be a valid index for the linear "
               "system coefficient vector b",
               std::domain_error);
  }
  return this->_b[i];
}

bool LinearSystem::has_solutions(const bool strict_inequality) const
{
  if (this->size() == 0) {
    return true;
  }

  if (!strict_inequality) {
    OptimizationResult<double> res = maximize(_A[0]);

    return res.status() != res.INFEASIBLE;
  }

  for (auto row_it = std::begin(_A); row_it != std::end(_A); ++row_it) {
    OptimizationResult<double> res = maximize(*row_it);
    if (res.status() == res.INFEASIBLE) {
      return false;
    }
    OptimizationResult<double> res2 = minimize(*row_it);
    if (res.status() == res.INFEASIBLE) {
      return false;
    }

    if (res.objective_value() == res2.objective_value()) {
      return false;
    }
  }

  return true;
}

/**
 * Minimize the linear system
 *
 * @param[in] symbols array of symbols
 * @param[in] obj_fun objective function
 * @return minimum
 */
OptimizationResult<double>
LinearSystem::minimize(const std::vector<SymbolicAlgebra::Symbol<>> &symbols,
                       SymbolicAlgebra::Expression<> obj_fun) const
{
  using namespace SymbolicAlgebra;

  obj_fun.expand();

  LinearAlgebra::Vector<double> obj_fun_coeffs;
  Expression<> const_term(obj_fun);

  // Extract the coefficient of the i-th variable (grade 1)
  for (auto s_it = begin(symbols); s_it != end(symbols); ++s_it) {
    if (obj_fun.degree(*s_it) > 1) {
      SAPO_ERROR("the objective must be linear", std::domain_error);
    }
    double coeff = (obj_fun.get_coeff(*s_it, 1)).evaluate();

    obj_fun_coeffs.push_back(coeff);
    const_term = const_term.get_coeff(*s_it, 0);
  }

  const double c = const_term.evaluate();
  auto res = minimize(obj_fun_coeffs);

  if (res.status() != res.OPTIMUM_AVAILABLE) {
    return res;
  }

  return OptimizationResult<double>(res.optimum(), res.objective_value() + c);
}

/**
 * Maximize the linear system
 *
 * @param[in] symbols array of symbols
 * @param[in] obj_fun objective function
 * @return maximum
 */
OptimizationResult<double>
LinearSystem::maximize(const std::vector<SymbolicAlgebra::Symbol<>> &symbols,
                       SymbolicAlgebra::Expression<> obj_fun) const
{
  using namespace SymbolicAlgebra;

  obj_fun.expand();

  LinearAlgebra::Vector<double> obj_fun_coeffs;
  Expression<> const_term(obj_fun);

  // Extract the coefficient of the i-th variable (grade 1)
  for (auto s_it = begin(symbols); s_it != end(symbols); ++s_it) {
    if (obj_fun.degree(*s_it) > 1) {
      SAPO_ERROR("the objective must be linear", std::domain_error);
    }
    const double coeff = obj_fun.get_coeff(*s_it, 1).evaluate();
    obj_fun_coeffs.push_back(coeff);
    const_term = const_term.get_coeff(*s_it, 0);
  }

  const double c = const_term.evaluate();
  auto res = maximize(obj_fun_coeffs);

  if (res.status() != res.OPTIMUM_AVAILABLE) {
    return res;
  }

  return OptimizationResult<double>(res.optimum(), res.objective_value() + c);
}

/**
 * @brief Check whether the solutions of a linear system satisfy a constraint
 *
 * This method establishes whether all the solutions \f$s\f$ of the linear
 * system satisfy a constraint \f$\textrm{Ai} \cdot s \leq \textrm{bi}\f$.
 *
 * @param[in] Ai is a value vector
 * @param[in] bi is a scalar value
 * @return `true` if and only if \f$\textrm{Ai} \cdot s \leq \textrm{bi}\f$
 *         for any the solution \f$s\f$ of the linear system
 */
bool LinearSystem::satisfies(const LinearAlgebra::Vector<double> &Ai,
                             const double &bi) const
{
  if (size() == 0) {
    return false;
  }

  OptimizationResult<double> res = this->maximize(Ai);
  if ((res.status() == res.OPTIMUM_AVAILABLE) && res.objective_value() <= bi) {
    return true;
  }

  return false;
}

bool LinearSystem::satisfies(const LinearSystem &ls) const
{
  if (!this->has_solutions()) {
    return true;
  }

  for (unsigned int i = 0; i < ls.size(); i++) {
    if (!this->satisfies(ls._A[i], ls._b[i])) {
      return false;
    }
  }

  return true;
}

bool LinearSystem::constraint_is_redundant(const unsigned int i) const
{
  LinearSystem tmp(*this);
  LinearAlgebra::Vector<double> Ai(dim(), 0);
  double bi(0);

  // replace the i-th constraint with the empty constraint
  std::swap(Ai, tmp._A[i]);
  std::swap(bi, tmp._b[i]);

  // check whether the i-th constraint is redundant
  bool result(tmp.satisfies(Ai, bi));

  // reinsert the i-th into the system
  std::swap(Ai, tmp._A[i]);
  std::swap(bi, tmp._b[i]);

  return result;
}

/**
 * Remove redundant constraints from a linear system
 *
 * This method removes redundant constraints from the system.
 * The order of the non-redundant constraints can be shuffled after
 * the call.
 *
 * @return A reference to this object after removing all the
 *         redundant constraints
 */
LinearSystem &LinearSystem::simplify()
{
  unsigned int i = 0;
  unsigned int last_index = size() - 1;

  while (i < last_index) { // for every unchecked constraint
    // if it is redundant
    if (constraint_is_redundant(i)) {
      // swap it with the last constraint
      std::swap(_A[i], _A[last_index]);
      std::swap(_b[i], _b[last_index]);

      // remove the redundant constraint
      _A.resize(last_index);
      _b.resize(last_index);

      --last_index;
    } else { // otherwise, i.e. if it is not redundant

      // increase the index of the next constraint to be checked
      i++;
    }
  }

  if (constraint_is_redundant(last_index)) {
    _A.resize(last_index);
    _b.resize(last_index);
  }

  return *this;
}

/**
 * Remove redundant constraints
 */
LinearSystem LinearSystem::get_simplified() const
{
  LinearSystem simpler(*this);

  return simpler.simplify();
}
