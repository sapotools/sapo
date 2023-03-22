/**
 * @file LinearSystem.h
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @author Alberto Casagrande <acasagrande@units.it>
 * @brief Represent and manipulate linear systems
 * @version 0.2
 * @date 2015-10-14
 *
 * @copyright Copyright (c) 2015-2022
 */

#ifndef LINEARSYSTEM_H_
#define LINEARSYSTEM_H_

#include <iostream>
#include <utility>
#include <vector>

#include "LinearAlgebra.h"
#include "SymbolicAlgebra.h"
#include "Simplex.h"

/**
 * @brief Linear systems of the form \f$A \cdot x \leq b\f$
 *
 * This class represent linear systems having the form
 * \f$A \cdot x \leq b\f$ where \f$x\f$ is a variable vector
 * and \f$A\f$ and \f$b\f$ are a matrix and a vector,
 * respectively, with values in `double`s.
 *
 * @tparam T is the type of the linear system coefficients
 * @tparam APPROX_TYPE is the type of approximation used in place of T
 *      during the computation
 */
template<typename T = double, typename APPROX_TYPE = T>
class LinearSystem
{

protected:
  std::vector<LinearAlgebra::Vector<T>> _A; //!< the matrix \f$A\f$
  LinearAlgebra::Vector<T> _b;              //!< the vector \f$b\f$

  /**
   * Check if a constraint belongs to the linear system
   *
   * @param[in] Ai is direction
   * @param[in] bi js an offset
   * @returns `true` if and only if
   *          \f$\textrm{Ai} \cdot x \leq \textrm{bi}\f$
   *          is one of the rows of the linear system
   */
  bool contains(const LinearAlgebra::Vector<T> &Ai, const T &bi) const;

  /**
   * @brief Check whether the solutions of a linear system satisfy a constraint
   *
   * This method establishes whether all the solutions \f$s\f$ of the linear
   * system satisfy a constraint \f$\textrm{Ai} \cdot s \leq \textrm{bi}\f$.
   *
   * @param[in] Ai is a value vector
   * @param[in] bi is a scalar value
   * @return `true` if the method can establish that
   *         \f$\textrm{Ai} \cdot s \leq \textrm{bi}\f$ holds
   *         for all the solutions \f$s\f$ of the linear system.
   *         `false` if the method can establish that
   *         \f$\textrm{Ai} \cdot s \leq \textrm{bi}\f$ does not
   *         hold for some of the solutions \f$s\f$ of the linear
   *         system. `uncertain` in the remaining cases
   */
  TriBool satisfies(const LinearAlgebra::Vector<T> &Ai, const T &bi) const;

  /**
   * Check whether one of the constraints in a linear system is redundant.
   *
   * This method establishes whether the i-th constraint of a linear
   * system is redundant.
   *
   * @param[in] i is the index of the constraint to be checked
   * @return `true` if the method can establish that the `i`-th
   *         constraint is redundant. `false` if it can establish
   *         that the `i`-th constraint is not redundant.
   *         `uncertain` is the remaining cases
   */
  TriBool constraint_is_redundant(const size_t i) const;

public:
  /**
   * Constructor that instantiates an empty linear system
   */
  LinearSystem();

  /**
   * Copy constructor
   *
   * @tparam APPROX2 is the approximation type of the
   *         original linear system
   * @param[in] ls is the original linear system
   */
  template<typename APPROX2>
  LinearSystem(const LinearSystem<T, APPROX2> &ls);

  /**
   * Swap constructor
   *
   * @tparam APPROX2 is the approximation type of the
   *         original linear system
   * @param[in] ls is the linear system
   */
  template<typename APPROX2>
  LinearSystem(LinearSystem<T, APPROX2> &&ls);

  /**
   * Constructor
   *
   * @param[in] A is the matrix
   * @param[in] b is the offset vector
   */
  LinearSystem(const std::vector<LinearAlgebra::Vector<T>> &A,
               const LinearAlgebra::Vector<T> &b);

  /**
   * Move constructor
   *
   * @param[in] A is the matrix
   * @param[in] b is the offset vector
   */
  LinearSystem(std::vector<LinearAlgebra::Vector<T>> &&A,
               LinearAlgebra::Vector<T> &&b);

  /**
   * Constructor from a set of symbolic expressions
   *
   * This constructor builds a linear system \f$A \cdot x \leq b\f$
   * by taking in input the variable vector \f$x\f$ and an array
   * `expressions` of symbolic expressions. The i-th expression
   * defines the i-th row in the matrix \f$A\f$ and the i-th
   * values of the vector \f$b\f$. For instance, if `x` is the
   * vector \f$[x_0, x_1, x_2, x_3]\f$ and the 3rd expression
   * in `expressions` is \f$3*x_2 +5 -3/4 * x_0\f$ when the
   * third row in \f$A\f$ will be \f$[-3/4, 0, 3, 0]\f$ and
   * the third value in `b` \f$5\f$.
   *
   * @param[in] x is the variable vector of the linear system
   * @param[in] expressions symbolic constraints
   * @todo Test expressions linearity
   */
  LinearSystem(const std::vector<SymbolicAlgebra::Symbol<T>> &x,
               const std::vector<SymbolicAlgebra::Expression<T>> &expressions);

  /**
   * Deep assignment method
   *
   * @tparam APPROX2 is the approximation type of the
   *         original linear system
   * @param[in] orig is the original model of the linear system
   * @return a reference to the modified object
   */
  template<typename APPROX2>
  LinearSystem<T, APPROX_TYPE> &
  operator=(const LinearSystem<T, APPROX2> &orig);

  /**
   * @brief Add a new constraint to the system
   *
   * This method add a new linear constraint of the type \f$v\cdot x <= b\f$ to
   * the system
   *
   * @param v is the vector of the variable coefficents
   * @param b is the constant term
   * @return a reference to the updated linear system
   */
  LinearSystem<T, APPROX_TYPE> &
  add_constraint(const LinearAlgebra::Vector<T> &v, const T &b);

  /**
   * @brief Add new constraints to the system
   *
   * This method add all the constraints of a linear system to the
   * current linear system.
   *
   * @param linear_system is the linear system whose constraints must be
   *      added to the current one
   * @return a reference to the updated linear system
   */
  LinearSystem<T, APPROX_TYPE> &
  add_constraints(const LinearSystem<T, APPROX_TYPE> &linear_system);

  /**
   * Return the template matrix
   *
   * @return template matrix
   */
  const std::vector<LinearAlgebra::Vector<T>> &A() const;

  /**
   * Return the offset vector
   *
   * @return offset vector
   */
  const LinearAlgebra::Vector<T> &b() const;

  /**
   * Return the i-th row of the matrix
   *
   * @param[in] i row index
   * @return the i-th row
   */
  const LinearAlgebra::Vector<T> &A(const size_t i) const;

  /**
   * Return the i-th element of the offset vector
   *
   * @param[in] i column index
   * @return i-th element
   */
  const T &b(const size_t i) const;

  // optimization functions

  /**
   * Optimize a linear system
   *
   * @param[in] obj_fun objective function
   * @param[in] maximize is a Boolean flag to establish whether
   *        `obj_fun` must be maximized (`true`) or minimized
   *        (`false`) over the system
   * @return optimum
   */
  OptimizationResult<APPROX_TYPE>
  optimize(const LinearAlgebra::Vector<T> &obj_fun, const bool maximize) const;

  /**
   * Minimize the linear system
   *
   * @param[in] obj_fun objective function
   * @return minimum
   */
  OptimizationResult<APPROX_TYPE>
  minimize(const LinearAlgebra::Vector<T> &obj_fun) const;

  /**
   * Maximize the linear system
   *
   * @param[in] obj_fun objective function
   * @return maximum
   */
  OptimizationResult<APPROX_TYPE>
  maximize(const LinearAlgebra::Vector<T> &obj_fun) const;

  /**
   * Minimize the linear system
   *
   * @param[in] symbols array of symbols
   * @param[in] obj_fun objective function
   * @return minimum
   */
  OptimizationResult<APPROX_TYPE>
  minimize(const std::vector<SymbolicAlgebra::Symbol<T>> &symbols,
           SymbolicAlgebra::Expression<T> obj_fun) const;

  /**
   * Maximize the linear system
   *
   * @param[in] symbols array of symbols
   * @param[in] obj_fun objective function
   * @return maximum
   */
  OptimizationResult<APPROX_TYPE>
  maximize(const std::vector<SymbolicAlgebra::Symbol<T>> &symbols,
           SymbolicAlgebra::Expression<T> obj_fun) const;

  /**
   * Establish whether a linear system has solutions
   *
   * @param[in] strict_inequality specifies whether the linear system is
   *         a strict inequality (i.e., \f$Ax < b\f$)
   * @return `true` when it can establish that the linear system has
   *         solutions. `false` when it can establish that the linear
   *         system does not admit any solution. `uncertain` when in
   *         doubt
   */
  TriBool has_solutions(const bool strict_inequality = false) const;

  /**
   * @brief Check whether all the solutions are also solutions for another
   * system
   *
   * @param ls is the considered linear system
   * @return `true` when the method can establish that the solutions of
   *         the current object are also solutions for `ls`. `false`
   *         when it can establish that not all the solutions of
   *         the current object are also solutions for `ls`.
   *         `uncertain` in the remaining cases
   */
  TriBool satisfies(const LinearSystem<T, APPROX_TYPE> &ls) const;

  /**
   * Remove redundant constraints from a linear system
   *
   * This method removes redundant constraints from the system. The order
   * of the non-redundant constraints can be shuffled after the call.
   *
   * @return A reference to this object after removing all the redundant
   *         constraints
   */
  LinearSystem<T, APPROX_TYPE> &simplify();

  /**
   * Create a new object and remove its redundant constraints
   *
   * @return a copy of this object after all the redundant constraints of
   *         the former have been removed
   */
  LinearSystem<T, APPROX_TYPE> get_simplified() const;

  /**
   * Return the number of variables
   */
  size_t dim() const;

  /**
   * Return the number of inequalities
   */
  size_t size() const;

  template<typename T2, typename APPROX2>
  friend void swap(LinearSystem<T2, APPROX2> &ls_1,
                   LinearSystem<T2, APPROX2> &ls_2);
};

/**
 * @brief Test whether two linear systems do not have common solutions
 *
 * @tparam T is the type of the linear system coefficients
 * @tparam APPROX_TYPE is the type of approximation used in place of T
 *      during the computation
 * @param ls1 is a linear system
 * @param ls2 is a linear system
 * @return `true` if and only if the solutions of `ls1` and `ls2` are
 *         disjoint
 * @return `true` when the method can establish that the solution sets
 *         of `ls1` and `ls2` are disjoint. `false` when it can
 *         establish that the solution sets are not dijoint.
 *         `uncertain` in the remaining cases
 */
template<typename T, typename APPROX_TYPE>
TriBool have_disjoint_solutions(const LinearSystem<T, APPROX_TYPE> &ls1,
                                const LinearSystem<T, APPROX_TYPE> &ls2);

/**
 * @brief Swap two linear systems
 *
 * @tparam T is the type of the linear system coefficients
 * @tparam APPROX_TYPE is the approximation type of the parameters
 * @param ls_1 is one of the linear system to be swapped
 * @param ls_2 is the second linear system to be swapped
 */
template<typename T, typename APPROX_TYPE>
void swap(LinearSystem<T, APPROX_TYPE> &ls_1,
          LinearSystem<T, APPROX_TYPE> &ls_2);

/**
 * @brief Test whether two linear systems are the same one
 *
 * @tparam T is the type of the linear system coefficients
 * @tparam APPROX_TYPE is the type of approximation used in place of T
 *      during the computation
 * @param P1 is the first linear system
 * @param P2 is the second linear system
 * @return `true` when the method can establish that the solution
 *         sets of two linear systems are the same set. `false`
 *         when it can establish that the solution sets are not
 *         the same. `uncertain` in the remaining cases
 */
template<typename T, typename APPROX_TYPE>
TriBool operator==(const LinearSystem<T, APPROX_TYPE> &P1,
                   const LinearSystem<T, APPROX_TYPE> &P2);

/**
 * @brief Test whether two linear systems differ
 *
 * @tparam T is the type of the linear system coefficients
 * @tparam APPROX_TYPE is the type of approximation used in place of T
 *      during the computation
 * @param P1 is the first linear system
 * @param P2 is the second linear system
 * @return `true` when the method can establish that the solution
 *         sets of two linear systems are not the same set.
 *         `false` when it can establish that the solution sets
 *         are the same. `uncertain` in the remaining cases
 */
template<typename T, typename APPROX_TYPE>
TriBool operator!=(const LinearSystem<T, APPROX_TYPE> &P1,
                   const LinearSystem<T, APPROX_TYPE> &P2);

/**
 * @brief Print a linear system in a stream
 *
 * @tparam T is the type of the linear system coefficients
 * @tparam APPROX_TYPE is the type of approximation used in place of T
 *      during the computation
 * @param out is the output stream
 * @param ls is the linear system to be print
 * @return a reference to the output stream
 */
template<typename T, typename APPROX_TYPE>
std::ostream &operator<<(std::ostream &out,
                         const LinearSystem<T, APPROX_TYPE> &ls);

template<typename T, typename APPROX_TYPE>
LinearSystem<T, APPROX_TYPE>::LinearSystem(
    const std::vector<LinearAlgebra::Vector<T>> &A,
    const LinearAlgebra::Vector<T> &b)
{
  if (A.size() != b.size()) {
    SAPO_ERROR("the number of rows in A must equals the "
               "number of elements in b",
               std::domain_error);
  }

  for (size_t i = 0; i < A.size(); i++) {
    if (A[i].size() != A[0].size()) {
      SAPO_ERROR("the " << (i + 1) << "th row"
                        << " and the first one"
                        << " differ in size",
                 std::domain_error);
    }
  }

#ifdef SMART_INSERT
  for (unsigned int i = 0; i < A.size(); i++) {
    if ((LinearAlgebra::norm_infinity(A[i]) > 0)
        && !is_true(this->satisfies(A[i], b[i]))) {
      this->_A.push_back(A[i]);
      this->_b.push_back(b[i]);
    }
  }
#else  // SMART_INSERT
  this->_A = A;
  this->_b = b;
#endif // SMART_INSERT
}

template<typename T, typename APPROX_TYPE>
LinearSystem<T, APPROX_TYPE>::LinearSystem(
    std::vector<LinearAlgebra::Vector<T>> &&A, LinearAlgebra::Vector<T> &&b):
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
      SAPO_ERROR("the " << (i + 1) << "th row"
                        << " and the first one"
                        << " differ in size",
                 std::domain_error);
    }
  }
}

template<typename T, typename APPROX_TYPE>
LinearSystem<T, APPROX_TYPE>::LinearSystem(): _A(), _b()
{
}

template<typename T, typename APPROX_TYPE>
template<typename APPROX2>
LinearSystem<T, APPROX_TYPE>::LinearSystem(const LinearSystem<T, APPROX2> &ls):
    _A(ls._A), _b(ls._b)
{
}

template<typename T, typename APPROX_TYPE>
template<typename APPROX2>
LinearSystem<T, APPROX_TYPE>::LinearSystem(LinearSystem<T, APPROX2> &&ls)
{
  swap(this->_A, ls._A);
  swap(this->_b, ls._b);
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
                     const LinearAlgebra::Vector<double> &A2,
                     const double &b2);

template<typename T, typename APPROX_TYPE>
bool LinearSystem<T, APPROX_TYPE>::contains(const LinearAlgebra::Vector<T> &Ai,
                                            const T &bi) const
{
  for (size_t i = 0; i < this->_A.size(); i++) {
    if ((this->_A[i] == Ai) && (this->_b[i] == bi)) {
      return true;
    }
  }
  return false;
}

template<typename T, typename APPROX_TYPE>
LinearSystem<T, APPROX_TYPE>::LinearSystem(
    const std::vector<SymbolicAlgebra::Symbol<T>> &x,
    const std::vector<SymbolicAlgebra::Expression<T>> &expressions)
{
  using namespace SymbolicAlgebra;
  std::vector<Expression<T>> lexpressions = expressions;
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

    T bi = const_term.evaluate();

    if (!this->contains(Ai, -bi)) {
      this->_A.push_back(Ai);
      this->_b.push_back(-bi);
    }
  }
}

template<typename T, typename APPROX_TYPE>
LinearSystem<T, APPROX_TYPE> &
LinearSystem<T, APPROX_TYPE>::add_constraint(const LinearAlgebra::Vector<T> &v,
                                             const T &b)
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

template<typename T, typename APPROX_TYPE>
LinearSystem<T, APPROX_TYPE> &LinearSystem<T, APPROX_TYPE>::add_constraints(
    const LinearSystem<T, APPROX_TYPE> &linear_system)
{
  if (linear_system.size() == 0) {
    return *this;
  }

  if (this->dim() != linear_system.dim()) {
    SAPO_ERROR("the two linear systems differ "
               "in the number of variables",
               std::domain_error);
  }

  for (const auto &direction: linear_system._A) {
    _A.push_back(direction);
  }

  for (const auto &b: linear_system._b) {
    _b.push_back(b);
  }

  return *this;
}

template<typename T, typename APPROX_TYPE>
inline const std::vector<LinearAlgebra::Vector<T>> &
LinearSystem<T, APPROX_TYPE>::A() const
{
  return this->_A;
}

template<typename T, typename APPROX_TYPE>
inline const LinearAlgebra::Vector<T> &LinearSystem<T, APPROX_TYPE>::b() const
{
  return this->_b;
}

template<typename T, typename APPROX_TYPE>
const LinearAlgebra::Vector<T> &
LinearSystem<T, APPROX_TYPE>::A(const size_t i) const
{
  if (i >= this->_A.size()) {
    SAPO_ERROR("the parameter must be a valid row index "
               "for the linear system matrix A",
               std::domain_error);
  }
  return this->_A[i];
}

template<typename T, typename APPROX_TYPE>
const T &LinearSystem<T, APPROX_TYPE>::b(const size_t i) const
{
  if (i >= this->_b.size()) {
    SAPO_ERROR("the parameter must be a valid index for the linear "
               "system coefficient vector b",
               std::domain_error);
  }
  return this->_b[i];
}

template<typename T, typename APPROX_TYPE>
OptimizationResult<APPROX_TYPE>
LinearSystem<T, APPROX_TYPE>::optimize(const LinearAlgebra::Vector<T> &obj_fun,
                                       const bool maximize) const
{
  SimplexMethodOptimizer<T, APPROX_TYPE> optimizer;

  return optimizer(this->_A, this->_b, obj_fun,
                   (maximize ? MAXIMIZE : MINIMIZE));
}

template<typename T, typename APPROX_TYPE>
inline OptimizationResult<APPROX_TYPE> LinearSystem<T, APPROX_TYPE>::minimize(
    const LinearAlgebra::Vector<T> &obj_fun) const
{
  return this->optimize(obj_fun, false);
}

template<typename T, typename APPROX_TYPE>
inline OptimizationResult<APPROX_TYPE> LinearSystem<T, APPROX_TYPE>::maximize(
    const LinearAlgebra::Vector<T> &obj_fun) const
{
  return this->optimize(obj_fun, true);
}

/**
 * @brief Extract the coeffiecients from a linear expressions
 *
 * @tparam T is the constant type
 * @param[in] symbols is the vector of the symbols to be processed
 * @param[in] function is the function whose coefficients are aimed
 * @param[out] coeffs is the vector of extracted coefficients from  `function`
 * @param[out] constant_term is the constant term in `function`
 */
template<typename T>
void extract_linear_coefficients(
    const std::vector<SymbolicAlgebra::Symbol<T>> &symbols,
    SymbolicAlgebra::Expression<T> function, LinearAlgebra::Vector<T> &coeffs,
    T &constant_term)
{
  using namespace SymbolicAlgebra;

  function.expand();

  coeffs.resize(0);
  Expression<T> const_term(function);

  // Extract the coefficient of the i-th variable (grade 1)
  for (const auto &symbol: symbols) {
    if (function.degree(symbol) > 1) {
      SAPO_ERROR("the function is assumed to be linear. The degree of \""
                     << symbol << "\" is " << function.degree(symbol) << ".",
                 std::domain_error);
    }

    coeffs.push_back((function.get_coeff(symbol, 1)).evaluate());
    const_term = const_term.get_coeff(symbol, 0);
  }

  constant_term = const_term.evaluate();
}

template<typename T, typename APPROX_TYPE>
OptimizationResult<APPROX_TYPE> LinearSystem<T, APPROX_TYPE>::minimize(
    const std::vector<SymbolicAlgebra::Symbol<T>> &symbols,
    SymbolicAlgebra::Expression<T> obj_fun) const
{
  LinearAlgebra::Vector<T> obj_fun_coeffs;
  T const_term;

  extract_linear_coefficients(symbols, obj_fun, obj_fun_coeffs, const_term);

  auto res = minimize(obj_fun_coeffs);

  if (res.status() != res.OPTIMUM_AVAILABLE) {
    return res;
  }

  return OptimizationResult<APPROX_TYPE>(res.optimum(),
                                         res.objective_value() + const_term);
}

template<typename T, typename APPROX_TYPE>
OptimizationResult<APPROX_TYPE> LinearSystem<T, APPROX_TYPE>::maximize(
    const std::vector<SymbolicAlgebra::Symbol<T>> &symbols,
    SymbolicAlgebra::Expression<T> obj_fun) const
{
  LinearAlgebra::Vector<T> obj_fun_coeffs;
  T const_term;

  extract_linear_coefficients(symbols, obj_fun, obj_fun_coeffs, const_term);

  auto res = maximize(obj_fun_coeffs);

  if (res.status() != res.OPTIMUM_AVAILABLE) {
    return res;
  }

  return OptimizationResult<APPROX_TYPE>(res.optimum(),
                                         res.objective_value() + const_term);
}

template<typename T, typename APPROX_TYPE>
TriBool
LinearSystem<T, APPROX_TYPE>::satisfies(const LinearAlgebra::Vector<T> &Ai,
                                        const T &bi) const
{
  if (size() == 0) {
    return false;
  }

  OptimizationResult<APPROX_TYPE> res = this->maximize(Ai);

  return ((res.status() == res.OPTIMUM_AVAILABLE)
          && res.objective_value() <= bi);
}

template<typename T, typename APPROX_TYPE>
TriBool LinearSystem<T, APPROX_TYPE>::satisfies(
    const LinearSystem<T, APPROX_TYPE> &ls) const
{
  if (is_false(this->has_solutions())) {
    return true;
  }

  TriBool sat{true};

  for (unsigned int i = 0; i < ls.size(); i++) {
    sat = sat && this->satisfies(ls._A[i], ls._b[i]);

    if (is_false(sat)) {
      return false;
    }
  }

  return sat;
}

template<typename T, typename APPROX_TYPE>
TriBool
LinearSystem<T, APPROX_TYPE>::has_solutions(const bool strict_inequality) const
{
  if (this->size() == 0) {
    return true;
  }

  SimplexMethodOptimizer<T, APPROX_TYPE> optimizer;

  auto feasible = optimizer.is_feasible(_A, _b);

  if (is_false(feasible) || !strict_inequality) {
    return feasible;
  }

  TriBool result{true};

  for (const auto &row: _A) {

    auto max_row = optimizer(_A, _b, row, MAXIMIZE);
    auto min_row = optimizer(_A, _b, row, MINIMIZE);

    if (min_row.feasible_set_is_empty() || max_row.feasible_set_is_empty()) {
      return false;
    }

    result = result && (min_row.objective_value() < max_row.objective_value());

    if (is_false(result)) {
      return result;
    }
  }

  return result;
}

template<typename T, typename APPROX_TYPE>
template<typename APPROX2>
LinearSystem<T, APPROX_TYPE> &
LinearSystem<T, APPROX_TYPE>::operator=(const LinearSystem<T, APPROX2> &orig)
{
  _A = orig._A;
  _b = orig._b;

  return *this;
}

template<typename T, typename APPROX_TYPE>
TriBool
LinearSystem<T, APPROX_TYPE>::constraint_is_redundant(const size_t i) const
{
  LinearSystem<T, APPROX_TYPE> tmp(*this);
  LinearAlgebra::Vector<T> Ai(dim(), 0);
  T bi(0);

  // replace the i-th constraint with the empty constraint
  std::swap(Ai, tmp._A[i]);
  std::swap(bi, tmp._b[i]);

  // check whether the i-th constraint is redundant
  TriBool result(tmp.satisfies(Ai, bi));

  // reinsert the i-th into the system
  std::swap(Ai, tmp._A[i]);
  std::swap(bi, tmp._b[i]);

  return result;
}

template<typename T, typename APPROX_TYPE>
LinearSystem<T, APPROX_TYPE> &LinearSystem<T, APPROX_TYPE>::simplify()
{
  unsigned int i = 0;
  unsigned int last_index = size() - 1;

  while (i < last_index) { // for every unchecked constraint
    // if it is redundant
    if (is_true(constraint_is_redundant(i))) {
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

  if (is_true(constraint_is_redundant(last_index))) {
    _A.resize(last_index);
    _b.resize(last_index);
  }

  return *this;
}

template<typename T, typename APPROX_TYPE>
inline LinearSystem<T, APPROX_TYPE>
LinearSystem<T, APPROX_TYPE>::get_simplified() const
{
  LinearSystem<T, APPROX_TYPE> simpler(*this);

  return simpler.simplify();
}

template<typename T, typename APPROX_TYPE>
inline size_t LinearSystem<T, APPROX_TYPE>::dim() const
{
  if (size() == 0) {
    return 0;
  }

  return this->_A[0].size();
}

template<typename T, typename APPROX_TYPE>
inline size_t LinearSystem<T, APPROX_TYPE>::size() const
{
  return this->_b.size();
}

template<typename T, typename APPROX_TYPE>
inline void swap(LinearSystem<T, APPROX_TYPE> &ls_1,
                 LinearSystem<T, APPROX_TYPE> &ls_2)
{
  std::swap(ls_1._A, ls_2._A);
  std::swap(ls_1._b, ls_2._b);
}

template<typename T, typename APPROX_TYPE>
TriBool have_disjoint_solutions(const LinearSystem<T, APPROX_TYPE> &ls1,
                                const LinearSystem<T, APPROX_TYPE> &ls2)
{
  if (ls1.dim() != ls2.dim()) {
    SAPO_ERROR("the two linear systems must have the "
               "same number of variables",
               std::domain_error);
  }

  if (ls1.dim() == 0) {
    return false;
  }

  auto A{ls1.A()};
  A.reserve(ls1.size() + ls2.size());
  std::copy(std::begin(ls2.A()), std::end(ls2.A()), std::back_inserter(A));

  auto b{ls1.b()};
  b.reserve(ls1.size() + ls2.size());
  std::copy(std::begin(ls2.b()), std::end(ls2.b()), std::back_inserter(b));

  SimplexMethodOptimizer<T, APPROX_TYPE> optimizer;

  return !optimizer.is_feasible(A, b);
}

template<typename T, typename APPROX_TYPE>
inline TriBool operator==(const LinearSystem<T, APPROX_TYPE> &P1,
                          const LinearSystem<T, APPROX_TYPE> &P2)
{
  return P1.satisfies(P2) && P2.satisfies(P1);
}

template<typename T, typename APPROX_TYPE>
inline TriBool operator!=(const LinearSystem<T, APPROX_TYPE> &P1,
                          const LinearSystem<T, APPROX_TYPE> &P2)
{
  return !(P1 == P2);
}

template<typename T, typename APPROX_TYPE>
std::ostream &operator<<(std::ostream &out,
                         const LinearSystem<T, APPROX_TYPE> &ls)
{
  for (size_t row_idx = 0; row_idx < ls.size(); row_idx++) {
    if (row_idx != 0) {
      out << std::endl;
    }

    for (size_t col_idx = 0; col_idx < ls.dim(); col_idx++) {
      out << ls.A(row_idx)[col_idx] << " ";
    }
    out << "<= " << ls.b(row_idx);
  }

  return out;
}

#endif /* LINEARSYSTEM_H_ */
