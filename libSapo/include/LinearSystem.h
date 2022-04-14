/**
 * @file Polytope.h
 * Represent and manipulate a linear system
 * It can be used to represent polytopes (reached states, parameters, etc.)
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#ifndef LINEARSYSTEM_H_
#define LINEARSYSTEM_H_

#include <iostream>
#include <utility>
#include <vector>

#include "LinearAlgebra.h"
#include "SymbolicAlgebra.h"

#include "JSONStreamer.h"

/**
 * @brief The result of an optimization process
 *
 * @tparam T is the type of the output value
 */
template<typename T>
class OptimizationResult
{
  T _optimum;  //!< The optimum value
  int _status; //!< The status of the process
public:
  /**
   * @brief Construct a new optimization result object
   *
   * @param optimum is the optimum value
   * @param status is the status of the optimization process
   */
  OptimizationResult(const T &optimum, const int &status):
      _optimum(optimum), _status(status)
  {
  }

  /**
   * @brief Copy constructor
   *
   * @param orig is the template object
   */
  OptimizationResult(const OptimizationResult<T> &orig):
      _optimum(orig._optimum), _status(orig._status)
  {
  }

  /**
   * @brief Get the optimum value
   *
   * @return The optimum value
   */
  const T &optimum() const
  {
    return _optimum;
  }

  /**
   * @brief Get the status of the optimization process
   *
   * @return The status of the optimization process
   */
  const int &status() const
  {
    return _status;
  }

  OptimizationResult<T> &operator=(const OptimizationResult<T> &orig)
  {
    _optimum = orig._optimum;
    _status = orig._status;

    return *this;
  }
};

/**
 * Optimize a linear system
 *
 * @param[in] A template matrix of the system to optimize
 * @param[in] b offset vector of the system to optimize
 * @param[in] obj_fun objective function
 * @param[in] maximize is a Boolean parameter to establish whether
 *        maximize (true) or minimize (false) `obj_fun` over the system
 * @return optimum
 */
OptimizationResult<double> optimize(const std::vector<Vector<double>> &A,
                                    const Vector<double> &b,
                                    const Vector<double> &obj_fun,
                                    const bool maximize);

class LinearSystem
{

protected:
  std::vector<Vector<double>> A; // matrix A
  Vector<double> b;              // vector b

  /**
   * Check if a constraint belongs to the linear system
   *
   * @param[in] Ai direction
   * @param[in] bi offset
   * @returns true is Ai x <= bi is in the linear system
   */
  bool is_in(const Vector<double> &Ai, const double &bi) const;

  /**
   * Check whether the solutions of a linear system satisfy a constraint.
   *
   * This method establishes whether all the solutions of a linear system
   * satisfy a constraint. Due to approximation errors, it may return
   * false even if this is the case. However, whenever it returns true,
   * the all the solutions of the linear system certainly satisfy the
   * inequality.
   *
   * @param[in] i is the index of the constraint to be checked
   * @return a Boolean value. When some of the solutions of the linear
   *     system do not satisfy the inequality, the returned value is
   *     false. When the method returns true, the constraint is certainly
   *     satisfied by any of the solutions of the system. There are cases
   *     in which the constraint is satisfied by all the solutions and
   *     this method returns false.
   */
  bool satisfies(const Vector<double> &Ai, const double bi) const;

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
  bool constraint_is_redundant(const unsigned int i) const;

public:
  /**
   * Constructor that instantiates an empty linear system
   */
  LinearSystem();

  /**
   * Copy constructor
   *
   * @param[in] orig the original linear system
   */
  LinearSystem(const LinearSystem &orig);

  /**
   * Move constructor
   *
   * @param[in] orig the original linear system
   */
  LinearSystem(LinearSystem &&orig);

  /**
   * Constructor
   *
   * @param[in] A template matrix
   * @param[in] b offset vector
   */
  LinearSystem(const std::vector<Vector<double>> &A, const Vector<double> &b);

  /**
   * Move constructor
   *
   * @param[in] A template matrix
   * @param[in] b offset vector
   */
  LinearSystem(std::vector<Vector<double>> &&A, Vector<double> &&b);

  /**
   * Constructor from a set of symbolic expressions
   *
   * @param[in] vars list of variables appearing in the constraints
   * @param[in] constraints symbolic constraints
   */
  LinearSystem(const std::vector<SymbolicAlgebra::Symbol<>> &vars,
               const std::vector<SymbolicAlgebra::Expression<>> &constraints);

  /**
   * Deep assignament method
   *
   * @param[in] orig is the original model of the linear system
   * @return a reference to the modified object
   */
  LinearSystem &operator=(const LinearSystem &orig)
  {
    A = orig.A;
    b = orig.b;

    return *this;
  }

  /**
   * Return the template matrix
   *
   * @return template matrix
   */
  const std::vector<Vector<double>> &getA() const
  {
    return this->A;
  }

  /**
   * Return the offset vector
   *
   * @return offset vector
   */
  const Vector<double> &getb() const
  {
    return this->b;
  }

  /**
   * Return the (i,j) element of the template matrix
   *
   * @param[in] i row index
   * @param[in] j column index
   * @return (i,j) element
   */
  const double &getA(const unsigned int i, const unsigned int j) const;

  /**
   * Return the i-th element of the offset vector
   *
   * @param[in] i column index
   * @return i-th element
   */
  const double &getb(const unsigned int i) const;

  // optimization functions

  /**
   * Optimize a linear system
   *
   * @param[in] obj_fun objective function
   * @param[in] maximize is a Boolean parameter to establish whether
   *        maximize (true) or minimize (false) `obj_fun` over the system
   * @return optimum
   */
  OptimizationResult<double> optimize(const Vector<double> &obj_fun,
                                      const bool maximize) const;

  /**
   * Minimize the linear system
   *
   * @param[in] obj_fun objective function
   * @return minimum
   */
  OptimizationResult<double> minimize(const Vector<double> &obj_fun) const;

  /**
   * Maximize the linear system
   *
   * @param[in] obj_fun objective function
   * @return maximum
   */
  OptimizationResult<double> maximize(const Vector<double> &obj_fun) const;

  /**
   * Minimize the linear system
   *
   * @param[in] symbols array of symbols
   * @param[in] obj_fun objective function
   * @return minimum
   */
  OptimizationResult<double>
  minimize(const std::vector<SymbolicAlgebra::Symbol<>> &symbols,
           const SymbolicAlgebra::Expression<> &obj_fun) const;

  /**
   * Maximize the linear system
   *
   * @param[in] symbols array of symbols
   * @param[in] obj_fun objective function
   * @return maximum
   */
  OptimizationResult<double>
  maximize(const std::vector<SymbolicAlgebra::Symbol<>> &symbols,
           const SymbolicAlgebra::Expression<> &obj_fun) const;

  // testing methods

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
  bool has_solutions(const bool strict_inequality = false) const;

  /**
   * Remove redundant constraints from a linear system.
   *
   * This method removes redundant constraints from the system. The order
   * of the non-redundant constraints can be shuffled after the call.
   *
   * @return A reference to this object after removing all the redundant
   * constraints.
   */
  LinearSystem &simplify();

  /**
   * Create a new object and remove its redundant constraints.
   *
   * @return a copy of this object after all the redundant constraints of
   *     the former have been removed.
   */
  LinearSystem get_simplified() const;

  /**
   * Return the number of variables
   */
  unsigned int dim() const
  {
    if (size() == 0) {
      return 0;
    }

    return this->A[0].size();
  }

  /**
   * Return the number of inequalities
   */
  unsigned int size() const
  {
    return this->b.size();
  }

  friend void swap(LinearSystem &ls_1, LinearSystem &ls_2);
};

inline void swap(LinearSystem &ls_1, LinearSystem &ls_2)
{
  std::swap(ls_1.A, ls_2.A);
  std::swap(ls_1.b, ls_2.b);
}

std::ostream &operator<<(std::ostream &out, const LinearSystem &ls);

template<typename T>
JSON::ostream &operator<<(JSON::ostream &out, const std::vector<T> &v)
{
  out << "[";
  for (auto el_it = std::begin(v); el_it != std::end(v); ++el_it) {
    if (el_it != std::begin(v)) {
      out << ",";
    }
    out << *el_it;
  }
  out << "]";

  return out;
}

template<typename T>
JSON::ostream &operator<<(JSON::ostream &out,
                          const DenseLinearAlgebra::Matrix<T> &A)
{
  out << "[";
  for (auto row_it = std::begin(A); row_it != std::end(A); ++row_it) {
    if (row_it != std::begin(A)) {
      out << ",";
    }

    out << *row_it;
  }
  out << "]";

  return out;
}

JSON::ostream &operator<<(JSON::ostream &out, const LinearSystem &ls);

#endif /* LINEARSYSTEM_H_ */
