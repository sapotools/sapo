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

  /**
   * @brief Assignment operator
   * 
   * @param orig is the orginal optimization result
   * @return a reference to the updated object
   */
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
 * @param[in] maximize is a Boolean flag to establish whether
 *        `obj_fun` must be maximized (`true`) or minimized 
 *        (`false`) over the system
 * @return optimum
 */
OptimizationResult<double> optimize(const std::vector<LinearAlgebra::Vector<double>> &A,
                                    const LinearAlgebra::Vector<double> &b,
                                    const LinearAlgebra::Vector<double> &obj_fun,
                                    const bool maximize);

/**
 * @brief Linear systems of the form \f$A \cdot x \leq b\f$
 * 
 * This class represent linear systems having the form 
 * \f$A \cdot x \leq b\f$ where \f$x\f$ is a variable vector
 * and \f$A\f$ and \f$b\f$ are a matrix and a vector, 
 * respectively, with values in `double`s.
 */
class LinearSystem
{

protected:
  std::vector<LinearAlgebra::Vector<double>> A; //!< the matrix \f$A\f$
  LinearAlgebra::Vector<double> b;              //!< the vector \f$b\f$

  /**
   * Check if a constraint belongs to the linear system
   *
   * @param[in] Ai is direction
   * @param[in] bi js an offset
   * @returns `true` if and only if \f$\textrm{Ai} \cdot x \leq \textrm{bi}\f$ 
   *          is one of the rows of the linear system
   */
  bool is_in(const LinearAlgebra::Vector<double> &Ai, const double &bi) const;

  /**
   * @brief Check whether the solutions of a linear system satisfy a constraint
   *
   * This method establishes whether all the solutions \f$s\f$ of the linear 
   * system satisfy a constraint \f$\textrm{Ai} \cdot s \leq \textrm{bi}\f$. 
   * Due to the approximation errors, it may return `false` even if this is 
   * the case. However, whenever it returns `true`, all the solutions of the 
   * linear system certainly satisfy the inequality.
   * 
   * @param[in] Ai is a value vector
   * @param[in] bi is a scalar value
   * @return When some of the solutions \f$s\f$ of the linear system do not 
   *     satisfy the inequality \f$\textrm{Ai} \cdot s \leq \textrm{bi}\f$, 
   *     the returned value is `false`. When the method returns `true`, 
   *     the inequality is certainly satisfied by any of the solutions 
   *     \f$s\f$ of the system. There are cases in which the inequality
   *     is satisfied by all the solutions and this method returns `false`
   */
  bool satisfies(const LinearAlgebra::Vector<double> &Ai, const double& bi) const;

  /**
   * Check whether one of the constraints in a linear system is redundant.
   *
   * This method establishes whether the i-th constraint of a linear
   * system is redundant. Due to approximation errors, it may return
   * `false` even if the constraint is redundant. However, whenever it
   * returns `true`, the constraint is certainly redundant.
   *
   * @param[in] i is the index of the constraint to be checked
   * @return a Boolean value. When the constraint is non-redundanct, the
   *     returned value is true. When it returns `true`, the constraint is
   *     certainly redundant. There are cases in which the constraint is
   *     redundant and this method returns `false`.
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
   * @param[in] ls is the original linear system
   */
  LinearSystem(const LinearSystem &ls);

  /**
   * Swap constructor
   *
   * @param[in] ls is the linear system
   */
  LinearSystem(LinearSystem &&ls);

  /**
   * Constructor
   *
   * @param[in] A is the matrix 
   * @param[in] b is the offset vector
   */
  LinearSystem(const std::vector<LinearAlgebra::Vector<double>> &A, const LinearAlgebra::Vector<double> &b);

  /**
   * Move constructor
   *
   * @param[in] A is the matrix 
   * @param[in] b is the offset vector
   */
  LinearSystem(std::vector<LinearAlgebra::Vector<double>> &&A, LinearAlgebra::Vector<double> &&b);

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
  LinearSystem(const std::vector<SymbolicAlgebra::Symbol<>> &x,
               const std::vector<SymbolicAlgebra::Expression<>> &expressions);

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
  const std::vector<LinearAlgebra::Vector<double>> &getA() const
  {
    return this->A;
  }

  /**
   * Return the offset vector
   *
   * @return offset vector
   */
  const LinearAlgebra::Vector<double> &getb() const
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
   * @param[in] maximize is a Boolean flag to establish whether
 *        `obj_fun` must be maximized (`true`) or minimized 
 *        (`false`) over the system
   * @return optimum
   */
  OptimizationResult<double> optimize(const LinearAlgebra::Vector<double> &obj_fun,
                                      const bool maximize) const;

  /**
   * Minimize the linear system
   *
   * @param[in] obj_fun objective function
   * @return minimum
   */
  OptimizationResult<double> minimize(const LinearAlgebra::Vector<double> &obj_fun) const;

  /**
   * Maximize the linear system
   *
   * @param[in] obj_fun objective function
   * @return maximum
   */
  OptimizationResult<double> maximize(const LinearAlgebra::Vector<double> &obj_fun) const;

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
   * Due to approximation errors, it may return `true` for some systems
   * having no solution too. However, when it returns `false`, the 
   * linear system certainly has no solution.
   *
   * @param[in] strict_inequality specifies whether the linear system is
   *         a strict inequality (i.e., \f$Ax < b\f$)
   * @return a Boolean value. If the returned value is false, then the
   *       linear system has no solution
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

  /**
   * @brief Swap two linear systems
   * 
   * @param ls_1 is one of the linear system to be swapped
   * @param ls_2 is the second linear system to be swapped
   */
  friend inline void swap(LinearSystem &ls_1, LinearSystem &ls_2)
  {
    std::swap(ls_1.A, ls_2.A);
    std::swap(ls_1.b, ls_2.b);
  }
};

/**
 * @brief Print a linear system in a stream
 * 
 * @param out is the output stream
 * @param ls is the linear system to be print
 * @return a reference to the output stream
 */
std::ostream &operator<<(std::ostream &out, const LinearSystem &ls);

/**
 * @brief Print a linear system in a JSON stream
 * 
 * @param out is the output JSON stream
 * @param ls is the linear system to be print
 * @return a reference to the output JSON stream
 */
JSON::ostream &operator<<(JSON::ostream &out, const LinearSystem &ls);

#endif /* LINEARSYSTEM_H_ */
