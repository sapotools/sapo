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
 */
class LinearSystem
{

protected:
  std::vector<LinearAlgebra::Vector<double>> _A; //!< the matrix \f$A\f$
  LinearAlgebra::Vector<double> _b;              //!< the vector \f$b\f$

  /**
   * Check if a constraint belongs to the linear system
   *
   * @param[in] Ai is direction
   * @param[in] bi js an offset
   * @returns `true` if and only if \f$\textrm{Ai} \cdot x \leq \textrm{bi}\f$
   *          is one of the rows of the linear system
   */
  bool contains(const LinearAlgebra::Vector<double> &Ai,
                const double &bi) const;

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
  bool satisfies(const LinearAlgebra::Vector<double> &Ai,
                 const double &bi) const;

  /**
   * Check whether one of the constraints in a linear system is redundant.
   *
   * This method establishes whether the i-th constraint of a linear
   * system is redundant.
   *
   * @param[in] i is the index of the constraint to be checked
   * @return `true` if and only if the `i`-th constraint is redundant
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
  LinearSystem(const std::vector<LinearAlgebra::Vector<double>> &A,
               const LinearAlgebra::Vector<double> &b);

  /**
   * Move constructor
   *
   * @param[in] A is the matrix
   * @param[in] b is the offset vector
   */
  LinearSystem(std::vector<LinearAlgebra::Vector<double>> &&A,
               LinearAlgebra::Vector<double> &&b);

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
   * Deep assignment method
   *
   * @param[in] orig is the original model of the linear system
   * @return a reference to the modified object
   */
  LinearSystem &operator=(const LinearSystem &orig);

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
  LinearSystem &add_constraint(const LinearAlgebra::Vector<double> &v,
                               const double &b);

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
  LinearSystem &add_constraints(const LinearSystem &linear_system);

  /**
   * Return the template matrix
   *
   * @return template matrix
   */
  const std::vector<LinearAlgebra::Vector<double>> &A() const
  {
    return this->_A;
  }

  /**
   * Return the offset vector
   *
   * @return offset vector
   */
  const LinearAlgebra::Vector<double> &b() const
  {
    return this->_b;
  }

  /**
   * Return the i-th row of the matrix
   *
   * @param[in] i row index
   * @return the i-th row
   */
  const LinearAlgebra::Vector<double> &A(const unsigned int i) const;

  /**
   * Return the i-th element of the offset vector
   *
   * @param[in] i column index
   * @return i-th element
   */
  const double &b(const unsigned int i) const;

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
  OptimizationResult<double>
  optimize(const LinearAlgebra::Vector<double> &obj_fun,
           const bool maximize) const;

  /**
   * Minimize the linear system
   *
   * @param[in] obj_fun objective function
   * @return minimum
   */
  OptimizationResult<double>
  minimize(const LinearAlgebra::Vector<double> &obj_fun) const;

  /**
   * Maximize the linear system
   *
   * @param[in] obj_fun objective function
   * @return maximum
   */
  OptimizationResult<double>
  maximize(const LinearAlgebra::Vector<double> &obj_fun) const;

  /**
   * Minimize the linear system
   *
   * @param[in] symbols array of symbols
   * @param[in] obj_fun objective function
   * @return minimum
   */
  OptimizationResult<double>
  minimize(const std::vector<SymbolicAlgebra::Symbol<>> &symbols,
           SymbolicAlgebra::Expression<> obj_fun) const;

  /**
   * Maximize the linear system
   *
   * @param[in] symbols array of symbols
   * @param[in] obj_fun objective function
   * @return maximum
   */
  OptimizationResult<double>
  maximize(const std::vector<SymbolicAlgebra::Symbol<>> &symbols,
           SymbolicAlgebra::Expression<> obj_fun) const;

  /**
   * Establish whether a linear system has solutions
   *
   * @param[in] strict_inequality specifies whether the linear system is
   *         a strict inequality (i.e., \f$Ax < b\f$)
   * @return `true` if and only if the linear system has solutions
   */
  bool has_solutions(const bool strict_inequality = false) const;

  /**
   * @brief Check whether all the solutions are also solutions for another
   * system
   *
   * @param ls is the considered linear system
   * @return `true` if and only if the solutions of the current object are
   *         also solutions for `ls`
   */
  bool satisfies(const LinearSystem &ls) const;

  /**
   * Remove redundant constraints from a linear system
   *
   * This method removes redundant constraints from the system. The order
   * of the non-redundant constraints can be shuffled after the call.
   *
   * @return A reference to this object after removing all the redundant
   *         constraints
   */
  LinearSystem &simplify();

  /**
   * Create a new object and remove its redundant constraints
   *
   * @return a copy of this object after all the redundant constraints of
   *         the former have been removed
   */
  LinearSystem get_simplified() const;

  /**
   * Return the number of variables
   */
  inline unsigned int dim() const
  {
    if (size() == 0) {
      return 0;
    }

    return this->_A[0].size();
  }

  /**
   * Return the number of inequalities
   */
  inline unsigned int size() const
  {
    return this->_b.size();
  }

  /**
   * @brief Swap two linear systems
   *
   * @param ls_1 is one of the linear system to be swapped
   * @param ls_2 is the second linear system to be swapped
   */
  friend inline void swap(LinearSystem &ls_1, LinearSystem &ls_2)
  {
    std::swap(ls_1._A, ls_2._A);
    std::swap(ls_1._b, ls_2._b);
  }
};

/**
 * @brief Test whether two linear systems do not have common solutions
 *
 * @param ls1 is a linear system
 * @param ls2 is a linear system
 * @return `true` if and only if the solutions of `ls1` and `ls2` are
 *         disjoint
 */
bool have_disjoint_solutions(const LinearSystem &ls1, const LinearSystem &ls2);

/**
 * @brief Test whether two linear systems are the same one
 *
 * @param P1 is the first linear system
 * @param P2 is the second linear system
 * @return `true` if and only if the solution sets of two linear systems
 *         are the same set
 */
inline bool operator==(const LinearSystem &P1, const LinearSystem &P2)
{
  return P1.satisfies(P2) && P2.satisfies(P1);
}

/**
 * @brief Test whether two linear systems differ
 *
 * @param P1 is the first linear system
 * @param P2 is the second linear system
 * @return `true` if and only if the solution sets of two
 *         linear systems differ
 */
inline bool operator!=(const LinearSystem &P1, const LinearSystem &P2)
{
  return !(P1 == P2);
}

/**
 * @brief Print a linear system in a stream
 *
 * @param out is the output stream
 * @param ls is the linear system to be print
 * @return a reference to the output stream
 */
std::ostream &operator<<(std::ostream &out, const LinearSystem &ls);

#endif /* LINEARSYSTEM_H_ */
