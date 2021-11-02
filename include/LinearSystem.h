/**
 * @file LinearSystem.h
 * Represent and manipulate a linear system
 * It can be used to represent polytopes (reached states, parameters, etc.)
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#ifndef LINEARSYSTEM_H_
#define LINEARSYSTEM_H_

#include <iostream>

#include <ginac/ginac.h>
#include <utility>
#include <vector>

class LinearSystemSet;

class LinearSystem
{

private:
  std::vector<std::vector<double>> A; // matrix A
  std::vector<double> b;              // vector b

  /**
   * Check if a constraint belongs to the linear system
   *
   * @param[in] Ai direction
   * @param[in] bi offset
   * @returns true is Ai x <= bi is in the linear system
   */
  bool isIn(std::vector<double> Ai, const double bi) const;

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
  bool satisfies(const std::vector<double> &Ai, const double bi) const;

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
  bool constraintIsRedundant(const unsigned int i) const;

  LinearSystemSet *get_a_finer_covering(const std::vector<bool> &bvect_base,
                                        const unsigned int cidx,
                                        LinearSystemSet *tmp_covering,
                                        std::vector<std::vector<double>> &A,
                                        std::vector<double> &b) const;

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
   * Swap constructor
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
  LinearSystem(const std::vector<std::vector<double>> &A,
               const std::vector<double> &b);

  /**
   * Constructor from a set of symbolic expressions
   *
   * @param[in] vars list of variables appearing in the constraints
   * @param[in] constraints symbolic constraints
   */
  LinearSystem(const GiNaC::lst &vars, const GiNaC::lst &constraints);

  /**
   * Return the template matrix
   *
   * @return template matrix
   */
  const std::vector<std::vector<double>> &getA() const
  {
    return this->A;
  }

  /**
   * Return the offset vector
   *
   * @return offset vector
   */
  const std::vector<double> &getb() const
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
   * Minimize the linear system
   *
   * @param[in] vars list of variables
   * @param[in] obj_fun objective function
   * @return minimum
   */
  double minLinearSystem(const GiNaC::lst &vars,
                         const GiNaC::ex &obj_fun) const;

  /**
   * Maximize the linear system
   *
   * @param[in] vars list of variables
   * @param[in] obj_fun objective function
   * @return maximum
   */
  double maxLinearSystem(const GiNaC::lst &vars,
                         const GiNaC::ex &obj_fun) const;

  /**
   * Minimize the linear system
   *
   * @param[in] obj_fun objective function
   * @return minimum
   */
  double minLinearSystem(const std::vector<double> &obj_fun_coeffs) const;

  /**
   * Maximize the linear system
   *
   * @param[in] obj_fun objective function
   * @return maximum
   */
  double maxLinearSystem(const std::vector<double> &obj_fun_coeffs) const;

  // testing methods

  /**
   * Establish whether a linear system has no solutions
   *
   * Due to approximation errors, it may return false for some empty
   * systems too. However, when it returns true, the set is certainly
   * empty.
   *
   * @param[in] strict_inequality specifies whether the linear system is
   *         a strict inequality (i.e., Ax < b).
   * @return a Boolean value. If the returned value is true, then the
   *       linear system is empty.
   */
  bool isEmpty(const bool strict_inequality = false) const;

  /**
   * Check whether all the solutions of a linear system are also solutions
   * for another linear system.
   *
   * This method establishes whether all the solutions of a linear system
   * are are also solutions for another linear system. Due to
   * approximation errors, it may return false even if this is the case.
   * However, whenever it returns true, the all the solutions of the
   * linear system are certainly solutions for the linear system passed
   * as parameter.
   *
   * @param[in] ls is the linear system whose solution set is compared to
   *     that of this linear system.
   * @return a Boolean value. When some of the solutions of this linear
   *     system are not solutions for the parameter, the returned value
   *     is false. When the method returns true, all the solution of the
   *     object are also solutions for the parameter. There are cases in
   *     which the set of object solutions is a subset of the parameter
   *     solutions and, still, this method returns false.
   */
  bool satisfies(const LinearSystem &ls) const;

  /**
   *  Split a linear system in set of linear systems.
   *
   *  This method splits a linear system in a set of linear systems such
   *  that the set union of their solutions equals the solution set of the
   *  original linear system.
   *
   *  @return A linear system set whose union of solution sets equals
   *      the solution set of the original linear system.
   */
  LinearSystemSet get_a_finer_covering() const;

  // operations on linear system
  /**
   * Update a linear system by joining the constraints of a linear system.
   *
   * This method works in-place and changes the calling object.
   *
   * @param[in] ls a linear system.
   * @return a reference to the updated object.
   */
  LinearSystem &intersectWith(const LinearSystem &ls);

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
   * Determine the volume of the bounding box of the linear system
   *
   * @return volume of the bounding box
   */
  double volBoundingBox();

  void plotRegion(std::ostream &os = std::cout, const char color = ' ') const;

  void plotRegionT(std::ostream &os, const double t) const;
  void plotRegion(std::ostream &os, const std::vector<int> &rows,
                  const std::vector<int> &cols) const;

  friend void swap(LinearSystem &ls_1, LinearSystem &ls_2);

  /**
   * Compute the intersection of two linear systems
   *
   * @param[in] A is a linear system
   * @param[in] B is a linear system
   * @return the linear system that represents the set of values
   * 	    satisfying both the parameters.
   */
  friend LinearSystem intersection(const LinearSystem &A,
                                   const LinearSystem &B);
};

inline void swap(LinearSystem &ls_1, LinearSystem &ls_2)
{
  std::swap(ls_1.A, ls_2.A);
  std::swap(ls_1.b, ls_2.b);
}

/**
 * Compute the complementary of a vector of values.
 *
 * @param[in] orig the vector of values to be complementated.
 * @return A vector containing the complementaries of the parameter values
 */
template<typename T>
std::vector<T> get_complementary(const std::vector<T> &orig)
{
  std::vector<T> res{orig};

  transform(res.begin(), res.end(), res.begin(), std::negate<T>());

  return res;
}

template<typename T>
std::ostream &operator<<(std::ostream &out, const std::vector<T> &v)
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
std::ostream &operator<<(std::ostream &out,
                         const std::vector<std::vector<T>> &A)
{
  for (auto row_it = std::begin(A); row_it != std::end(A); ++row_it) {
    if (row_it != std::begin(A)) {
      out << std::endl;
    }

    out << *row_it;
  }

  return out;
}

std::ostream &operator<<(std::ostream &out, const LinearSystem &ls);

#endif /* LINEARSYSTEM_H_ */
