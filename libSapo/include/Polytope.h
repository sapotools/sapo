/**
 * @file Polytope.h
 * @author <acasagrande@units.it>
 * @brief  Represent and manipulate polytopes
 * @version 0.1
 * @date 2021-11-16
 *
 * @copyright Copyright (c) 2021-2022
 */

#ifndef POLYTOPE_H_
#define POLYTOPE_H_

#include "LinearSystem.h"

/**
 * @brief Unions of closed sets
 *
 * This class represents unions of closed sets.
 *
 * @tparam BASIC_SET_TYPE is the basic set type
 */
template<class BASIC_SET_TYPE>
class SetsUnion;

/**
 * @brief A polytope representation class
 *
 * Polytopes are bounded convex sets representable
 * by a linear system \f[A \cdot x \leq b\f].
 * This class extends the `LinearSystem` class and
 * provides some set-based methods and functions
 * such as `is_empty` or `intersect`.
 */
class Polytope : public LinearSystem<double>
{

private:
  std::list<Polytope>
  split(const std::vector<bool> &bvect_base, const unsigned int cidx,
        std::list<Polytope> &tmp_covering, std::vector<std::vector<double>> &A,
        std::vector<double> &b, const unsigned int num_of_splits) const;

public:
  /**
   * Constructor that instantiates a unbounded polytope
   */
  Polytope(): LinearSystem<double>() {}

  /**
   * Copy constructor
   *
   * @param[in] orig the original polytope
   */
  Polytope(const Polytope &orig): LinearSystem<double>(orig) {}

  /**
   * Swap constructor
   *
   * @param[in] orig the original polytope
   */
  Polytope(Polytope &&orig): LinearSystem<double>(orig) {}

  /**
   * Constructor
   *
   * @param[in] A template matrix
   * @param[in] b offset vector
   * @todo Test polytope closeness
   */
  Polytope(const std::vector<std::vector<double>> &A,
           const std::vector<double> &b):
      LinearSystem<double>(A, b)
  {
  }

  /**
   * Move Constructor
   *
   * @param[in] A template matrix
   * @param[in] b offset vector
   */
  Polytope(std::vector<std::vector<double>> &&A, std::vector<double> &&b):
      LinearSystem<double>(A, b)
  {
  }

  /**
   * Constructor from a set of symbolic expressions
   *
   * @param[in] vars list of variables appearing in the constraints
   * @param[in] constraints symbolic constraints
   * @todo Test polytope closeness
   */
  Polytope(const std::vector<SymbolicAlgebra::Symbol<>> &vars,
           const std::vector<SymbolicAlgebra::Expression<>> &constraints):
      LinearSystem<double>(vars, constraints)
  {
  }

  /**
   * @brief Copy operator
   *
   * @param orig is the original Polytope
   * @return a reference to the updated object
   */
  Polytope &operator=(const Polytope &orig)
  {
    static_cast<LinearSystem<double> *>(this)->operator=(orig);

    return *this;
  }

  /**
   * @brief Copy operator
   *
   * @param orig is the original Polytope
   * @return a reference to the updated object
   */
  Polytope &operator=(Polytope &&orig)
  {
    static_cast<LinearSystem<double> *>(this)->operator=(orig);

    return *this;
  }

  /**
   * Establish whether a polytope is empty
   *
   * @param[in] strict_inequality specifies whether the polytope is
   *         defined by a strict inequality (i.e., Ax < b).
   * @return `true` when it can establish that the polytope is
   *         empty. `false` when it can establish that the polytope
   *         is not empty. `uncertain` in the other cases
   */
  inline TriBool is_empty(const bool strict_inequality = false) const
  {
    return !this->has_solutions(strict_inequality);
  }

  /**
   * Establish whether a polytope interior is empty
   *
   * @return `true` when it can establish that the polytope interior
   *         is empty. `false` when it can establish that the polytope
   *         interior is not empty. `uncertain` in the other cases
   */
  inline TriBool is_interior_empty() const
  {
    return !this->has_solutions(true);
  }

  /**
   * Check whether one polytope is a subset of another polytope
   *
   * This method establishes whether the current polytope is a
   * subset of another polytope.
   *
   * @param[in] P is the polytope that are compared to the current
   *             object
   * @return `true` if the method can establish that this polytope
   *         is a subset of `P`. `false` when it can establish that
   *         this polytope is not a subset of `P`. `uncertain` in
   *         the remaining cases
   */
  inline TriBool is_subset_of(const Polytope &P) const
  {
    return this->satisfies(P);
  }

  /**
   * Check whether one polytope is a superset of another polytope
   *
   * This method establishes whether the current polytope is a
   * superset of another polytope.
   *
   * @param[in] P is the polytope that are compared to the current
   *             object
   * @return `true` if the method can establish that this polytope
   *         is a superset of `P`. `false` when it can establish that
   *         this polytope is not a superset of `P`. `uncertain` in
   *         the remaining cases
   */
  inline TriBool includes(const Polytope &P) const
  {
    return P.satisfies(*this);
  }

  /**
   * @brief Expand the polytope
   *
   * This method expands the polytope so that each of its boundaries
   * is moved by a value `delta`.
   *
   * @param delta is the aimed expansion
   * @return a reference to the updated polytope
   */
  Polytope &expand_by(const double delta);

  /**
   *  Split a polytope in a list of polytopes.
   *
   *  This method splits a polytope in a list of polytopes such
   *  that their set union equals the original polytope.
   *
   *  @return A list of polytopes such that their union equals
   *      the original polytope.
   */
  inline std::list<Polytope> split() const
  {
    return this->split(this->_A.size());
  }

  /**
   *  Split a polytope in a list of polytopes.
   *
   *  This method splits a polytope in a list of polytopes such
   *  that their set union equals the original polytope. Each
   *  polytope in the original list is split in
   *  \f$2^\textrm{num_of_splits}\f$ polytopes at most.
   *
   * @param num_of_splits is the number of splits to performed.
   * @return A list of polytopes such that their union equals
   *      the original polytope.
   */
  std::list<Polytope> split(const unsigned int num_of_splits) const;

  /**
   * Update a polytope by intersecting it with another one.
   *
   * This method works in-place and changes the calling object.
   *
   * @param[in] P is a polytope.
   * @return a reference to the updated object.
   */
  Polytope &intersect_with(const Polytope &P);

  /**
   * Determine the volume of the bounding box of the polytope
   *
   * @return volume of the bounding box
   */
  double bounding_box_volume() const;

  /**
   * @brief Swap two polytopes
   *
   * This method swaps two polytope objects.
   *
   * @param P1 is a polytope
   * @param P2 is a polytope
   */
  friend inline void swap(Polytope &P1, Polytope &P2)
  {
    swap(*(static_cast<LinearSystem *>(&P1)),
         *(static_cast<LinearSystem *>(&P2)));
  }

  /**
   * Compute the intersection of two polytopes
   *
   * @param[in] P1 is a polytope
   * @param[in] P2 is a polytope
   * @return the intersection of the two parameters.
   */
  friend Polytope intersect(const Polytope &P1, const Polytope &P2);

  /**
   * Compute and over-approximation of the union of two polytopes
   *
   * @param[in] P1 is a polytope
   * @param[in] P2 is a polytope
   * @return an over-approximation of the union of `P1` and `P2`
   */
  friend Polytope over_approximate_union(const Polytope &P1,
                                         const Polytope &P2);

  /**
   * @brief Subtract two polytopes and close the result
   *
   * @param[in] P1 is a polytope
   * @param[in] P2 is a polytope
   * @return a union of polytopes obtained by closing the set
   *         \f$P1\setminus P2\f$
   */
  friend SetsUnion<Polytope> subtract_and_close(const Polytope &P1,
                                                const Polytope &P2);
};

/**
 * @brief Test whether two polytopes are disjoint
 *
 * @param P1 is a polytope
 * @param P2 is a polytope
 * @return `true` if and only if `P1` and `P2` are disjoint
 * @return `true` if the method can establish that `P1` and `P2`
 *         are disjoint. `false` when it can establish that
 *         `P1` and `P2` are not disjoint. `uncertain` in
 *         the remaining cases
 */
inline TriBool are_disjoint(const Polytope &P1, const Polytope &P2)
{
  return have_disjoint_solutions(P1, P2);
}

/**
 * @brief Get the expansion of a linear set
 *
 * This method expands a linear set so that each of its boundaries
 * is moved by a value `epsilon`.
 *
 * @tparam BASIC_SET_TYPE is the basic set type
 * @param S is the set to be expanded
 * @param epsilon is the aimed expansion
 * @return an expanded version of `S`
 */
template<class BASIC_SET_TYPE>
BASIC_SET_TYPE expand(const BASIC_SET_TYPE &S, const double epsilon)
{
  BASIC_SET_TYPE expS(S);

  expS.expand_by(epsilon);

  return expS;
}

/**
 * Compute and over-approximation of the union of two polytopes
 *
 * @param[in] P1 is a polytope
 * @param[in] P2 is a polytope
 * @return an over-approximation of the union of `P1` and `P2`
 */
Polytope over_approximate_union(const Polytope &P1, const Polytope &P2);

/**
 * @brief Subtract two polytopes and close the result
 *
 * @param[in] P1 is a polytope
 * @param[in] P2 is a polytope
 * @return a union of polytopes obtained by closing the set
 *         \f$P1\setminus P2\f$
 */
SetsUnion<Polytope> subtract_and_close(const Polytope &P1, const Polytope &P2);

#endif /* POLYTOPE_H_ */
