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
 * @brief A polytope representation class
 *
 * Polytopes are bounded convex sets representable
 * by a linear system \f[A \cdot x \leq b\f].
 * This class extends the `LinearSystem` class and
 * provides some set-based methods and functions
 * such as `is_empty` or `intersect`.
 */
class Polytope : public LinearSystem
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
  Polytope(): LinearSystem() {}

  /**
   * Copy constructor
   *
   * @param[in] orig the original polytope
   */
  Polytope(const Polytope &orig): LinearSystem(orig) {}

  /**
   * Swap constructor
   *
   * @param[in] orig the original polytope
   */
  Polytope(Polytope &&orig): LinearSystem(orig) {}

  /**
   * Constructor
   *
   * @param[in] A template matrix
   * @param[in] b offset vector
   * @todo Test polytope closeness
   */
  Polytope(const std::vector<std::vector<double>> &A,
           const std::vector<double> &b):
      LinearSystem(A, b)
  {
  }

  /**
   * Move Constructor
   *
   * @param[in] A template matrix
   * @param[in] b offset vector
   */
  Polytope(std::vector<std::vector<double>> &&A, std::vector<double> &&b):
      LinearSystem(A, b)
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
      LinearSystem(vars, constraints)
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
    static_cast<LinearSystem *>(this)->operator=(orig);

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
    static_cast<LinearSystem *>(this)->operator=(orig);

    return *this;
  }

  /**
   * Establish whether a polytope is empty
   *
   * @param[in] strict_inequality specifies whether the polytope is
   *         defined by a strict inequality (i.e., Ax < b).
   * @return `true` if and only if the polytope is empty
   */
  bool is_empty(const bool strict_inequality = false) const
  {
    return !this->has_solutions(strict_inequality);
  }

  /**
   * Check whether one polytope is a subset of another polytope
   *
   * This method establishes whether the current polytope is a
   * subset of another polytope.
   *
   * @param[in] P is the polytope that are compared to the current
   *             object
   * @return `true` if and only if this polytope is a subset of `P`
   */
  inline bool is_subset_of(const Polytope &P) const
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
   * @return `true` if and only if this polytope is a superset of `P`
   */
  inline bool includes(const Polytope &P) const
  {
    return P.satisfies(*this);
  }

  /**
   * @brief Expand the polytope
   *
   * This method expands the polytope so that each of its boundaries
   * is moved by a value `epsilon`.
   *
   * @param epsilon is the aimed expansion
   * @return a reference to the updated polytope
   */
  Polytope &expand_by(const double epsilon);

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
};

/**
 * @brief Test whether two polytopes are disjoint
 *
 * @param P1 is a polytope
 * @param P2 is a polytope
 * @return `true` if and only if `P1` and `P2` are disjoint
 */
inline bool are_disjoint(const Polytope &P1, const Polytope &P2)
{
  return disjoint_solutions(P1, P2);
}

/**
 * @brief Get the expansion of a linear set
 *
 * This method expands a linear set so that each of its boundaries
 * is moved by a value `epsilon`.
 *
 * @tparam SET_TYPE is the linear set type
 * @param S is the set to be expanded
 * @param epsilon is the aimed expansion
 * @return an expanded version of `S`
 */
template<typename SET_TYPE>
SET_TYPE expand(const SET_TYPE &S, const double epsilon)
{
  SET_TYPE expS(S);

  expS.expand_by(epsilon);

  return expS;
}

Polytope over_approximate_union(const Polytope &P1, const Polytope &P2);

#endif /* POLYTOPE_H_ */
