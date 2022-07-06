/**
 * @file StickyUnion.h
 * @author Alberto Casagrande <acasagrande@units.it>
 * @brief Define and handle sticky unions
 * @version 0.1
 * @date 2022-06-02
 *
 * @copyright Copyright (c) 2022
 */

#ifndef STICKY_UNION_H_
#define STICKY_UNION_H_

#include <list>

#include "SetsUnion.h"

/**
 * @brief Sets unions over-approximating groups of intersecting sets
 *
 * The "sticky-intersecting" relation obtained by extending by
 * transitivity the "are intersecting" relation between sets is an
 * equivalence relation, \f$\approx_{\textrm{in}}\f$. The relation
 * \f$\approx_{\textrm{in}}\f$ partitions any collection of sets
 * into classes such that two sets \f$S_1\f$ and \f$S_2\f$ belong to
 * the same class if and only if there exists a subsets of original
 * collection \f$\{s_0, \ldots, s_k\}\f$ such that \f$S_0=s_0\f$,
 * \f$S_1=s_k\f$, and \f$s_i\f$ does intersect \f$s_{i+1}\f$ for all
 * \f$i \in [0, k-1]\f$.
 *
 * This C++ class represents union of sets in which each of the
 * \f$\approx_{\textrm{in}}\f$-classes is over-approximated by a
 * single basic set, e.g., bundle or polytope.
 *
 * This class is such that, for all unions of sets `U` and for all
 * sets `S`, if \f$S \subseteq U\f$, then
 * `StickyUnion(U).includes(S)`.
 *
 * @tparam BASIC_SET_TYPE is the set type
 */
template<class BASIC_SET_TYPE>
class StickyUnion
{
private:
  /**
   * @brief \f$\approx_{\textrm{in}}\f$-class type
   *
   * This type represents the "sticky-intersecting" classes. It
   * maintains both the union of sets in the class and the
   * over-approximation of the union to speed up the computation.
   */
  typedef struct {
    BASIC_SET_TYPE over_approx;           //!< over-approximation of the union
    SetsUnion<BASIC_SET_TYPE> sets_union; //!< sets in the class
  } class_type;

  std::list<class_type> _classes; //!< "sticky-intersecting" classes

public:
  /**
   * @brief Constructor
   */
  StickyUnion(): _classes() {}

  /**
   * @brief Copy constructor
   *
   * @tparam CONTAINER is the type of the basic set container
   * @param container is a container
   */
  template<template<class> class CONTAINER>
  StickyUnion(const CONTAINER<BASIC_SET_TYPE> &container): _classes()
  {
    for (auto it = std::begin(container); it != std::end(container); ++it) {
      add(*it);
    }
  }

  /**
   * @brief Copy constructor
   *
   * @param set_obj is a set
   */
  StickyUnion(const BASIC_SET_TYPE &set_obj): _classes()
  {
    add(set_obj);
  }

  /**
   * @brief Copy constructor
   *
   * @param container is a container
   */
  StickyUnion(const std::list<BASIC_SET_TYPE> &container): _classes()
  {
    for (auto it = std::begin(container); it != std::end(container); ++it) {
      add(*it);
    }
  }

  /**
   * @brief Copy constructor
   *
   * @param orig is a sticky union
   */
  StickyUnion(const StickyUnion<BASIC_SET_TYPE> &orig)
  {
    _classes = orig._classes;
  }

  /**
   * @brief A swap constructor for sticky unions
   *
   * @param[in] orig is a sticky union
   */
  StickyUnion(StickyUnion<BASIC_SET_TYPE> &&orig)
  {
    std::swap(_classes, orig._classes);
  }

  /**
   * @brief An assignment operator for sticky unions
   *
   * @param[in] orig is a sticky unions
   */
  StickyUnion<BASIC_SET_TYPE> &
  operator=(const StickyUnion<BASIC_SET_TYPE> &orig)
  {
    _classes = orig._classes;

    return *this;
  }

  /**
   * @brief A swap assignment for sticky unions
   *
   * @param[in] orig is a sticky unions
   */
  StickyUnion<BASIC_SET_TYPE> &operator=(StickyUnion<BASIC_SET_TYPE> &&orig)
  {
    std::swap(_classes, orig._classes);

    return *this;
  }

  /**
   * @brief Add a new set to the union
   *
   * This method adds a set to the union, computes the new
   * "sticky-intersecting" partition, and over-approximates
   * its classes by using a single basic set.
   *
   * The number of set operation required by this method is
   * \f$O(n_s)\f$ where \f$n_s\f$ is the number of the
   * original exact sets in the sticky union.
   *
   * @param set_obj is the set to be added
   * @return `true` if and only if `set_obj` was not already
   *         contained in the object
   */
  bool add(const BASIC_SET_TYPE &set_obj)
  {
    if (set_obj.is_empty()) {
      return false;
    }

    // `jc_union` will contain the sets union
    // of the joint classes
    SetsUnion<BASIC_SET_TYPE> jc_union;

    // `jc_approx` is devoted to the over-approximation
    // of `jc_union`
    BASIC_SET_TYPE jc_approx = set_obj;

    // for each class
    for (auto c_it = std::begin(_classes); c_it != std::end(_classes);) {

      // if it is not disjoint
      if (!are_disjoint(c_it->sets_union, set_obj)) {

        // reverse the class sets union in `jc_union`.
        // Since all the classes are disjoint by definition
        // checking for inclusion as done by `SetsUnion::add`
        // can be avoided
        jc_union.splice(jc_union.end(), c_it->sets_union);

        // update `jc_approx`
        jc_approx = over_approximate_union(c_it->over_approx, jc_approx);

        // remove the class
        c_it = _classes.erase(c_it);
      } else {

        // if it is disjoint, go to the next class
        ++c_it;
      }
    }

    // add the new set to the sets union and check for
    // changed in the union itself
    bool changed = jc_union.add(set_obj);

    // add the joint class to the partition
    _classes.push_back({std::move(jc_approx), std::move(jc_union)});

    return changed;
  }

  /**
   * @brief Update a sticky union by joining a sets union
   *
   * This method works in-place and changes the calling object.
   *
   * @param[in] sets_union is a sets union
   * @return a reference to the updated object
   */
  bool update(const SetsUnion<BASIC_SET_TYPE> &sets_union)
  {
    bool changed = false;

    for (auto it = std::cbegin(sets_union); it != std::cend(sets_union);
         ++it) {
      changed = changed || add(*it);
    }

    return changed;
  }

  /**
   * @brief Test whether any of the sets in the union includes a set
   *
   * @param set_obj is the set whose inclusion has to be tested
   * @return `true` if and only if any of the set in the sticky union
   *     includes `set_obj`
   */
  bool any_includes(const BASIC_SET_TYPE &set_obj) const
  {
    if (set_obj.is_empty()) {
      return true;
    }

    // for each class
    for (auto c_it = std::begin(_classes); c_it != std::end(_classes);
         ++c_it) {

      // if the class over-approximation includes
      if (c_it->sets_union.any_includes(set_obj)) {
        return true;
      }
    }

    return false;
  }

  /**
   * @brief Test whether a sticky union includes a sets union
   *
   * @param sets_union is the sets union whose inclusion has to be tested
   * @return `true` if and only if for any set in `sets_union` there exists a
   *     set in the current object that includes it
   */
  bool any_includes(const SetsUnion<BASIC_SET_TYPE> &sets_union) const
  {
    for (auto s_it = std::begin(sets_union); s_it != std::end(sets_union);
         ++s_it) {
      if (!any_includes(*s_it)) {
        return false;
      }
    }

    return true;
  }

  /**
   * @brief Test whether the over-approximation of the classes includes a set
   *
   * @param set_obj is the set whose inclusion has to be tested
   * @return `true` if and only if the over-approximation of the union of
   *  one of the \f$\approx_{\textrm{in}}\f$-classes includes `set_obj`
   */
  bool includes(const BASIC_SET_TYPE &set_obj) const
  {
    if (set_obj.is_empty()) {
      return true;
    }

    // for each class
    for (auto c_it = std::begin(_classes); c_it != std::end(_classes);
         ++c_it) {

      // if the class over-approximation includes
      if (c_it->over_approx.includes(set_obj)) {
        return true;
      }
    }

    return false;
  }

  /**
   * @brief Test whether the over-approximation of the classes includes a sets
   * union
   *
   * @param sets_union is the sets union whose inclusion has to be tested
   * @return `true` if and only if the over-approximation of the union of
   *  one of the \f$\approx_{\textrm{in}}\f$-classes includes `sets_union`
   */
  bool includes(const SetsUnion<BASIC_SET_TYPE> &sets_union) const
  {
    for (auto s_it = std::begin(sets_union); s_it != std::end(sets_union);
         ++s_it) {
      if (!includes(*s_it)) {
        return false;
      }
    }

    return true;
  }

  /**
   * @brief Check that the over-approximation of a sticky union satisfies a
   * linear system
   *
   * This method checks whether the over-approximations of the classes of a
   * sticky union satisfy a linear system
   *
   * @param ls is the considered linear system
   * @return `true` if and only if all the points in the over-approximation of
   *          the classes of the current object are solutions for `ls`
   */
  bool satisfies(const LinearSystem &ls) const
  {
    if (ls.size() == 0) {
      return true;
    }

    for (auto c_it = std::begin(_classes); c_it != std::end(_classes);
         ++c_it) {
      if (!c_it->over_approx.satisfies(ls)) {
        return false;
      }
    }

    return true;
  }

  /**
   * @brief Check that a sets in a sticky union satisfy a linear system
   *
   * This method checks whether all the sets in a sticky union satisfy a
   * linear system.
   *
   * @param ls is the considered linear system
   * @return `true` if and only if all the points in the basic sets of the
   *         current object are solutions for `ls`
   */
  bool exact_satisfies(const LinearSystem &ls) const
  {
    if (ls.size() == 0) {
      return true;
    }

    for (auto c_it = std::begin(_classes); c_it != std::end(_classes);
         ++c_it) {
      if (!c_it->sets_union.satisfies(ls)) {
        return false;
      }
    }

    return true;
  }

  /**
   * @brief Convert a sticky union in the corresponding sets union
   *
   * @return the sets union that produced the current object
   */
  operator SetsUnion<BASIC_SET_TYPE>() const
  {
    SetsUnion<BASIC_SET_TYPE> sets_union;

    // for each class
    for (auto c_it = std::begin(_classes); c_it != std::end(_classes);
         ++c_it) {
      // Since all the classes are disjoint by definition
      // checking for inclusion as done by `SetsUnion::add`
      // can be avoided
      for (auto s_it = std::begin(c_it->sets_union);
           s_it != std::end(c_it->sets_union); ++s_it) {
        sets_union.push_back(*s_it);
      }
    }

    return sets_union;
  }

  /**
   * @brief Compute the intersection between two sticky unions
   *
   * @tparam SET_TYPE is a set type, e.g., SetsUnion or Bundle
   * @param sticky_union is a sticky union
   * @return the intersection of the current object and `sticky_union`
   */
  template<class SET_TYPE>
  SetsUnion<BASIC_SET_TYPE>
  get_intersection_with(const StickyUnion<BASIC_SET_TYPE> &sticky_union) const
  {
    SetsUnion<BASIC_SET_TYPE> sets_union;

    // for each class
    for (auto c_it = std::begin(_classes); c_it != std::end(_classes);
         ++c_it) {

      // append the intersection of `sticky_union` and the over-approximation
      // of the class pointed by `c_it`
      sets_union.update(sticky_union.get_intersection_with(c_it->over_approx));
    }

    return sets_union;
  }

  /**
   * @brief Compute the intersection of a sticky union and a sets union
   *
   * @tparam SET_TYPE is a set type, e.g., SetsUnion or Bundle
   * @param set_obj is a set
   * @return the intersection of the current object and `set_obj`
   */
  template<class SET_TYPE>
  SetsUnion<BASIC_SET_TYPE>
  get_intersection_with(const SET_TYPE &set_obj) const
  {
    SetsUnion<BASIC_SET_TYPE> sets_union;

    // for each class
    for (auto c_it = std::begin(_classes); c_it != std::end(_classes);
         ++c_it) {

      // append the intersection of `set_obj` and the over-approximation of the
      // class pointed by `c_it`
      sets_union.add(intersect(c_it->over_approx, set_obj));
    }

    return sets_union;
  }

  /**
   * @brief Test whether the sticky union is empty
   *
   * @return `true` if and only if the sticky union is empty
   */
  inline bool is_empty() const
  {
    return _classes.size() == 0;
  }

  /**
   * @brief Get the number of sets in the sticky union
   *
   * @returns the number of sets in the sticky union
   */
  size_t size() const
  {
    size_t res = 0;

    // for each class
    for (auto c_it = std::begin(_classes); c_it != std::end(_classes);
         ++c_it) {

      res += c_it->sets_union.size();
    }

    return res;
  }

  /**
   * @brief Get the number of classes in the
   * \f$\approx_{\textrm{in}}\f$-partition
   *
   * @returns the number of classes in the
   * \f$\approx_{\textrm{in}}\f$-partition
   */
  inline size_t number_of_classes() const
  {
    return _classes.size();
  }

  /**
   * @brief Get the space dimension of the sets
   *
   * @returns the space dimension of the sets or
   *          0 if the union is the empty
   */
  inline size_t dim() const
  {
    if (is_empty()) {
      return 0;
    }

    return _classes.front().over_approx.dim();
  }

  /**
   * @brief Compute the intersection of two sticky unions
   *
   * @tparam BASIC_SET_TYPE2 is the basic set type, e.g., Polytope or
   * Bundle
   * @param A is a sticky union of sets
   * @param B is a sticky union of sets
   * @return the intersection of `A` and `B`
   */
  template<class BASIC_SET_TYPE2>
  friend SetsUnion<BASIC_SET_TYPE2>
  intersect(const StickyUnion<BASIC_SET_TYPE2> &A,
            const StickyUnion<BASIC_SET_TYPE2> &B);

  /**
   * @brief Compute the intersection of a set and a sticky union
   *
   * @tparam BASIC_SET_TYPE2 is the basic set type, e.g., Polytope or
   * Bundle
   * @tparam SET_TYPE is a set type, e.g., SetsUnion or Bundle
   * @param sticky_union is a sticky union of sets
   * @param set_obj is a set
   * @return the intersection of the two parameters
   */
  template<class BASIC_SET_TYPE2, class SET_TYPE>
  friend SetsUnion<BASIC_SET_TYPE2>
  intersect(const StickyUnion<BASIC_SET_TYPE2> &sticky_union,
            const SET_TYPE &set_obj);

  /**
   * @brief Compute the intersection of a sticky union and a set
   *
   * @tparam BASIC_SET_TYPE is the basic set type, e.g., Polytope or
   * Bundle
   * @tparam SET_TYPE is a set type, e.g., SetsUnion or Bundle
   * @param sticky_union is a sticky union of sets
   * @param set_obj is a set
   * @return the intersection of the two parameters
   */
  template<class BASIC_SET_TYPE2, class SET_TYPE>
  friend SetsUnion<BASIC_SET_TYPE2>
  intersect(const StickyUnion<BASIC_SET_TYPE2> &sticky_union,
            const BASIC_SET_TYPE2 &set_obj);
};

/**
 * @brief Compute the intersection of two sticky unions
 *
 * @tparam BASIC_SET_TYPE is the basic set type, e.g., Polytope or
 * Bundle
 * @param A is a sticky union of sets
 * @param B is a sticky union of sets
 * @return the intersection of `A` and `B`
 */
template<class BASIC_SET_TYPE>
SetsUnion<BASIC_SET_TYPE> intersect(const StickyUnion<BASIC_SET_TYPE> &A,
                                    const StickyUnion<BASIC_SET_TYPE> &B)
{
  SetsUnion<BASIC_SET_TYPE> sets_union;

  // for each class
  for (auto c_it = std::begin(A._classes); c_it != std::end(A._classes);
       ++c_it) {

    sets_union.update(intersect(B, c_it->over_approx));
  }

  return sets_union;
}

/**
 * @brief Compute the intersection of a sticky union and a set
 *
 * @tparam BASIC_SET_TYPE is the basic set type, e.g., Polytope or
 * Bundle
 * @tparam SET_TYPE is a set type, e.g., SetsUnion or Bundle
 * @param sticky_union is a sticky union of sets
 * @param set_obj is a set
 * @return the intersection of the two parameters
 */
template<class BASIC_SET_TYPE, class SET_TYPE>
SetsUnion<BASIC_SET_TYPE>
intersect(const StickyUnion<BASIC_SET_TYPE> &sticky_union,
          const SET_TYPE &set_obj)
{
  SetsUnion<BASIC_SET_TYPE> sets_union;

  // for each class
  for (auto c_it = std::begin(sticky_union._classes);
       c_it != std::end(sticky_union._classes); ++c_it) {

    // append the intersection of `set_obj` and the over-approximation of the
    // class pointed by `c_it`
    sets_union.update(intersect(c_it->over_approx, set_obj));
  }

  return sets_union;
}

/**
 * @brief Compute the intersection of a set and a sticky union
 *
 * @tparam BASIC_SET_TYPE is the basic set type, e.g., Polytope or
 * Bundle
 * @tparam SET_TYPE is a set type, e.g., SetsUnion or Bundle
 * @param set_obj is a sets union
 * @param sticky_union is a sticky union of sets
 * @return the intersection of the two parameters
 */
template<class BASIC_SET_TYPE, class SET_TYPE>
inline SetsUnion<BASIC_SET_TYPE>
intersect(const SET_TYPE &set_obj,
          const StickyUnion<BASIC_SET_TYPE> &sticky_union)
{
  return intersect(sticky_union, set_obj);
}

/**
 * @brief Compute the intersection of a sticky union and a set
 *
 * @tparam BASIC_SET_TYPE is the basic set type, e.g., Polytope or
 * Bundle
 * @param sticky_union is a sticky union of sets
 * @param set_obj is a set
 * @return the intersection of the two parameters
 */
template<class BASIC_SET_TYPE, class SET_TYPE>
SetsUnion<BASIC_SET_TYPE>
intersect(const StickyUnion<BASIC_SET_TYPE> &sticky_union,
          const BASIC_SET_TYPE &set_obj)
{
  SetsUnion<BASIC_SET_TYPE> sets_union;

  // for each class
  for (auto c_it = std::begin(sticky_union._classes);
       c_it != std::end(sticky_union._classes); ++c_it) {

    // append the intersection of `set_obj` and the over-approximation of the
    // class pointed by `c_it`
    sets_union.add(intersect(c_it->over_approx, set_obj));
  }

  return sets_union;
}

/**
 * @brief Compute the intersection of a set and a sticky union
 *
 * @tparam BASIC_SET_TYPE is the basic set type, e.g., Polytope or
 * Bundle
 * @tparam SET_TYPE is a set type, e.g., SetsUnion or Bundle
 * @param set_obj is a sets union
 * @param sticky_union is a sticky union of sets
 * @return the intersection of the two parameters
 */
template<class BASIC_SET_TYPE, class SET_TYPE>
inline SetsUnion<BASIC_SET_TYPE>
intersect(const BASIC_SET_TYPE &set_obj,
          const StickyUnion<BASIC_SET_TYPE> &sticky_union)
{
  return intersect(sticky_union, set_obj);
}

#endif // STICKY_UNION_H_