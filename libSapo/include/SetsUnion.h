/**
 * @file SetsUnion.h
 * @author Alberto Casagrande <acasagrande@units.it>
 * @brief Representing and handling unions of closed sets
 * @version 0.1
 * @date 2022-05-27
 *
 * @copyright Copyright (c) 2022
 */

#ifndef SETSUNION_H_
#define SETSUNION_H_

#include <list>
#include <memory>

#ifdef WITH_THREADS
#include <mutex>
#include <shared_mutex>

#include "SapoThreads.h"
#endif // WITH_THREADS

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
 * @brief Unions of closed sets
 *
 * This template represents unions of closed sets as a
 * list of non-comparable sets: any two sets \f$A\f$ and
 * \f$B\f$ in the list are such that
 * \f$A \not\subseteq B\f$ and \f$A \not\supseteq B\f$.
 * Whenever a new sets \f$S\f$ is added to the union,
 * if \f$S\f$ is a subset for any set in the list,
 * the list does not change. Otherwise, all the sets in
 * the list included by \f$S\f$ are removed from the
 * list and \f$S\f$ itself is appended at its end.
 *
 * @tparam BASIC_SET_TYPE is the basic set type
 */
template<class BASIC_SET_TYPE>
class SetsUnion : private std::list<BASIC_SET_TYPE>
{
private:
  /**
   * @brief Append a list of closed sets at the end of the current sets union
   *
   * @param list is the list of closed sets to append
   */
  inline void append(std::list<BASIC_SET_TYPE> &list)
  {
    this->splice(end(), list);
  }

  /**
   * @brief Append a list of closed sets at the end of the current sets union
   *
   * @param list is the list of closed sets to append
   */
  inline void append(std::list<BASIC_SET_TYPE> &&list)
  {
    this->splice(end(), std::move(list));
  }

  /**
   * @brief Add a set to a sets union
   *
   * This method adds a set `set_obj` to a union of closed sets.
   * First of all the method, tests whether the new set
   * is empty. If this is the case, no set is added to the
   * union. No set is added to the union also in the case
   * in which `set_obj` is subset of some of the first
   * `sets_to_cmp` elements in the list of the united sets.
   * Otherwise, all the first `sets_to_cmp` elements in the
   * list of the united sets that are subsets of `set_obj`
   * are removed from the list and `set_obj` is added
   * to the set.
   * The method returns a Boolean value that testifies
   * whether the computation has added `set_obj` at the
   * end of the sets list.
   *
   * @param[in] set_obj is the set to be added
   * @param[in] sets_to_cmp is the number of closed sets in the list
   *            to be compared to `set_obj`
   * @return `true` if and only if the computation has added
   *         `set_obj` at the end of the sets list
   */
  bool add(const BASIC_SET_TYPE &set_obj, size_t sets_to_cmp)
  {
    if (size() != 0 && (this->front().dim() != set_obj.dim())) {
      SAPO_ERROR("adding a set to a union of closed sets "
                 "that has different dimension",
                 std::domain_error);
    }

    if (set_obj.is_empty()) {
      return false;
    }

    // for any of the first `sets_to_cmp` sets in the list
    auto it = std::begin(*this);
    while (it != std::end(*this) && sets_to_cmp-- > 0) {
      // if the set includes `set_obj`, then
      // `set_obj` is already included in the union
      // and the current object can be returned
      if (it->includes(set_obj)) {
        return false;
      }

      // if the set is a subset of `set_obj`, then
      // this set can be removed as `set_obj` is
      // going to be added to the list
      if (it->is_subset_of(set_obj)) {
        it = std::list<BASIC_SET_TYPE>::erase(it);
      } else {
        ++it;
      }
    }

    // if `set_obj` is not a subset of any
    // of the sets in the head, then
    // append it to the union itself
    this->push_back(set_obj);

    return true;
  }

  /**
   * @brief Add a set to a sets union
   *
   * This method adds a set `set_obj` to a union of closed sets.
   * First of all the method, tests whether the new set
   * is empty. If this is the case, no set is added to the
   * union. No set is added to the union also in the case
   * in which `set_obj` is subset of some of the first
   * `sets_to_cmp` elements in the list of the united sets.
   * Otherwise, all the first `sets_to_cmp` elements in the
   * list of the united sets that are subsets of `set_obj`
   * are removed from the list and `set_obj` is added
   * to the set.
   * The method returns a Boolean value that testifies
   * whether the computation has added `set_obj` at the
   * end of the sets list.
   *
   * @param[in] set_obj is the set to be added
   * @param[in] sets_to_cmp is the number of closed sets in the list
   *            to be compared to `set_obj`
   * @return `true` if and only if `set_obj` was not already
   *         contained in the object
   */
  bool add(BASIC_SET_TYPE &&set_obj, size_t sets_to_cmp)
  {
    if (size() != 0 && (this->front().dim() != set_obj.dim())) {
      SAPO_ERROR("adding a set to a union of closed sets "
                 "that has different dimension",
                 std::domain_error);
    }

    if (set_obj.is_empty()) {
      return false;
    }

    // for any of the first `sets_to_cmp` sets in the list
    auto it = std::begin(*this);
    while (it != std::end(*this) && sets_to_cmp-- > 0) {
      // if the set includes `set_obj`, then
      // `set_obj` is already included in the union
      // and the current object can be returned
      if (it->includes(set_obj)) {
        return false;
      }

      // if the set is a subset of `set_obj`, then
      // this set can be removed as `set_obj` is
      // going to be added to the list
      if (it->is_subset_of(set_obj)) {
        it = std::list<BASIC_SET_TYPE>::erase(it);
      } else {
        ++it;
      }
    }

    // if `set_obj` is not a subset of any
    // of the sets in the head, then
    // append it to the union itself
    this->push_back(std::move(set_obj));

    return true;
  }

public:
  /**
   * @brief A constant iterator alias
   */
  using const_iterator = typename std::list<BASIC_SET_TYPE>::const_iterator;

  /**
   * @brief A iterator alias
   */
  using iterator = typename std::list<BASIC_SET_TYPE>::iterator;

  /**
   * @brief Constructor
   *
   */
  SetsUnion() {}

  /**
   * @brief Constructor
   *
   * @param[in] set_obj is the set to be included in
   *            the union
   */
  SetsUnion(const BASIC_SET_TYPE &set_obj)
  {
    this->add(set_obj);
  }

  /**
   * @brief Constructor
   *
   * @param[in] set_obj is the set to be included in
   *            the union
   */
  SetsUnion(BASIC_SET_TYPE &&set_obj)
  {
    this->add(std::move(set_obj));
  }

  /**
   * @brief A copy constructor for a union of closed sets
   *
   * @param[in] orig is a union of closed sets
   */
  SetsUnion(const SetsUnion<BASIC_SET_TYPE> &orig):
      std::list<BASIC_SET_TYPE>(orig)
  {
  }

  /**
   * @brief Constructor
   *
   * @param[in] sets is a container of closed sets
   */
  SetsUnion(const std::list<BASIC_SET_TYPE> &sets)
  {
    for (auto it = std::begin(sets); it != std::end(sets); ++it) {
      this->add(*it);
    }
  }

  /**
   * @brief Constructor
   *
   * @param[in] sets is a container of closed sets
   */
  SetsUnion(std::list<BASIC_SET_TYPE> &&sets)
  {
    for (auto it = std::begin(sets); it != std::end(sets); ++it) {
      this->add(std::move(*it));
    }
  }

  /**
   * @brief A swap constructor for a union of closed sets
   *
   * @param[in] orig is a union of closed sets
   */
  SetsUnion(SetsUnion<BASIC_SET_TYPE> &&orig)
  {
    std::swap(*(static_cast<std::list<BASIC_SET_TYPE> *>(this)),
              *(static_cast<std::list<BASIC_SET_TYPE> *>(&orig)));
  }

  /**
   * @brief An assignment operator for unions of closed sets
   *
   * @param[in] orig is a union of closed sets
   */
  SetsUnion<BASIC_SET_TYPE> &operator=(const SetsUnion<BASIC_SET_TYPE> &orig)
  {
    this->clear();

    for (auto it = std::cbegin(orig); it != std::cend(orig); ++it) {
      this->push_back(*it);
    }

    return *this;
  }

  /**
   * @brief A swap assignment for unions of closed sets
   *
   * @param[in] orig is a union of closed sets
   */
  SetsUnion<BASIC_SET_TYPE> &operator=(SetsUnion<BASIC_SET_TYPE> &&orig)
  {
    std::swap(*(static_cast<std::list<BASIC_SET_TYPE> *>(this)),
              *(static_cast<std::list<BASIC_SET_TYPE> *>(&orig)));

    return *this;
  }

  /**
   * @brief Get the first set in the union
   *
   * @return A reference to the first set in the union
   */
  inline BASIC_SET_TYPE &front()
  {
    return static_cast<std::list<BASIC_SET_TYPE> *>(this)->front();
  }

  /**
   * @brief Get the first set in the union
   *
   * @return A constant reference to the first set in the union
   */
  inline const BASIC_SET_TYPE &front() const
  {
    return static_cast<const std::list<BASIC_SET_TYPE> *>(this)->front();
  }

  /**
   * @brief Add a set to the union
   *
   * This method adds a set `set_obj` to a union of closed sets.
   * First of all the method, tests whether the new set
   * is empty. If this is the case, no set is added to the
   * union. No set is added to the union also in the case
   * in which `set_obj` is subset of some of the elements
   * in the list of the united sets.
   * Otherwise, the elements in the list of the united
   * sets that are subsets of `set_obj` are removed from
   * the list and `set_obj` is added to the set.
   *
   * @param[in] set_obj is the set to be added
   * @return `true` if and only if `set_obj` was not already
   *         contained in the object
   */
  inline bool add(const BASIC_SET_TYPE &set_obj)
  {
    return add(set_obj, size());
  }

  /**
   * @brief Add a set to the union
   *
   * This method adds a set `set_obj` to a union of closed sets.
   * First of all the method, tests whether the new set
   * is empty. If this is the case, no set is added to the
   * union. No set is added to the union also in the case
   * in which `set_obj` is subset of some of the elements
   * in the list of the united sets.
   * Otherwise, the elements in the list of the united
   * sets that are subsets of `set_obj` are removed from
   * the list and `set_obj` is added to the set.
   *
   * @param[in] set_obj is the set to be added
   * @return `true` if and only if `set_obj` was not already
   *         contained in the object
   */
  inline bool add(BASIC_SET_TYPE &&set_obj)
  {
    return add(std::move(set_obj), size());
  }

  /**
   * @brief Get the last set that changed the union
   *
   * @return the last set added to the the union that
   * changed it
   */
  inline const BASIC_SET_TYPE &last_added_set() const
  {
    return std::list<BASIC_SET_TYPE>::back();
  }

  /**
   * @brief Update a sets union by joining another sets union
   *
   * This method works in-place and changes the calling object.
   *
   * @param[in] sets_union is a sets union
   * @return a reference to the updated object
   */
  SetsUnion<BASIC_SET_TYPE> &
  update(const SetsUnion<BASIC_SET_TYPE> &sets_union)
  {
    unsigned int from_sets_union = 0;

    for (auto it = std::cbegin(sets_union); it != std::cend(sets_union);
         ++it) {
      if (this->add(*it, size() - from_sets_union)) {
        ++from_sets_union;
      }
    }

    return *this;
  }

  /**
   * @brief Update a sets union by joining another sets union
   *
   * This method works in-place and changes the calling object.
   *
   * @param[in] sets_union is a sets union
   * @return a reference to the updated object
   */
  SetsUnion<BASIC_SET_TYPE> &update(SetsUnion<BASIC_SET_TYPE> &&sets_union)
  {
    unsigned int from_sets_union = 0;

    for (auto it = std::begin(sets_union); it != std::end(sets_union); ++it) {
      if (this->add(std::move(*it), size() - from_sets_union)) {
        ++from_sets_union;
      }
    }

    return *this;
  }

  /**
   * @brief Check whether a sets union satisfies a linear system
   *
   * This method checks whether all the points in the
   * current object are solutions for a linear system.
   *
   * @param ls is the considered linear system
   * @return `true` if and only if all the points of
   *          the current union of closed sets are solutions
   *          for `ls`
   */
  bool satisfies(const LinearSystem &ls) const
  {
    for (auto it = begin(); it != end(); ++it) {
      if (!it->satisfies(ls)) {
        return false;
      }
    }

    return true;
  }

  /**
   * @brief Test whether a sets union is subset of a set
   *
   * This method tests whether the current object is subset
   * for a set.
   *
   * @tparam BASIC_SET_TYPE2 is the type of the tested set
   * @param[in] set_obj is the set
   * @return `true` if and only if the current bundle is a
   *         subset of `set_obj`
   */
  template<class BASIC_SET_TYPE2>
  bool is_subset_of(const BASIC_SET_TYPE2 &set_obj) const
  {
    for (auto it = begin(); it != end(); ++it) {
      if (!it->is_subset_of(set_obj)) {
        return false;
      }
    }

    return true;
  }

  /**
   * @brief Test whether a sets union is subset of a sets union
   *
   * @tparam BASIC_SET_TYPE2 is the type of the tested set
   * @param sets_union is the tested sets union
   * @return `true` if and only if the current bundle is a
   *         subset of `sets_union`
   */
  template<class BASIC_SET_TYPE2>
  inline bool is_subset_of(const SetsUnion<BASIC_SET_TYPE2> &sets_union) const
  {
    return sets_union.includes(*this);
  }

  /**
   * @brief Check whether a sets union includes a closed set
   *
   * @tparam BASIC_SET_TYPE2 is the type of the tested set
   * @param set_obj is the set whose inclusion must be tested
   * @return `true` if and only if `set_obj` is a subset of the
   *         current object
   */
  template<class BASIC_SET_TYPE2>
  inline bool includes(const BASIC_SET_TYPE2 &set_obj) const
  {
    return subtract_and_close(set_obj, *this).is_empty();
  }

  /**
   * @brief Check whether a sets union includes another sets union
   *
   * @tparam BASIC_SET_TYPE2 is the type of the tested sets union
   * @param sets_union is the set whose inclusion must be tested
   * @return `true` if and only if `sets_union` is a subset of the
   *         current object
   */
  template<class BASIC_SET_TYPE2>
  inline bool includes(const SetsUnion<BASIC_SET_TYPE2> &sets_union) const
  {
    return subtract_and_close(sets_union, *this).is_empty();
  }

  /**
   * @brief Check whether one of the sets in a union includes a set
   *
   * @tparam BASIC_SET_TYPE2 is the type of the tested set
   * @param[in] set_obj is the set whose inclusion must be tested
   * @return `true` if and only if `set_obj` is a subset of some of
   *         the sets in the union
   */
  template<class BASIC_SET_TYPE2>
  bool any_includes(const BASIC_SET_TYPE2 &set_obj) const
  {

    if (set_obj.is_empty()) {
      return true;
    }

#ifdef WITH_THREADS
    class ThreadResult
    {
      mutable std::shared_timed_mutex mutex;
      bool value;

    public:
      ThreadResult(): value(false) {}

      bool get() const
      {
        std::shared_lock<std::shared_timed_mutex> rlock(mutex);

        return value;
      }

      void set(const bool &value)
      {
        std::unique_lock<std::shared_timed_mutex> wlock(mutex);

        this->value = value;
      }
    };

    ThreadResult result;

    auto check_and_update = [&result, &set_obj](const BASIC_SET_TYPE &s) {
      if (!result.get() && set_obj.is_subset_of(s)) {
        result.set(true);
      }
    };

    ThreadPool::BatchId batch_id = thread_pool.create_batch();

    for (auto it = std::cbegin(*this); it != std::cend(*this); ++it) {
      // submit the task to the thread pool
      thread_pool.submit_to_batch(batch_id, check_and_update, std::ref(*it));
    }

    // join to the pool threads
    thread_pool.join_threads(batch_id);

    // close the batch
    thread_pool.close_batch(batch_id);

    return result.get();
#else  // WITH_THREADS
    for (auto it = std::cbegin(*this); it != std::cend(*this); ++it) {
      if (set_obj.is_subset_of(*it)) {
        return true;
      }
    }

    return false;
#endif // WITH_THREADS
  }

  /**
   * @brief Get the number of closed sets in the union
   *
   * @returns the number of closed sets in the union
   */
  size_t size() const
  {
    return std::list<BASIC_SET_TYPE>::size();
  }

  /**
   * @brief Get the space dimension of the sets
   *
   * @returns the space dimension of the sets or
   *          0 if the union is the empty
   */
  size_t dim() const
  {
    if (this->empty()) {
      return 0;
    }

    return this->front().dim();
  }

  /**
   * @brief Expand the union of closed sets
   *
   * This method expands the union so that each of its boundaries
   * is moved by a value `delta`.
   *
   * @param delta is the aimed expansion
   * @return a reference to the updated object
   */
  SetsUnion<BASIC_SET_TYPE> &expand_by(const double delta)
  {
    for (auto it = std::begin(*this); it != std::end(*this); ++it) {
      it->expand_by(delta);
    }

    return *this;
  }

  /**
   * @brief The begin iterator
   *
   * @return the begin iterator
   */
  iterator begin()
  {
    return std::list<BASIC_SET_TYPE>::begin();
  }

  /**
   * @brief The end iterator
   *
   * @return the end iterator
   */
  iterator end()
  {
    return std::list<BASIC_SET_TYPE>::end();
  }

  /**
   * @brief The constant begin iterator
   *
   * @return the constant begin iterator
   */
  const_iterator begin() const
  {
    return std::list<BASIC_SET_TYPE>::begin();
  }

  /**
   * @brief The constant end iterator
   *
   * @return the constant end iterator
   */
  const_iterator end() const
  {
    return std::list<BASIC_SET_TYPE>::end();
  }

  /**
   * @brief Test whether the union of closed sets is empty
   *
   * @return `true` if and only if the union of closed sets is empty
   */
  inline bool is_empty() const
  {
    // since sets are added to the union exclusively by using
    // add methods and those methods exclusively push in the
    // list non-empty set, the union is empty if and only if
    // the list is empty

    return this->empty();
  }

  /**
   * @brief Get the list of closed sets in the union
   *
   * @return the list of closed sets in the union
   */
  std::list<BASIC_SET_TYPE> get_list() const
  {
    std::list<BASIC_SET_TYPE> set_list;

    std::copy(begin(), end(), std::back_inserter(set_list));

    return set_list;
  }

  template<class BASIC_SET_TYPE2>
  friend class StickyUnion;
};

/**
 * @brief Test whether two unions of closed sets are disjoint
 *
 * @tparam BASIC_SET_TYPE is the basic set type
 * @param A is a sets union
 * @param B is a sets union
 * @return `true` if and only if `A` and `B` are disjoint
 */
template<class BASIC_SET_TYPE>
bool are_disjoint(const SetsUnion<BASIC_SET_TYPE> &A,
                  const SetsUnion<BASIC_SET_TYPE> &B)
{
  for (auto A_it = std::begin(A); A_it != std::end(A); ++A_it) {
    for (auto B_it = std::begin(B); B_it != std::end(B); ++B_it) {
      if (!are_disjoint(*A_it, *B_it)) {
        return false;
      }
    }
  }

  return true;
}

/**
 * @brief Test whether a sets union and a set are disjoint
 *
 * @tparam BASIC_SET_TYPE is the basic set type
 * @param A is a sets union
 * @param B is a set
 * @return `true` if and only if `A` and `B` are disjoint
 */
template<class BASIC_SET_TYPE>
bool are_disjoint(const SetsUnion<BASIC_SET_TYPE> &A, const BASIC_SET_TYPE &B)
{
  for (auto A_it = std::begin(A); A_it != std::end(A); ++A_it) {
    if (!are_disjoint(*A_it, B)) {
      return false;
    }
  }

  return true;
}

/**
 * @brief Test whether a set and a sets union are disjoint
 *
 * @tparam BASIC_SET_TYPE is the basic set type
 * @param A is a set
 * @param B is a sets union
 * @return `true` if and only if `A` and `B` are disjoint
 */
template<class BASIC_SET_TYPE>
inline bool are_disjoint(const BASIC_SET_TYPE &A,
                         const SetsUnion<BASIC_SET_TYPE> &B)
{
  return are_disjoint(B, A);
}

/**
 * Compute the intersection of two sets unions
 *
 * @tparam BASIC_SET_TYPE is the basic set type
 * @param[in] A is a sets union
 * @param[in] B is a sets union
 * @return the sets union representing \f$A \cap B\f$
 */
template<class BASIC_SET_TYPE>
SetsUnion<BASIC_SET_TYPE> intersect(const SetsUnion<BASIC_SET_TYPE> &A,
                                    const SetsUnion<BASIC_SET_TYPE> &B)
{
  SetsUnion<BASIC_SET_TYPE> result;

  for (auto A_it = std::begin(A); A_it != std::end(A); ++A_it) {
    for (auto B_it = std::begin(B); B_it != std::end(B); ++B_it) {
      result.add(intersect(*A_it, *B_it));
    }
  }

  return result;
}

/**
 * Compute the intersection between a sets unions and a set
 *
 * @tparam BASIC_SET_TYPE is the basic set type
 * @param[in] A is a sets union
 * @param[in] B is a set
 * @return the sets union representing \f$A \cap B\f$
 */
template<class BASIC_SET_TYPE>
inline SetsUnion<BASIC_SET_TYPE> intersect(const SetsUnion<BASIC_SET_TYPE> &A,
                                           const BASIC_SET_TYPE &B)
{
  return intersect(A, SetsUnion<BASIC_SET_TYPE>(B));
}

/**
 * Compute the intersection between a sets unions and a set
 *
 * @tparam BASIC_SET_TYPE is the basic set type
 * @param[in] A is a set
 * @param[in] B is a sets union
 * @return the bundles union representing \f$A \cap B\f$
 */
template<class BASIC_SET_TYPE>
inline SetsUnion<BASIC_SET_TYPE> intersect(const BASIC_SET_TYPE &A,
                                           const SetsUnion<BASIC_SET_TYPE> &B)
{
  return intersect(B, A);
}

/**
 * Compute the union of two sets unions
 *
 * @tparam BASIC_SET_TYPE is the basic set type
 * @param[in] A is a sets union
 * @param[in] B is a sets union
 * @return the sets union representing \f$A \cup B\f$
 */
template<class BASIC_SET_TYPE>
inline SetsUnion<BASIC_SET_TYPE> make_union(const SetsUnion<BASIC_SET_TYPE> &A,
                                            const SetsUnion<BASIC_SET_TYPE> &B)
{
  return SetsUnion<BASIC_SET_TYPE>(A).update(B);
}

/**
 * Compute the union of two sets
 *
 * @tparam BASIC_SET_TYPE is the basic set type
 * @param[in] A is a union of closed sets
 * @param[in] B is a set
 * @return the sets union representing \f$A \cup B\f$
 */
template<class BASIC_SET_TYPE>
inline SetsUnion<BASIC_SET_TYPE> make_union(const SetsUnion<BASIC_SET_TYPE> &A,
                                            const BASIC_SET_TYPE &B)
{
  return make_union(A, SetsUnion<BASIC_SET_TYPE>(B));
}

/**
 * Compute the union of two sets
 *
 * @tparam BASIC_SET_TYPE is the basic set type
 * @param[in] A is a set
 * @param[in] B is a set
 * @return the sets union representing \f$A \cup B\f$
 */
template<class BASIC_SET_TYPE>
inline SetsUnion<BASIC_SET_TYPE> make_union(const BASIC_SET_TYPE &A,
                                            const BASIC_SET_TYPE &B)
{
  SetsUnion<BASIC_SET_TYPE> result(A);

  result.add(B);

  return result;
}

/**
 * @brief Subtract a sets union to a set and close the result
 *
 * @param[in] set_obj is a closed set
 * @param[in] sets_union is a union of closed sets
 * @return a union of sets obtained by closing the difference
 *         between `set_obj` and `sets_union`
 */
template<class BASIC_SET_TYPE>
SetsUnion<BASIC_SET_TYPE>
subtract_and_close(const BASIC_SET_TYPE &set_obj,
                   const SetsUnion<BASIC_SET_TYPE> &sets_union)
{
  SetsUnion<BASIC_SET_TYPE> res = set_obj;

  for (auto u_it = std::begin(sets_union); u_it != std::end(sets_union);
       ++u_it) {
    SetsUnion<BASIC_SET_TYPE> tmp;
    for (auto s_it = std::begin(res); s_it != std::end(res); ++s_it) {
      tmp.update(subtract_and_close(*s_it, *u_it));
    }

    res = std::move(tmp);
  }

  return res;
}

/**
 * @brief Subtract a set to a sets union and close the result
 *
 * @param[in] sets_union is a union of closed sets
 * @param[in] set_obj is a closed set
 * @return a union of sets obtained by closing the difference
 *         between `sets_union` and `set_obj`
 */
template<class BASIC_SET_TYPE>
SetsUnion<BASIC_SET_TYPE>
subtract_and_close(const SetsUnion<BASIC_SET_TYPE> &sets_union,
                   const BASIC_SET_TYPE &set_obj)
{
  SetsUnion<BASIC_SET_TYPE> res;

  for (auto u_it = std::begin(sets_union); u_it != std::end(sets_union);
       ++u_it) {
    res.update(subtract_and_close(*u_it, set_obj));
  }

  return res;
}

/**
 * @brief Subtract a set to a sets union and close the result
 *
 * @param[in] minuend is a union of closed sets
 * @param[in] subtrahend is a union of closed sets
 * @return a union of sets obtained by closing the difference
 *         between `minuend` and `subtrahend`
 */
template<class BASIC_SET_TYPE>
SetsUnion<BASIC_SET_TYPE>
subtract_and_close(const SetsUnion<BASIC_SET_TYPE> &minuend,
                   const SetsUnion<BASIC_SET_TYPE> &subtrahend)
{
  SetsUnion<BASIC_SET_TYPE> res;

  for (auto s_it = std::begin(minuend); s_it != std::end(minuend); ++s_it) {
    res.update(subtract_and_close(*s_it, subtrahend));
  }

  return res;
}

/**
 * @brief Over-approximate the union of a union of sets
 *
 * @tparam BASIC_SET_TYPE is the type of the closed sets
 * @param sets_union is the union of basic sets
 * @return an over-approximation of the union of the
 *    sets in `sets_union`
 */
template<class BASIC_SET_TYPE>
BASIC_SET_TYPE
over_approximate_union(const SetsUnion<BASIC_SET_TYPE> &sets_union)
{
  BASIC_SET_TYPE union_set;
  for (auto it = std::begin(sets_union); it != std::end(sets_union); ++it) {
    if (union_set.dim() == 0) {
      union_set = *it;
    } else {
      union_set = over_approximate_union(union_set, *it);
    }
  }

  return union_set;
}

/**
 * @brief Over-approximate the union of a collection of sets
 *
 * @tparam BASIC_SET_TYPE is the type of the closed sets
 * @param container is the container of the closed sets
 * @return an over-approximation of the union of the
 *    sets in `container`
 */
template<class BASIC_SET_TYPE>
BASIC_SET_TYPE
over_approximate_union(const std::list<BASIC_SET_TYPE> &container)
{
  BASIC_SET_TYPE union_set;
  for (auto it = std::begin(container); it != std::end(container); ++it) {
    if (union_set.dim() == 0) {
      union_set = *it;
    } else {
      union_set = over_approximate_union(union_set, *it);
    }
  }

  return union_set;
}

/**
 * @brief Compute the chain-join of a set in a list of sets
 *
 * The chain-join of a basic set \f$S\f$ in a list \f$L\f$ of
 * sets is the smallest basic set \f$C\f$ that
 * over-approximates
 * \f$S \cup \bigcup_{i=1}^n L_i\f$, where \f$L_i \in L\f$ for
 * all $i \in [1,n]$, and does not intersect the remaining
 * sets in \f$L\f$.
 *
 * @tparam BASIC_SET_TYPE is the basic set type
 * @param sets_list is a list of basic sets
 * @param set_obj is a basic set
 * @return the chain-join of `set_obj` in `sets_list`
 */
template<class BASIC_SET_TYPE>
BASIC_SET_TYPE chain_join(std::list<BASIC_SET_TYPE> sets_list,
                          const BASIC_SET_TYPE &set_obj)
{
  BASIC_SET_TYPE joint(set_obj);
  bool fix_point_reached;

  do {
    std::list<BASIC_SET_TYPE> intersecting;
    fix_point_reached = true;

    for (auto it = std::begin(sets_list); it != std::end(sets_list);) {
      if (!are_disjoint(*it, set_obj)) {
        intersecting.push_back(std::move(*it));
        it = sets_list.erase(it);

        fix_point_reached = false;
      } else {
        ++it;
      }
    }

    if (!fix_point_reached) {
      intersecting.push_back(std::move(joint));
      joint = over_approximate_union(intersecting);
    }
  } while (fix_point_reached);

  return joint;
}

#endif /* SETSUNION_H_ */
