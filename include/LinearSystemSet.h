/**
 * @file LinearSystemSet.h
 * Represent and manipulate a set of linear systems
 * It can be used to represent a symbolic union of polytopes
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#ifndef LINEARSYSTEMSET_H_
#define LINEARSYSTEMSET_H_

#include <list>
#include <memory>

#include "LinearSystem.h"
#include "OutputFormatter.h"

#define MINIMIZE_LS_SET_REPRESENTATION true

class LinearSystemSet : private std::vector<LinearSystem>
{
public:
  using const_iterator = std::vector<LinearSystem>::const_iterator;
  using iterator = std::vector<LinearSystem>::iterator;

  /**
   * Constructor
   *
   */
  LinearSystemSet();

  /**
   * Constructor that instantiates a singleton set
   *
   * @param[in] ls element of the set
   */
  LinearSystemSet(LinearSystem &&ls);

  /**
   * Constructor that instantiates a singleton set
   *
   * @param[in] ls element of the set
   */
  LinearSystemSet(const LinearSystem &ls);

  /**
   * A copy constructor for a linear system set
   *
   * @param[in] orig a linear system set
   */
  LinearSystemSet(const LinearSystemSet &orig);

  /**
   * A swap constructor for a linear system set
   *
   * @param[in] orig is a rvalue linear system set
   */
  LinearSystemSet(LinearSystemSet &&orig);

  /**
   * A copy assignment for a linear system set
   *
   * @param[in] orig a linear system set
   */
  LinearSystemSet &operator=(const LinearSystemSet &orig);

  /**
   * A swap assignment for a linear system set
   *
   * @param[in] orig is a rvalue linear system set
   */
  LinearSystemSet &operator=(LinearSystemSet &&orig);

  /**
   * Add a linear system to the set
   *
   * @param[in] ls linear system to be added
   * @return a reference to the new linear system set
   */
  LinearSystemSet &add(const LinearSystem &ls);

  /**
   * Add a linear system to the set
   *
   * @param[in] ls linear system to be added
   * @return a reference to the new linear system set
   */
  LinearSystemSet &add(LinearSystem &&ls);

  LinearSystemSet &simplify();

  std::list<LinearSystem> get_a_finer_covering() const;

  // inplace operations on set
  /**
   * Update a linear system set by joining the linear systems in a set.
   *
   * This method works in-place and changes the calling object.
   *
   * @param[in] lss a linear system set.
   * @return a reference to the updated object.
   */
  LinearSystemSet &unionWith(const LinearSystemSet &LSset);

  /**
   * Update a linear system set by joining the linear systems in a set.
   *
   * This method works in-place and changes the calling object.
   *
   * @param[in] lss a linear system set.
   * @return a reference to the updated object.
   */
  LinearSystemSet &unionWith(LinearSystemSet &&LSset);

  /**
   * Update a linear system set by joining a linear system.
   *
   * This method works in-place and changes the calling object.
   *
   * @param[in] lss a linear system set.
   * @return a reference to the updated object.
   */
  LinearSystemSet &unionWith(const LinearSystem &ls)
  {
    this->add(ls);

    return *this;
  }

  /**
   * Update a linear system set by joining a linear system.
   *
   * This method works in-place and changes the calling object.
   *
   * @param[in] lss a linear system set.
   * @return a reference to the updated object.
   */
  LinearSystemSet &unionWith(LinearSystem &&ls)
  {
    this->add(ls);

    return *this;
  }

  double boundingVol() const;

  /**
   * Get the size of this set, i.e,
   * the number of linear systems
   *
   * @returns number of linear systems in the set
   */
  unsigned int size() const
  {
    return std::vector<LinearSystem>::size();
  }

  /**
   * Get the number of variables
   *
   * @returns number of columns of linear systems in the set
   */
  unsigned int dim() const;

  iterator begin()
  {
    return std::vector<LinearSystem>::begin();
  }

  iterator end()
  {
    return std::vector<LinearSystem>::end();
  }

  const_iterator begin() const
  {
    return std::vector<LinearSystem>::begin();
  }

  const_iterator end() const
  {
    return std::vector<LinearSystem>::end();
  }

  bool isEmpty() const;

  /**
   * Print the set of linear systems
   */
  void print() const;

  /**
   * Print the linear system set in Matlab format (for plotregion script)
   *
   * @param[in] os is the output stream
   * @param[in] color color of the polytope to plot
   */
  void plotRegion(std::ostream &os = std::cout, const char color = ' ') const;

  /**
   * Compute the intersection of two linear system sets
   *
   * @param[in] A is a linear system set
   * @param[in] B is a linear system set
   * @return the linear system set that represents the set of values satisfying
   * both the parameters
   */
  friend LinearSystemSet intersection(const LinearSystemSet &A,
                                      const LinearSystemSet &B);

  /**
   * Compute the intersection between a linear system sets and a linear system
   *
   * @param[in] A is a linear system set
   * @param[in] B is a linear system
   * @return the linear system set that represents the set of values satisfying
   * both the parameters
   */
  friend LinearSystemSet intersection(const LinearSystemSet &A,
                                      const LinearSystem &B);

  ~LinearSystemSet();
};

/**
 * Compute the intersection between a linear system sets and a linear system
 *
 * @param[in] A is a linear system
 * @param[in] B is a linear system set
 * @return the linear system set that represents the set of values satisfying
 * both the parameters
 */
inline LinearSystemSet intersection(const LinearSystem &A,
                                    const LinearSystemSet &B)
{
  return intersection(B, A);
}

/**
 * Compute the union of two linear system sets
 *
 * @param[in] A is a linear system set
 * @param[in] B is a linear system set
 * @return the linear system set that represents the set of values satisfying
 *          at least one of the parameters
 */
LinearSystemSet unionset(const LinearSystemSet &A, const LinearSystemSet &B);

/**
 * Compute the union of two linear systems
 *
 * @param[in] A is a linear system
 * @param[in] B is a linear system
 * @return the linear system set that represents the set of values satisfying
 *          at least one of the parameters.
 */
LinearSystemSet unionset(const LinearSystem &A, const LinearSystem &B);

/**
 * Test whether all the LinearSystemSet in a list are empty
 *
 * @param[in] lss_list a list of linear system set
 * @returns false if one of the linear system set in the list is not empty
 */
bool every_set_is_empty(const std::list<LinearSystemSet> &lss_list);

template<typename OSTREAM>
OSTREAM &operator<<(OSTREAM &os, const LinearSystemSet &ls)
{
  using OF = OutputFormatter<OSTREAM>;

  if (ls.size() == 0) {
    os << OF::empty_set();
  } else {
    os << OF::set_begin();
    for (auto it = ls.begin(); it != ls.end(); ++it) {
      if (it != ls.begin()) {
        os << OF::set_separator();
      }
      os << *it;
    }

    os << OF::set_end();
  }

  return os;
}

#endif /* LINEARSYSTEMSET_H_ */
