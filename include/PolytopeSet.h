/**
 * @file PolytopeSet.h
 * Represent and manipulate a set of polytopes
 * It can be used to represent a symbolic union of polytopes
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#ifndef LINEARSYSTEMSET_H_
#define LINEARSYSTEMSET_H_

#include <list>
#include <memory>

#include "Polytope.h"
#include "OutputFormatter.h"

#define MINIMIZE_LS_SET_REPRESENTATION true

class PolytopeSet : private std::vector<Polytope>
{
public:
  using const_iterator = std::vector<Polytope>::const_iterator;
  using iterator = std::vector<Polytope>::iterator;

  /**
   * Constructor
   *
   */
  PolytopeSet();

  /**
   * Constructor that instantiates a singleton set
   *
   * @param[in] P is the polytope to be included in
   *      the set.
   */
  PolytopeSet(Polytope &&P);

  /**
   * Constructor that instantiates a singleton set
   *
   * @param[in] P is the polytope to be included in
   *      the set.
   */
  PolytopeSet(const Polytope &P);

  /**
   * A copy constructor for a polytope set
   *
   * @param[in] orig a polytope set
   */
  PolytopeSet(const PolytopeSet &orig);

  /**
   * A swap constructor for a polytope set
   *
   * @param[in] orig is a rvalue polytope set
   */
  PolytopeSet(PolytopeSet &&orig);

  /**
   * A copy assignment for a polytope set
   *
   * @param[in] orig a polytope set
   */
  PolytopeSet &operator=(const PolytopeSet &orig);

  /**
   * A swap assignment for a polytope set
   *
   * @param[in] orig is a rvalue polytope set
   */
  PolytopeSet &operator=(PolytopeSet &&orig);

  /**
   * Add a polytope to the set
   *
   * @param[in] P is the polytope to be added
   * @return a reference to the new polytope set
   */
  PolytopeSet &add(const Polytope &P);

  /**
   * Add a polytope to the set
   *
   * @param[in] ls polytope to be added
   * @return a reference to the new polytope set
   */
  PolytopeSet &add(Polytope &&ls);

  PolytopeSet &simplify();

  std::list<Polytope> get_a_finer_covering() const;

  // inplace operations on set
  /**
   * Update a polytope set by joining the polytopes of another set.
   *
   * This method works in-place and changes the calling object.
   *
   * @param[in] lss a polytope set.
   * @return a reference to the updated object.
   */
  PolytopeSet &unionWith(const PolytopeSet &LSset);

  /**
   * Update a polytope set by joining the polytopes of another set.
   *
   * This method works in-place and changes the calling object.
   *
   * @param[in] Ps is a polytope set.
   * @return a reference to the updated object.
   */
  PolytopeSet &unionWith(PolytopeSet &&Ps);

  /**
   * Update a polytope set by joining a polytope.
   *
   * This method works in-place and changes the calling object.
   *
   * @param[in] Ps is a polytope.
   * @return a reference to the updated object.
   */
  PolytopeSet &unionWith(const Polytope &P)
  {
    this->add(P);

    return *this;
  }

  /**
   * Update a polytope set by joining a polytope.
   *
   * This method works in-place and changes the calling object.
   *
   * @param[in] lss a polytope set.
   * @return a reference to the updated object.
   */
  PolytopeSet &unionWith(Polytope &&ls)
  {
    this->add(ls);

    return *this;
  }

  double volume_of_bounding_boxes() const;

  /**
   * Get the size of this set, i.e,
   * the number of polytopes
   *
   * @returns number of polytopes in the set
   */
  unsigned int size() const
  {
    return std::vector<Polytope>::size();
  }

  /**
   * Get the number of variables
   *
   * @returns number of columns of polytopes in the set
   */
  unsigned int dim() const;

  iterator begin()
  {
    return std::vector<Polytope>::begin();
  }

  iterator end()
  {
    return std::vector<Polytope>::end();
  }

  const_iterator begin() const
  {
    return std::vector<Polytope>::begin();
  }

  const_iterator end() const
  {
    return std::vector<Polytope>::end();
  }

  bool is_empty() const;

  /**
   * Print the set of polytopes
   */
  void print() const;

  /**
   * Print the polytope set in Matlab format (for plotregion script)
   *
   * @param[in] os is the output stream
   * @param[in] color color of the polytope to plot
   */
  void plotRegion(std::ostream &os = std::cout, const char color = ' ') const;

  /**
   * Compute the intersection of two polytope sets
   *
   * @param[in] A is a polytope set
   * @param[in] B is a polytope set
   * @return the polytope set that represents the set of values satisfying
   * both the parameters
   */
  friend PolytopeSet intersection(const PolytopeSet &A,
                                      const PolytopeSet &B);

  /**
   * Compute the intersection between a polytope sets and a polytope
   *
   * @param[in] A is a polytope set
   * @param[in] B is a polytope
   * @return the polytope set that represents the set of values satisfying
   * both the parameters
   */
  friend PolytopeSet intersection(const PolytopeSet &A,
                                      const Polytope &B);

  ~PolytopeSet();
};

/**
 * Compute the intersection between a polytope sets and a polytope
 *
 * @param[in] A is a polytope
 * @param[in] B is a polytope set
 * @return the polytope set that represents the set of values satisfying
 * both the parameters
 */
inline PolytopeSet intersection(const Polytope &A,
                                    const PolytopeSet &B)
{
  return intersection(B, A);
}

/**
 * Compute the union of two polytope sets
 *
 * @param[in] A is a polytope set
 * @param[in] B is a polytope set
 * @return the polytope set that represents the set of values satisfying
 *          at least one of the parameters
 */
PolytopeSet unionset(const PolytopeSet &A, const PolytopeSet &B);

/**
 * Compute the union of two polytopes
 *
 * @param[in] A is a polytope
 * @param[in] B is a polytope
 * @return the polytope set that represents the set of values satisfying
 *          at least one of the parameters.
 */
PolytopeSet unionset(const Polytope &A, const Polytope &B);

/**
 * Test whether all the PolytopeSet in a list are empty
 *
 * @param[in] ps_list a list of polytope set
 * @returns false if one of the polytope set in the list is not empty
 */
bool every_set_is_empty(const std::list<PolytopeSet> &ps_list);

template<typename OSTREAM>
OSTREAM &operator<<(OSTREAM &os, const PolytopeSet &Ps)
{
  using OF = OutputFormatter<OSTREAM>;

  if (Ps.size() == 0) {
    os << OF::empty_set();
  } else {
    os << OF::set_begin();
    for (auto it = std::cbegin(Ps); it != std::cend(Ps); ++it) {
      if (it != std::cbegin(Ps)) {
        os << OF::set_separator();
      }
      os << *it;
    }

    os << OF::set_end();
  }

  return os;
}

#endif /* LINEARSYSTEMSET_H_ */
