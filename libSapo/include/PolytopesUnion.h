/**
 * @file PolytopesUnion.h
 * Represent and manipulate a union of polytopes.
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#ifndef POLYTOPESUNION_H_
#define POLYTOPESUNION_H_

#include <list>
#include <memory>

#include "Polytope.h"
#include "OutputFormatter.h"

/**
 * @brief Unions of polytopes
 *
 * This class represents unions of polytopes.
 */
class PolytopesUnion : private std::vector<Polytope>
{
public:
  /**
   * @brief A constant iterator alias
   */
  using const_iterator = std::vector<Polytope>::const_iterator;

  /**
   * @brief A iterator alias
   */
  using iterator = std::vector<Polytope>::iterator;

  /**
   * Constructor
   *
   */
  PolytopesUnion();

  /**
   * Constructor
   *
   * @param[in] P is the polytope to be included in
   *      the union.
   */
  PolytopesUnion(Polytope &&P);

  /**
   * Constructor
   *
   * @param[in] P is the polytope to be included in
   *      the union.
   */
  PolytopesUnion(const Polytope &P);

  /**
   * A copy constructor for a polytopes union
   *
   * @param[in] orig a polytopes union
   */
  PolytopesUnion(const PolytopesUnion &orig);

  /**
   * A swap constructor for a polytopes union
   *
   * @param[in] orig is a rvalue polytopes union
   */
  PolytopesUnion(PolytopesUnion &&orig);

  /**
   * A copy assignment for a polytopes union
   *
   * @param[in] orig a polytopes union
   */
  PolytopesUnion &operator=(const PolytopesUnion &orig);

  /**
   * A swap assignment for a polytopes union
   *
   * @param[in] orig is a rvalue polytopes union
   */
  PolytopesUnion &operator=(PolytopesUnion &&orig);

  /**
   * Add a polytope to the union
   *
   * @param[in] P is the polytope to be added
   * @return a reference to the new polytopes union
   */
  PolytopesUnion &add(const Polytope &P);

  /**
   * Add a polytope to the union
   *
   * @param[in] P is the polytope to be added
   * @return a reference to the new polytopes union
   */
  PolytopesUnion &add(Polytope &&P);

  /**
   * Simplify the representation of the polytopes union
   *
   * @return a reference to the simplified polytopes union
   */
  PolytopesUnion &simplify();

  /**
   * Split the polytopes union in a list of polytopes
   *
   * @return a list of polytopes
   */
  std::list<Polytope> split() const;

  // inplace operations on polytopes union
  /**
   * Update a polytopes union by joining another polytopes union.
   *
   * This method works in-place and changes the calling object.
   *
   * @param[in] Pu a polytopes union.
   * @return a reference to the updated object.
   */
  PolytopesUnion &add(const PolytopesUnion &Pu);

  /**
   * Update a polytopes union by joining another polytopes union.
   *
   * This method works in-place and changes the calling object.
   *
   * @param[in] Pu is a polytopes union.
   * @return a reference to the updated object.
   */
  PolytopesUnion &add(PolytopesUnion &&Pu);

  /**
   * Compute the sums of each polytope bounding box
   *
   * @return the sums of each polytope bounding box.
   */
  double volume_of_bounding_boxes() const;

  /**
   * Check whether one of the polytopes in a union includes a polytope
   *
   * @param[in] P is the polytopes whose inclusion must be tested
   * @return `true` if the parameter is included by one of the
   *         polytopes in the union; it returns `false`, otherwise
   */
  bool any_includes(const Polytope &P) const;

  /**
   * Get the number of polytopes in the union
   *
   * @returns the number of polytopes in the union
   */
  unsigned int size() const
  {
    return std::vector<Polytope>::size();
  }

  /**
   * Get the space dimension of the polytopes
   *
   * @returns the space dimension of the polytopes
   */
  unsigned int dim() const;

  /**
   * @brief The begin iterator
   *
   * @return the begin iterator
   */
  iterator begin()
  {
    return std::vector<Polytope>::begin();
  }

  /**
   * @brief The end iterator
   *
   * @return the end iterator
   */
  iterator end()
  {
    return std::vector<Polytope>::end();
  }

  /**
   * @brief The constant begin iterator
   *
   * @return the constant begin iterator
   */
  const_iterator begin() const
  {
    return std::vector<Polytope>::begin();
  }

  /**
   * @brief The constant end iterator
   *
   * @return the constant end iterator
   */
  const_iterator end() const
  {
    return std::vector<Polytope>::end();
  }

  /**
   * Test whether the union of polytopes is empty.
   *
   * @returns true, if the union of polytopes is empty; false, otherwise.
   */
  bool is_empty() const;

  /**
   * Compute the intersection of two polytopes unions
   *
   * @param[in] A is a polytopes union
   * @param[in] B is a polytopes union
   * @return the polytopes union representing the intersection of the
   *         parameters
   */
  friend PolytopesUnion intersect(const PolytopesUnion &A,
                                  const PolytopesUnion &B);

  /**
   * Compute the intersection between a polytopes unions and a polytope
   *
   * @param[in] A is a polytopes union
   * @param[in] B is a polytope
   * @return the polytopes union representing the intersection of the
   *         parameters
   */
  friend PolytopesUnion intersect(const PolytopesUnion &A, const Polytope &B);
};

/**
 * Compute the intersection between a polytopes unions and a polytope
 *
 * @param[in] A is a polytopes union
 * @param[in] B is a polytope
 * @return the polytopes union representing the intersection of the
 * parameters
 */
PolytopesUnion intersect(const PolytopesUnion &A, const Polytope &B);

/**
 * Compute the intersection between a polytopes unions and a polytope
 *
 * @param[in] A is a polytope
 * @param[in] B is a polytopes union
 * @return the polytopes union representing the intersection of the two
 * parameters
 */
inline PolytopesUnion intersect(const Polytope &A, const PolytopesUnion &B)
{
  return intersect(B, A);
}

/**
 * Compute the union of two polytopes unions
 *
 * @param[in] A is a polytopes union
 * @param[in] B is a polytopes union
 * @return the union of the two polytopes unions
 */
PolytopesUnion unite(const PolytopesUnion &A, const PolytopesUnion &B);

/**
 * Compute the union of two polytopes
 *
 * @param[in] A is a polytope
 * @param[in] B is a polytope
 * @return the union of the two polytopes
 */
PolytopesUnion unite(const Polytope &A, const Polytope &B);

/**
 * Test whether all the PolytopesUnion in a list are empty
 *
 * @param[in] ps_list a list of polytopes union
 * @returns false if one of the polytopes union in the list is not empty
 */
bool every_set_is_empty(const std::list<PolytopesUnion> &ps_list);

/**
 * @brief Print a polytope union in a stream
 *
 * @tparam OSTREAM is the output stream type
 * @param os is the output stream
 * @param Pu is the polytope union to be printed
 * @return a reference to the output stream
 */
template<typename OSTREAM>
OSTREAM &operator<<(OSTREAM &os, const PolytopesUnion &Pu)
{
  using OF = OutputFormatter<OSTREAM>;

  if (Pu.size() == 0) {
    os << OF::empty_set();
  } else {
    os << OF::set_begin();
    for (auto it = std::cbegin(Pu); it != std::cend(Pu); ++it) {
      if (it != std::cbegin(Pu)) {
        os << OF::set_separator();
      }
      os << *it;
    }

    os << OF::set_end();
  }

  return os;
}

#endif /* POLYTOPESUNION_H_ */
