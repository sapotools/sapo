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

#include <memory>

#include "LinearSystem.h"

#include <list>

#define MINIMIZE_LS_SET_REPRESENTATION true

class LinearSystemSet
{

private:
  using pointer = std::shared_ptr<LinearSystem>;
  using container = std::vector<pointer>;
  container set; // set of linear systems

public:
  /**
   * An iterator class for linear system set
   */
  class iterator
  {

  private:
    container::iterator iit;

  public:
    using iterator_category = container::iterator::iterator_category;
    using difference_type = container::iterator::difference_type;
    using value_type = LinearSystem;
    using pointer = std::shared_ptr<LinearSystem>;
    using reference = LinearSystem &;

    /**
     *  Constructor
     */
    iterator(const container::iterator &it): iit(it) {}

    /**
     *  Copy constructor
     */
    iterator(const iterator &orig): iit(orig.iit) {}

    /**
     *  Swap constructor
     */
    iterator(iterator &&orig)
    {
      std::swap(iit, orig.iit);
    }

    /**
     *  Reference operator
     *
     * @return a reference to the pointed linear system.
     */
    reference operator*() const
    {
      return *(*iit);
    }

    /**
     *  Pointer operator
     *
     * @return a pointer to the pointed linear system.
     */
    pointer operator->()
    {
      return *iit;
    }

    /**
     *  Prefix increment
     *
     * @return this iterator after incrementing it.
     */
    iterator &operator++()
    {
      iit++;
      return *this;
    }

    /**
     *  Postfix increment
     *
     * @return this iterator before incrementing it.
     */
    iterator operator++(int)
    {
      return iterator(iit++);
    }

    /**
     * Increment
     *
     * @param[in] a is the increment.
     * @return this iterator after incrementing it.
     */
    iterator operator+(const int a)
    {
      return iterator(iit + a);
    }

    /**
     * Decrement
     *
     * @param[in] a is the decrement.
     * @return this iterator after decrementing it.
     */
    iterator operator-(const int a)
    {
      return iterator(iit - a);
    }

    /**
     *  Check whether two iterators point the same object.
     *
     * @param[in] a is an iterator
     * @param[in] b is an iterator
     * @return true if and only if the two iterators point the same object.
     */
    friend inline bool operator==(const iterator &a, const iterator &b)
    {
      return a.iit == b.iit;
    }

    /**
     *  Check whether two iterators point the same object.
     *
     * @param[in] a is an iterator
     * @param[in] b is an iterator
     * @return true if and only if the two iterators point different objects.
     */
    friend inline bool operator!=(const iterator &a, const iterator &b)
    {
      return a.iit != b.iit;
    }
  };

  /**
   * A constant iterator class for linear system set
   */
  class const_iterator
  {
  private:
    container::const_iterator iit;

  public:
    using iterator_category = container::const_iterator::iterator_category;
    using difference_type = container::const_iterator::difference_type;
    using value_type = LinearSystem;
    using pointer = std::shared_ptr<LinearSystem>;
    using reference = LinearSystem &;

    /**
     *  Constructor
     */
    const_iterator(const container::const_iterator &it): iit(it) {}

    /**
     *  Copy constructor
     */
    const_iterator(const const_iterator &orig): iit(orig.iit) {}

    /**
     *  Swap constructor
     */
    const_iterator(const_iterator &&orig)
    {
      std::swap(iit, orig.iit);
    }

    /**
     *  Reference operator
     *
     * @return a reference to the pointed linear system.
     */
    reference operator*() const
    {
      return *(*iit);
    }

    /**
     *  Pointer operator
     *
     * @return a pointer to the pointed linear system.
     */
    pointer operator->()
    {
      return *iit;
    }

    /**
     *  Prefix increment
     *
     * @return this iterator after incrementing it.
     */
    const_iterator &operator++()
    {
      iit++;
      return *this;
    }

    /**
     * Increment
     *
     * @param[in] a is the increment.
     * @return this iterator after incrementing it.
     */
    const_iterator operator+(const int a)
    {
      return const_iterator(iit + a);
    }

    /**
     * Decrement
     *
     * @param[in] a is the decrement.
     * @return this iterator after decrementing it.
     */
    const_iterator operator-(const int a)
    {
      return const_iterator(iit - a);
    }

    /**
     *  Postfix increment
     *
     * @return this iterator before incrementing it.
     */
    const_iterator operator++(int)
    {
      return const_iterator(iit++);
    }

    /**
     *  Check whether two iterators point the same object.
     *
     * @param[in] a a constant iterator.
     * @param[in] b a constant iterator.
     * @return true if and only if the two iterators point the same object.
     */
    friend inline bool operator==(const const_iterator &a,
                                  const const_iterator &b)
    {
      return a.iit == b.iit;
    }

    /**
     *  Check whether two constant iterators point the same object.
     *
     * @param[in] a is a constant iterator.
     * @param[in] b is a constant iterator.
     * @return true if and only if the two iterators point different objects.
     */
    friend inline bool operator!=(const const_iterator &a,
                                  const const_iterator &b)
    {
      return a.iit != b.iit;
    }
  };

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
   * Constructor that instantiates a singleton set
   *
   * @param[in] ls a pointer to the element of the set
   */
  LinearSystemSet(pointer ls);

  /**
   * Constructor that instantiates a set from a vector of sets
   *
   * @param[in] set vector of linear systems
   */
  LinearSystemSet(const container &set);

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
   * @param[in] ls is a pointer to the linear system to be added
   * @return a reference to the new linear system set
   */
  LinearSystemSet &add(pointer ls);

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
    return this->set.size();
  }

  /**
   * Get the number of variables
   *
   * @returns number of columns of linear systems in the set
   */
  unsigned int dim() const;

  iterator begin()
  {
    return iterator(std::begin(set));
  }

  iterator end()
  {
    return iterator(std::end(set));
  }

  const_iterator begin() const
  {
    return const_iterator(std::cbegin(set));
  }

  const_iterator end() const
  {
    return const_iterator(std::cend(set));
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
bool every_set_is_empty(const std::list<LinearSystemSet>& lss_list);

std::ostream &operator<<(std::ostream &out, const LinearSystemSet &ls);

#endif /* LINEARSYSTEMSET_H_ */
