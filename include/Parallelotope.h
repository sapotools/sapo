/**
 * @file Parallelotope.h
 * Represent and manipulate a parallelotope
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#ifndef PARALLELOTOPE_H_
#define PARALLELOTOPE_H_

#include <vector>

#include "Polytope.h"

/**
 * @brief A class to represent Parallelotope.
 *
 * Parallelotopes are internally represented as a base vertex,
 * a set of versors, and the lengths of the corresponding
 * tensors.
 */
class Parallelotope
{
  using Vector = std::vector<double>;
  using Matrix = std::vector<Vector>;

private:
  Matrix _versors; //!< versors

  std::vector<double> _base_vertex; //!< base vertex
  std::vector<double> _lengths;     //!< tensors original lengths

public:
  /**
   * Constructor
   *
   * @param[in] template_matrix is the template matrix of the parallelotope
   * @param[in] lower_bound is the lower bound offsets of the parallelotope
   * @param[in] upper_bound is the upper bound offsets of the parallelotope
   */
  Parallelotope(const Matrix &template_matrix, const Vector &lower_bound,
                const Vector &upper_bound);

  /**
   * Get the parallelotope dimension
   *
   * @returns parallelotope dimension
   */
  unsigned int dim() const
  {
    return (_versors.size() == 0 ? 0 : _versors.front().size());
  }

  // TODO: test the following method
  /**
   * Build a polytope representing the parallelotope.
   *
   * @returns a polytope representing the parallelotope
   */
  operator Polytope() const;

  const std::vector<double> &base_vertex() const
  {
    return this->_base_vertex;
  }

  const std::vector<double> &lengths() const
  {
    return this->_lengths;
  }

  // TODO: this method does not return versors. Why is it
  //       called in this way? Are we assuming that the
  //       returned values are versors?
  const Matrix &versors() const
  {
    return this->_versors;
  }

  friend void swap(Parallelotope &A, Parallelotope &B);
};

#endif /* PARALLELOTOPE_H_ */
