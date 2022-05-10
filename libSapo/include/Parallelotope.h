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
private:
  std::vector<LinearAlgebra::Vector<double>> _generators; //!< generators

  LinearAlgebra::Vector<double> _base_vertex; //!< base vertex
  std::vector<double> _lengths;               //!< tensors original lengths

public:
  /**
   * Constructor
   *
   * @param[in] directions is the vector of the parallelotope directions
   * @param[in] lower_bound is the lower bound offsets of the parallelotope
   * @param[in] upper_bound is the upper bound offsets of the parallelotope
   */
  Parallelotope(const std::vector<LinearAlgebra::Vector<double>> &directions, 
                const LinearAlgebra::Vector<double> &lower_bound,
                const LinearAlgebra::Vector<double> &upper_bound);

  /**
   * Get the parallelotope dimension
   *
   * @returns parallelotope dimension
   */
  unsigned int dim() const
  {
    return (_generators.size() == 0 ? 0 : _generators.front().size());
  }

  // TODO: test the following method
  /**
   * Build a polytope representing the parallelotope.
   *
   * @returns a polytope representing the parallelotope
   */
  operator Polytope() const;

  /**
   * @brief Get the base vertex
   * 
   * @return a reference to the base vertex
   */
  const LinearAlgebra::Vector<double> &base_vertex() const
  {
    return this->_base_vertex;
  }

  /**
   * @brief Get the generator lengths
   * 
   * @return a reference to the generator lengths
   */
  const std::vector<double> &lengths() const
  {
    return this->_lengths;
  }

  /**
   * @brief Get the parallelotope generators
   * 
   * @return a reference to the parallelotope generators
   */
  const std::vector<LinearAlgebra::Vector<double>> &generators() const
  {
    return this->_generators;
  }

  /**
   * @brief Swap two parallelotope
   * 
   * This method swap two parallelotope objects.
   * 
   * @param A is a parallelotope
   * @param B is a parallelotope
   */
  friend void swap(Parallelotope &A, Parallelotope &B);
};

#endif /* PARALLELOTOPE_H_ */
