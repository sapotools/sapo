/**
 * @file Bundle.h
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @author Alberto Casagrande <acasagrande@units.it>
 * @brief Represent and manipulate bundles of parallelotopes
 * @version 0.2
 * @date 2022-05-04
 *
 * @copyright Copyright (c) 2016-2022
 */

#ifndef BUNDLE_H_
#define BUNDLE_H_

#include <set>

#include "LinearAlgebra.h"
#include "BaseConverter.h"
#include "Polytope.h"
#include "Parallelotope.h"

#include "STL/Atom.h"

#define SPLIT_MAGNITUDE_RATIO                                                 \
  0.75 //!< define the default versor magnitude multiplier for bundle splits

/**
 * @brief A class for parallelotope bundles
 *
 * A parallelotope bundle represents the intersection between
 * different non-singular parallelotopes.
 * This class stores all of the parallelotope directions/axes
 * in one single array and constraints each of them between
 * an upper and a lower bound.
 * The original parallelotopes can be rebuild by using the
 * a array of templates. Each template is a Natural-valued
 * array corresponding to a parallelotope. It stores which
 * directions are involved in the corresponding parallelotope.
 */
class Bundle
{
private:
  std::vector<LinearAlgebra::Vector<double>>
      _directions;                                //!< the vector of directions
  LinearAlgebra::Vector<double> _lower_bounds;    //!< direction upper bounds
  LinearAlgebra::Vector<double> _upper_bounds;    //!< direction lower bounds
  std::set<std::vector<unsigned int>> _templates; //!< templates matrix

  /**
   * Compute the edge lengths
   *
   * This method compute the edge lengths, i.e., the distances between the
   * half-spaces bounding the bundle.
   *
   * @returns the vector of the edge lengths
   */
  LinearAlgebra::Vector<double> edge_lengths();

  /**
   * @brief Add to the bundle templates for some directions
   *
   * @param missing_template_dirs the indices of the directions whose
   *                              templates we are missing
   */
  void add_templates_for(std::set<unsigned int> &missing_template_dirs);

public:
  /**
   * @brief A constructor
   *
   */
  Bundle();

  /**
   * @brief A copy constructor
   *
   * @param orig is the model for the new object
   */
  Bundle(const Bundle &orig);

  /**
   * @brief A move constructor
   *
   * @param orig is the model for the new object
   */
  Bundle(Bundle &&orig);

  /**
   * @brief A constructor
   *
   * @param[in] directions is the direction vector
   * @param[in] lower_bounds is the vector of direction lower bounds
   * @param[in] upper_bounds is the vector of direction upper bounds
   * @param[in] templates is the set of the templates
   */
  Bundle(const std::vector<LinearAlgebra::Vector<double>> &directions,
         const LinearAlgebra::Vector<double> &lower_bounds,
         const LinearAlgebra::Vector<double> &upper_bounds,
         const std::set<std::vector<unsigned int>> &templates);

  /**
   * @brief A move constructor
   *
   * @param[in] directions is the direction vector
   * @param[in] lower_bounds is the vector of direction lower bounds
   * @param[in] upper_bounds is the vector of direction upper bounds
   * @param[in] templates is the set of the templates
   */
  Bundle(std::vector<LinearAlgebra::Vector<double>> &&directions,
         LinearAlgebra::Vector<double> &&lower_bounds,
         LinearAlgebra::Vector<double> &&upper_bounds,
         const std::set<std::vector<unsigned int>> &templates);

  /**
   * @brief A constructor
   *
   * Whenever the templates are not specified at all, we assume that
   * all the directions are relevant and the templates are computed
   * automatically.
   *
   * @param[in] directions is the direction vector
   * @param[in] lower_bounds is the vector of direction lower bounds
   * @param[in] upper_bounds is the vector of direction upper bounds
   */
  Bundle(const std::vector<LinearAlgebra::Vector<double>> &directions,
         const LinearAlgebra::Vector<double> &lower_bounds,
         const LinearAlgebra::Vector<double> &upper_bounds);

  /**
   * @brief A move constructor
   *
   * Whenever the templates are not specified at all, we assume that
   * all the directions are relevant and the templates are computed
   * automatically.
   *
   * @param[in] directions is the direction vector
   * @param[in] lower_bounds is the vector of direction lower bounds
   * @param[in] upper_bounds is the vector of direction upper bounds
   */
  Bundle(std::vector<LinearAlgebra::Vector<double>> &&directions,
         LinearAlgebra::Vector<double> &&lower_bounds,
         LinearAlgebra::Vector<double> &&upper_bounds);

  /**
   * @brief Assignment operator
   *
   * @param orig is the original object to be copied
   * @return a reference to the updated object
   */
  Bundle &operator=(const Bundle &orig);

  /**
   * @brief Assignment operator
   *
   * @param orig is the original object to be copied
   * @return a reference to the updated object
   */
  Bundle &operator=(Bundle &&orig);

  /**
   * @brief Get the dimension of the bundle space
   *
   * @return the number of dimensions of the bundle space
   */
  unsigned int dim() const
  {
    return (_directions.size() == 0 ? 0 : _directions.front().size());
  }

  /**
   * @brief Get the number of templates in the bundle
   *
   * @return the number of templates in the bundle
   */
  unsigned int num_of_templates() const
  {
    return _templates.size();
  }

  /**
   * @brief Get the number of directions in the bundle
   *
   * @return the number of directions in the bundle
   */
  unsigned int size() const
  {
    return _directions.size();
  }

  /**
   * @brief Get the template vector
   *
   * @return a reference to the template vector
   */
  const std::set<std::vector<unsigned int>> &templates() const
  {
    return _templates;
  }

  /**
   * @brief Get the vector of directions
   *
   * @return a reference to the vector of directions
   */
  const std::vector<LinearAlgebra::Vector<double>> &directions() const
  {
    return _directions;
  }

  /**
   * @brief Get the i-th direction in the bundle
   *
   * @param i is the index of the aimed direction
   * @return  a reference to the i-th direction in the bundle
   */
  const LinearAlgebra::Vector<double> &get_direction(unsigned int i) const
  {
    return _directions[i];
  }

  /**
   * @brief Get the upper bound of the i-th direction
   *
   * @param i is the index of the direction whose
   *        upper bound must be returned
   * @return a reference to the upper bound of the i-th direction
   */
  const double &get_upper_bound(unsigned int i) const
  {
    return _upper_bounds[i];
  }

  /**
   * @brief Get the lower bound of the i-th direction
   *
   * @param i is the index of the direction whose
   *        lower bound must be returned
   * @return a reference to the lower bound of the i-th direction
   */
  const double &get_lower_bound(unsigned int i) const
  {
    return _lower_bounds[i];
  }

  /**
   * @brief Get the vector of the direction upper bounds
   *
   * @return a reference to the direction upper bounds
   */
  const std::vector<double> &upper_bounds() const
  {
    return _upper_bounds;
  }

  /**
   * @brief Get the vector of the direction lower bounds
   *
   * @return a reference to the direction lower bounds
   */
  const std::vector<double> &lower_bounds() const
  {
    return _lower_bounds;
  }

  /**
   * @brief Test whether the bundle is empty
   *
   * @return `true` if and only if the bundle is empty
   */
  bool is_empty() const;

  /**
   * @brief Test whether a bundle is subset of another bundle
   *
   * This method tests whether the current object is subset
   * of a bundle.
   *
   * @param[in] bundle is the tested bundle
   * @return `true` if and only if the current bundle is a
   *         subset of `bundle`
   */
  bool is_subset_of(const Bundle &bundle) const;

  /**
   * @brief Test whether a bundle is subset of a polytope
   *
   * This method tests whether the current object is subset
   * of a polytope.
   *
   * @param[in] P is a polytope
   * @return `true` if and only if the current bundle is a
   *         subset of `P`
   */
  inline bool is_subset_of(const Polytope &P) const
  {
    return this->satisfies(P);
  }

  /**
   * @brief Test whether a bundle includes another bundle
   *
   * This method tests whether a bundle is subset of the
   * current object.
   *
   * @param[in] bundle is the bundle whose inclusion is tested
   * @return `true` if and only if `bundle` is a subset of
   *         the current bundle
   */
  inline bool includes(const Bundle &bundle) const
  {
    return bundle.is_subset_of(*this);
  }

  /**
   * @brief Test whether a bundle includes a polytope
   *
   * @param P is the polytope whose inclusion is tested
   * @return `true` if and only if `P` is a subset of
   *         the current bundle
   */
  inline bool includes(const Polytope &P) const
  {
    return ((Polytope) * this).includes(P);
  }

  /**
   * @brief Check whether a bundle satisfies a linear system
   *
   * This method checks whether all the points in the
   * current object are solutions for a linear system.
   *
   * @param ls is the considered linear system
   * @return `true` if and only if all the points of
   *          the current bundle are solutions for `ls`
   */
  bool satisfies(const LinearSystem &ls) const;

  /**
   * @brief Generate the polytope represented by the bundle
   *
   * @returns polytope represented by the bundle
   */
  operator Polytope() const;

  /**
   * @brief Get a parallelotope of the bundle
   *
   * @param[in] bundle_template is the template for the bundle
   * @returns the parallelotope of the bundle associated to `bundle_template`
   */
  Parallelotope
  get_parallelotope(const std::vector<unsigned int> &bundle_template) const;

  /**
   * @brief Get a canonical bundle equivalent to the current one
   *
   * Get the canonical form for this bundle by minimizing the
   * difference between lower and upper bounds over all the
   * directions.
   *
   * @returns canonized bundle
   */
  Bundle get_canonical() const;

  /**
   * @brief Canonize the current bundle
   *
   * Turn this bundle in canonical form by minimizing the
   * difference between lower and upper bounds over all the
   * directions.
   *
   * @returns a reference to the canonized bundle
   */
  Bundle &canonize();

  /**
   * @brief Split the bundle in smaller sub-bundles
   *
   * This method splits the bundles whose maximal magnitude,
   * the maximal length of its generators, is greater than
   * `max_magnitude` into a list of sub-bundles whose
   * maximal magnitude is
   * \f$\textrm{max_magnitude}*\textrm{split_ratio}\f$.
   *
   * If \f$m\f$ in the maximal magnitude of the input bundle
   * and the bundle itself has \f$d\f$ directions,
   * this method may produces upto
   * \f$(\frac{m}{\textrm{max_magnitude}})^d\f$
   * sub-bundles.
   *
   * @param[in] max_magnitude is the maximal magnitude that
   *                          triggers the split
   * @param[in] split_ratio is the ratio between maximal
   *                        magnitude of the output bundles and
   *                        that that triggers the split
   * @return a list of sub-bundles whose union covers the
   *         input bundle
   */
  std::list<Bundle> split(const double max_magnitude,
                          const double split_ratio
                          = SPLIT_MAGNITUDE_RATIO) const;

  /**
   * Decompose the current symbolic polytope
   *
   * @param[in] alpha weight parameter in [0,1] for decomposition (0 for
   * distance, 1 for orthogonality)
   * @param[in] max_iters maximum number of randomly generated templates
   * @todo method does not work; it should be reviewed and fixed
   * @returns new bundle decomposing current symbolic polytope
   */
  Bundle decompose(double alpha, int max_iters);

  /**
   * @brief Get the intersection between two bundles
   *
   * This method intersects the current object and another bundle.
   * The result is stored in the current object and a reference
   * to it is returned.
   *
   * @param A is the intersecting bundle
   * @return a reference to the updated object
   */
  Bundle &intersect_with(const Bundle &A);

  /**
   * @brief Get the intersection between this bundle and a linear set
   *
   * This method intersects the current instance of the `Bundle` class
   * and a possible unbounded linear set provided as a `LinearSystem`.
   * The result is stored in the current object and a reference to it
   * is returned.
   *
   * @param ls is the intersecting linear system
   * @return a reference to the updated object
   */
  Bundle &intersect_with(const LinearSystem &ls);

  /**
   * @brief Expand the bundle
   *
   * This method expands the bundle so that each of its boundaries
   * is moved by a value `epsilon`.
   *
   * @param epsilon is the aimed expansion
   * @return a reference to the updated bundle
   */
  Bundle &expand_by(const double epsilon);

  virtual ~Bundle();

  friend void swap(Bundle &A, Bundle &B);

  friend Bundle over_approximate_union(const Bundle &b1, const Bundle &b2);
};

/**
 * @brief Swap the content of two bundles
 *
 * @param A is the first bundle to be swapped
 * @param B is the second bundle to be swapped
 */
void swap(Bundle &A, Bundle &B);

/**
 * @brief Equality between bundles
 *
 * This method tests whether the space represented by two bundles
 * is the same. Please, notice that this method does not consider
 * bundle templates; thus, event if two bundles are considered
 * equivalent by it, their evolutions according to the very same
 * dynamic laws may differ.
 *
 * @param b1 is a bundle
 * @param b2 is a bundle
 * @return `true` if and only if `b1` and `b2` represent the
 *         very same space set
 */
inline bool operator==(const Bundle &b1, const Bundle &b2)
{
  return b1.is_subset_of(b2) && b2.is_subset_of(b1);
}

/**
 * @brief Get the intersection between two bundles
 *
 * @param b1 is a bundle
 * @param b2 is a bundle
 * @return A bundle representing the intersection between
 *         the bundles `b1` and `b2`
 */
inline Bundle intersect(const Bundle &b1, const Bundle &b2)
{
  Bundle res(b1);

  res.intersect_with(b2);

  return res;
}

/**
 * @brief Get the intersection between a bundle and a linear system
 *
 * @param bundle is a bundle
 * @param linear_system is a linear system
 * @return A bundle representing the intersection between
 *         `bundle` and `linear_system`
 */
inline Bundle intersect(const Bundle &bundle,
                        const LinearSystem &linear_system)
{
  Bundle res(bundle);

  res.intersect_with(linear_system);

  return res;
}

/**
 * @brief Get the intersection between a linear system and a bundle
 *
 * @param linear_system is a linear system
 * @param bundle is a bundle
 * @return A bundle representing the intersection between
 *         `linear_system` and `bundle`
 */
inline Bundle intersect(const LinearSystem &linear_system,
                        const Bundle &bundle)
{
  return intersect(bundle, linear_system);
}

/**
 * @brief Compute the over-approximation of the union of two bundles
 *
 * @param b1 is a bundle
 * @param b2 is a bundle
 * @return A bundle over-approximating the union of `b1` and `b2`
 */
Bundle over_approximate_union(const Bundle &b1, const Bundle &b2);

#endif /* BUNDLE_H_ */
