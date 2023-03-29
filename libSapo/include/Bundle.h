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

#include "BundleTemplate.h"

#include "LinearAlgebra.h"
#include "Polytope.h"
#include "Parallelotope.h"
#include "Approximation.h"

#include "LinearAlgebraIO.h"

#define SPLIT_MAGNITUDE_RATIO                                                 \
  0.75 //!< define the default versor magnitude multiplier for bundle splits

/**
 * @brief Avoid \f$-0\f$
 *
 * This macro avoids \f$-0\f$ by replacing it by \f$0\f$.
 */
#define AVOID_NEG_ZERO(value) ((value) == 0 ? 0 : (value))

template<typename T, typename APPROX_TYPE>
class Bundle;

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
 * @brief Subtract two bundle and close the result
 *
 * @tparam T is the numeric type used to represent bundle
 * @tparam APPROX_TYPE is the type of approximation used in place of T
 *      during the computation
 * @param b1 is a bundle
 * @param b2 is a bundle
 * @return a union of bundles obtained by closing the set
 *         \f$b1\setminus b2\f$
 */
template<typename T, typename APPROX_TYPE>
SetsUnion<Bundle<T, APPROX_TYPE>>
subtract_and_close(const Bundle<T, APPROX_TYPE> &b1,
                   const Bundle<T, APPROX_TYPE> &b2);

/**
 * @brief A class for parallelotope bundles
 *
 * A parallelotope bundle represents the intersection between different
 * non-singular parallelotopes. This class stores all of the parallelotope
 * directions/axes in one single array and constraints each of them between an
 * upper and a lower bound. The parallelotopes are represented by means of
 * bundle templates. A template is the vector of the indices of the bundle
 * directions involved in the corresponding parallelotope.
 *
 * @tparam T is the numeric type used to represent bundle
 * @tparam APPROX_TYPE is the type of approximation used in place of T
 *      during the computation
 */
template<typename T = double, typename APPROX_TYPE = T>
class Bundle
{
  std::vector<LinearAlgebra::Vector<T>>
      _directions; //!< the vector of directions
  std::set<size_t>
      _adaptive_directions; //!< the set of adaptive direction indices
  LinearAlgebra::Vector<T> _lower_bounds; //!< direction upper bounds
  LinearAlgebra::Vector<T> _upper_bounds; //!< direction lower bounds
  std::set<BundleTemplate> _templates;    //!< bundle templates

  /**
   * @brief A draft horse function to split a bundle
   *
   * This recursive function is a private draft horse to
   * split a bundle whose maximal magnitude, the maximal
   * length of its generators, is greater than
   * `max_magnitude`. The bundle is split into a
   * list of sub-bundles whose maximal magnitude is
   * \f$\textrm{max_magnitude}*\textrm{split_ratio}\f$.
   * Each sub-bundle in the output list is built recursively
   * by considering one bundle direction per time and:
   * 1. splitting the direction in the opportune number of
   *    chunks, so to satisfy the sub-bundle maximal
   *    magnitude request;
   * 2. for each of the chunks, setting the sub-bundle
   *    boundaries of the selected direction according
   *    with the considered chunk and recursively
   *    consider the remaining directions.
   *
   * If \f$m\f$ in the maximal magnitude of the input bundle
   * and the bundle itself has \f$d\f$ directions,
   * this method may produces upto
   * \f$(\frac{m}{\textrm{max_magnitude}})^d\f$
   * sub-bundles.
   *
   * @param[out] res is the list of the resulting sub-bundles
   * @param[in, out] lower_bounds is the lower boundary vector
   *                              of the in-building sub-bundles
   * @param[in, out] upper_bounds is the upper boundary vector
   *                              of the in-building sub-bundles
   * @param[in] idx is the index of the direction to be split
   * @param[in] max_magnitude is the maximal magnitude that
   *                          triggers the split
   * @param[in] split_ratio is the ratio between maximal
   *                        magnitude of the output bundles and
   *                        that that triggers the split
   * @return a reference to the list of sub-bundles produced
   *         splitting the input bundle
   */
  std::list<Bundle<T, APPROX_TYPE>> &
  split_bundle(std::list<Bundle<T, APPROX_TYPE>> &res,
               LinearAlgebra::Vector<T> &lower_bounds,
               LinearAlgebra::Vector<T> &upper_bounds, const size_t idx,
               const T &max_magnitude, const T &split_ratio) const;

  /**
   * @private
   * @brief Validate directions and bounds
   *
   * This method validates directions and bounds of the current
   * bundle. It checks that all the directions have the same
   * dimension, that the number of directions corresponds to both
   * upper and lower bounds dimensions, and that the adaptive
   * direction set contains valid directions for the current bundle.
   * Whenever these conditions are not satisfied, the method thows a
   * `std::domain_error` exception.
   *
   * @throw `std::domain_error` if at least one of the following
   *        conditions holds: 1) `_directions` is empty, 2) two
   *        directions in `_directions` have a different size, 3)
   *        the size of `_directions` does not match the that of
   *        either `_upper_bounds` or `_lower_bounds`, or 4)
   *        `_adaptive_directions` contains a value that is not
   *        a valid index for `_directions`
   */
  void validate_directions_and_bounds() const;

  /**
   * @private
   * @brief Test whether a template is a complete basis for the space
   *
   * @param bundle_template is the template to be tested
   * @return true if and only if the set of the directions associated to
   *       bundle_template is a complete basis for the space
   */
  bool
  is_a_complete_basis(const std::vector<unsigned int> &bundle_template) const;

  /**
   * @brief Filter duplicate directions
   *
   * This method filter out duplicate directions from the current bundle.
   * Both lower and upper bound vectors and the adaptive direction set
   * are also updated to maintain their consistency with respect to the
   * direction vector. The method returns a vector that mapping the
   * old direction indices in the indices of the new direction vector.
   *
   * @return a vector mapping each index of the input direction vector
   *       in the corresponding index in the final direction vector
   */
  std::vector<size_t> filter_duplicated_directions();

  /**
   * @brief Validate a set of templates for the bundle
   *
   * This method validates a set of template for the current
   * bundle. It checks whether each of the templates corresponds
   * to a set linear independent directions whose cardinality
   * is the space dimension. If this is not the case, the
   * method throw a `std::domain_error` exception.
   *
   * @param templates is a set of possible template for the
   *        current bundle
   * @throw `std::domain_error` if not all the templates are
   *        valid templates for the current bundle
   */
  void validate_templates(
      const std::set<std::vector<unsigned int>> &templates) const;

  /**
   * @brief Duplicate multiple occurrencies of adaptive directions
   *
   * Adaptive directions change during the computation by definition.
   * Their changes depend on the template they belong to. Thus, even
   * though the same adaptive direction may appears in different
   * templates in the specification phase, that direction may change
   * in differently from template to template during the computation.
   * Because of this reason, every adaptive direction must appears
   * exacly once among all the templates. Whenever it was mentioned
   * more than once during the specification phase, every of its
   * occurences must be distinct from the others by duplicating the
   * original direction.
   * This method takes care of duplicating multiple occurences of
   * the same direction and updates `_directions`,
   * `_adaptive_directions`, `_lower_bounds` and `_upper_bounds`.
   * The method also returns an updated version of the template set
   * passed as parameter.
   *
   * @param templates is the set of current bundle templates
   *                  represented as sorted vectors
   * @return an version of `templates` updated according to the
   *         direction duplication
   */
  std::set<std::vector<unsigned int>> duplicate_adaptive_directions(
      std::set<std::vector<unsigned int>> &templates);

  /**
   * @brief Reorganize directions according to a set of new positions
   *
   * This function reorganizes the directions according to a new set
   * of positions and removes all the directions that do not have a
   * new position.
   *
   * @param new_positions is the map of the new positions
   */
  void
  resort_directions(const std::map<unsigned int, unsigned int> &new_positions);

  /**
   * @private
   * @brief Find new positions for all the positions in a set
   *
   * This method finds a new position for a set of mentioned
   * directions and gets rid of the non-metioned positions.
   *
   * @param positions is the set of position which must be replaced
   * @return a map assigning a new position to each of the old
   *        positions
   */
  static std::map<unsigned int, unsigned int>
  find_new_position_for(const std::set<unsigned int> &positions);

  /**
   * @private
   * @brief Search for linearly dependent bundle directions
   *
   * This function search for a bundle directions linearly
   * dependent to a vector. If such a direction exists, then
   * this method returns its index in the direction vector.
   * Otherwise, the method returns the first invalid index
   * of the direction vector, i.e., its size.
   *
   * @param v is a linear vector
   * @return if there exists a bundle directions linearly
   *         dependent to a vector `v`, then this method
   *         returns its index in the direction vector.
   *         Otherwise, the size of the direction vector
   */
  size_t
  get_a_linearly_dependent_direction(const LinearAlgebra::Vector<T> &v) const;

  /**
   * @brief Add missing templates
   *
   * This method adds some templates mapping old direction
   * indices to new indices through a provided map.
   *
   * @param new_templates is new template set
   * @param new_direction_indices is the map for the new
   *                              direction indices
   */
  void add_templates(const std::set<BundleTemplate> &new_templates,
                     const std::vector<unsigned int> &new_direction_indices);

  /**
   * @private
   * @brief Mark the directions linearly dependent to the i-th direction
   *
   * This method identifies all the bundle directions that are linearly
   * dependent to the `i`-th directions and whose index in the direction
   * vector is greater than `i`. This is done by assing the `j`-th
   * position of the `indices` vector to `i` every time `j>i` and
   * `_directions[j]` and `_directions[i]` are linearly dependent.
   *
   * @param[in, out] indices is a vector of valid indices for the bundle
   *        direction vector such that `indices[j]<=j` and, for any
   *        \f$j<i\f$, `indices[j]==k!=j` if and only if `k` is the
   *        smallest index such that and `_directions[j]` and
   *        `_direction[k]` are linearly dependent
   * @param[in] i is a valid index for the direction vector
   */
  void mark_linearly_dependent_directions(std::vector<size_t> &indices,
                                          const size_t &i) const;

  /**
   * @brief A constructor
   *
   * This constructor does not perform any check on the parameter
   * consistency.
   *
   * @param[in] adaptive_directions is the set of adaptive direction indices
   * @param[in] directions is the direction vector
   * @param[in] lower_bounds is the vector of direction lower bounds
   * @param[in] upper_bounds is the vector of direction upper bounds
   * @param[in] templates is the set of the bundle templates
   */
  Bundle(const std::set<size_t> &adaptive_directions,
         const std::vector<LinearAlgebra::Vector<T>> &directions,
         const LinearAlgebra::Vector<T> &lower_bounds,
         const LinearAlgebra::Vector<T> &upper_bounds,
         const std::set<BundleTemplate> &templates);

  /**
   * @brief A constructor
   *
   * This constructor does not perform any check on the parameter
   * consistency.
   *
   * @param[in] adaptive_directions is the set of adaptive direction indices
   * @param[in] directions is the direction vector
   * @param[in] lower_bounds is the vector of direction lower bounds
   * @param[in] upper_bounds is the vector of direction upper bounds
   * @param[in] templates is the set of the bundle templates
   */
  Bundle(const std::set<size_t> &adaptive_directions,
         std::vector<LinearAlgebra::Vector<T>> &&directions,
         LinearAlgebra::Vector<T> &&lower_bounds,
         LinearAlgebra::Vector<T> &&upper_bounds,
         const std::set<BundleTemplate> &templates);

public:
  /**
   * @brief A constructor
   *
   */
  Bundle();

  /**
   * @brief A copy constructor
   *
   * @tparam APPROX2 is the approximation type of the
   *         original polytope
   * @param orig is the model for the new object
   */
  template<typename APPROX2 = APPROX_TYPE>
  Bundle(const Bundle<T, APPROX2> &orig);

  /**
   * @brief A move constructor
   *
   * @tparam APPROX2 is the approximation type of the
   *         original polytope
   * @param orig is the model for the new object
   */
  template<typename APPROX2 = APPROX_TYPE>
  Bundle(Bundle<T, APPROX2> &&orig);

  /**
   * @brief A constructor
   *
   * @param[in] directions is the direction vector
   * @param[in] lower_bounds is the vector of direction lower bounds
   * @param[in] upper_bounds is the vector of direction upper bounds
   * @param[in] templates is the set of the bundle templates
   * @param[in] adaptive_direction is the set of adaptive directions
   * @param[in] remove_unused_directions is a flag to remove
   *      directions not belonging to any template
   */
  Bundle(const std::vector<LinearAlgebra::Vector<T>> &directions,
         const LinearAlgebra::Vector<T> &lower_bounds,
         const LinearAlgebra::Vector<T> &upper_bounds,
         std::set<std::vector<unsigned int>> templates,
         const std::set<size_t> &adaptive_directions,
         const bool remove_unused_directions = false);

  /**
   * @brief A constructor
   *
   * @param[in] directions is the direction vector
   * @param[in] lower_bounds is the vector of direction lower bounds
   * @param[in] upper_bounds is the vector of direction upper bounds
   * @param[in] templates is the set of the bundle templates
   * @param[in] remove_unused_directions is a flag to remove
   *      directions not belonging to any template
   */
  Bundle(const std::vector<LinearAlgebra::Vector<T>> &directions,
         const LinearAlgebra::Vector<T> &lower_bounds,
         const LinearAlgebra::Vector<T> &upper_bounds,
         const std::set<std::vector<unsigned int>> &templates,
         const bool remove_unused_directions = false);

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
   * @param[in] remove_unused_directions is a flag to remove
   *      directions not belonging to any template
   */
  Bundle(const std::vector<LinearAlgebra::Vector<T>> &directions,
         const LinearAlgebra::Vector<T> &lower_bounds,
         const LinearAlgebra::Vector<T> &upper_bounds,
         const bool remove_unused_directions = false);

  /**
   * @brief A constructor
   *
   * Build the bundle representing a polytope
   *
   * @tparam APPROX2 is the approximation type of the
   *         original polytope
   * @param P is the polytope whose bundle representation is aimed
   */
  template<typename APPROX2 = APPROX_TYPE>
  explicit Bundle(const Polytope<T, APPROX2> &P);

  /**
   * @brief Assignment operator
   *
   * @tparam APPROX2 is the approximation type of the
   *         original polytope
   * @param orig is the original object to be copied
   * @return a reference to the updated object
   */
  template<typename APPROX2 = APPROX_TYPE>
  Bundle<T, APPROX_TYPE> &operator=(const Bundle<T, APPROX2> &orig);

  /**
   * @brief Assignment operator
   *
   * @tparam APPROX2 is the approximation type of the
   *         original polytope
   * @param orig is the original object to be copied
   * @return a reference to the updated object
   */
  template<typename APPROX2 = APPROX_TYPE>
  Bundle<T, APPROX_TYPE> &operator=(Bundle<T, APPROX2> &&orig);

  /**
   * @brief Get the dimension of the bundle space
   *
   * @return the number of dimensions of the bundle space
   */
  inline size_t dim() const
  {
    return (_directions.size() == 0 ? 0 : _directions.front().size());
  }

  /**
   * @brief Get the number of templates in the bundle
   *
   * @return the number of templates in the bundle
   */
  inline size_t num_of_templates() const
  {
    return _templates.size();
  }

  /**
   * @brief Get the number of directions in the bundle
   *
   * @return the number of directions in the bundle
   */
  inline size_t size() const
  {
    return _directions.size();
  }

  /**
   * @brief Get the template vector
   *
   * @return a reference to the template vector
   */
  inline const std::set<BundleTemplate> &templates() const
  {
    return _templates;
  }

  /**
   * @brief Get the vector of directions
   *
   * @return a reference to the vector of directions
   */
  inline const std::vector<LinearAlgebra::Vector<T>> &directions() const
  {
    return _directions;
  }

  /**
   * @brief Get the i-th direction in the bundle
   *
   * @param i is the index of the aimed direction
   * @return  a reference to the i-th direction in the bundle
   */
  inline const LinearAlgebra::Vector<T> &get_direction(const size_t &i) const
  {
    return _directions[i];
  }

  /**
   * @brief Check whether the i-th direction is adaptive
   *
   * @param i is the index of the considered direction
   * @return `true` if and only if the i-th direction
   *      is adaptive
   */
  inline bool is_direction_adaptive(const size_t &i) const
  {
    return _adaptive_directions.count(i) > 0;
  }

  /**
   * @brief Get the set of the bundle adaptive direction indices
   *
   * @return the set of bundle adaptive direction indices
   */
  inline const std::set<size_t> &adaptive_directions() const
  {
    return _adaptive_directions;
  }

  /**
   * @brief Get the upper bound of the i-th direction
   *
   * @param i is the index of the direction whose
   *        upper bound must be returned
   * @return a reference to the upper bound of the i-th direction
   */
  const double &get_upper_bound(const size_t &i) const
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
  const double &get_lower_bound(const size_t &i) const
  {
    return _lower_bounds[i];
  }

  /**
   * @brief Get the vector of the direction upper bounds
   *
   * @return a reference to the direction upper bounds
   */
  const std::vector<T> &upper_bounds() const
  {
    return _upper_bounds;
  }

  /**
   * @brief Get the vector of the direction lower bounds
   *
   * @return a reference to the direction lower bounds
   */
  const std::vector<T> &lower_bounds() const
  {
    return _lower_bounds;
  }

  /**
   * @brief Test whether the bundle is empty
   *
   * @return `true` when it can establish that the bundle is
   *         empty. `false` when it can establish that the bundle
   *         is not empty. `uncertain` in the other cases
   */
  inline TriBool is_empty() const
  {
    return static_cast<Polytope<T, APPROX_TYPE>>(*this).is_empty();
  }

  /**
   * @brief Test whether the bundle interior is empty
   *
   * @return `true` when it can establish that the bundle interior
   *         is empty. `false` when it can establish that the bundle
   *         interior is not empty. `uncertain` in the other cases
   */
  inline TriBool is_interior_empty() const
  {
    return static_cast<Polytope<T, APPROX_TYPE>>(*this).is_interior_empty();
  }

  /**
   * @brief Test whether a bundle is subset of another bundle
   *
   * This method tests whether the current object is subset
   * of a bundle.
   *
   * @param[in] bundle is the tested bundle
   * @return `true` if the method can establish that this bundle
   *         is a subset of `P`. `false` when it can establish that
   *         this bundle is not a subset of `P`. `uncertain` in
   *         the remaining cases
   */
  TriBool is_subset_of(const Bundle<T, APPROX_TYPE> &bundle) const;

  /**
   * @brief Test whether a bundle is subset of a polytope
   *
   * This method tests whether the current object is subset
   * of a polytope.
   *
   * @param[in] P is a polytope
   * @return `true` if the method can establish that this bundle
   *         is a subset of `P`. `false` when it can establish that
   *         this bundle is not a subset of `P`. `uncertain` in
   *         the remaining cases
   */
  inline TriBool is_subset_of(const Polytope<T, APPROX_TYPE> &P) const
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
   * @return `true` if the method can establish that this bundle
   *         is a superset of `P`. `false` when it can establish that
   *         this bundle is not a superset of `P`. `uncertain` in
   *         the remaining cases
   */
  inline TriBool includes(const Bundle<T, APPROX_TYPE> &bundle) const
  {
    return bundle.is_subset_of(*this);
  }

  /**
   * @brief Test whether a bundle includes a polytope
   *
   * @param P is the polytope whose inclusion is tested
   * @return `true` if the method can establish that this bundle
   *         is a superset of `P`. `false` when it can establish that
   *         this bundle is not a superset of `P`. `uncertain` in
   *         the remaining cases
   */
  inline TriBool includes(const Polytope<T, APPROX_TYPE> &P) const
  {
    return static_cast<Polytope<T, APPROX_TYPE>>(*this).includes(P);
  }

  /**
   * @brief Check whether a bundle satisfies a linear system
   *
   * This method checks whether all the points in the
   * current object are solutions for a linear system.
   *
   * @param ls is the considered linear system
   * @return `true` if the method can establish that all the
   *          points in the current bundle are solutions for
   *          `ls`. `false` if it can establish that some of
   *          of the points in the current bundle are not
   *          solutions for `ls`. `uncertain` in the
   *          remaining cases
   */
  TriBool satisfies(const LinearSystem<T, APPROX_TYPE> &ls) const;

  /**
   * @brief Generate the polytope represented by the bundle
   *
   * @returns polytope represented by the bundle
   */
  operator Polytope<T, APPROX_TYPE>() const;

  /**
   * @brief Get a parallelotope of the bundle
   *
   * @param[in] bundle_template is the bundle template
   * @returns the parallelotope of the bundle associated to `bundle_template`
   */
  Parallelotope<T, APPROX_TYPE>
  get_parallelotope(const BundleTemplate &bundle_template) const;

  /**
   * @brief Get a canonical bundle equivalent to the current one
   *
   * Get the canonical form for this bundle by minimizing the
   * difference between lower and upper bounds over all the
   * directions.
   *
   * @returns canonized bundle
   */
  Bundle<T, APPROX_TYPE> get_canonical() const;

  /**
   * @brief Canonize the current bundle
   *
   * Turn this bundle in canonical form by minimizing the
   * difference between lower and upper bounds over all the
   * directions.
   *
   * @returns a reference to the canonized bundle
   */
  Bundle<T, APPROX_TYPE> &canonize();

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
  std::list<Bundle<T, APPROX_TYPE>> split(const T max_magnitude,
                                          const T split_ratio
                                          = SPLIT_MAGNITUDE_RATIO) const;

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
  Bundle<T, APPROX_TYPE> &intersect_with(const Bundle<T, APPROX_TYPE> &A);

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
  Bundle<T, APPROX_TYPE> &
  intersect_with(const LinearSystem<T, APPROX_TYPE> &ls);

  /**
   * @brief Expand the bundle
   *
   * This method expands the bundle so that each of its boundaries
   * is moved by a value `delta`.
   *
   * @param delta is the aimed expansion
   * @return a reference to the updated bundle
   */
  Bundle<T, APPROX_TYPE> &expand_by(const T delta);

  /**
   * Compute the edge lengths
   *
   * This method compute the edge lengths, i.e., the distances between the
   * half-spaces bounding the bundle.
   *
   * @returns the vector of the edge lengths
   */
  LinearAlgebra::Vector<T> edge_lengths() const;

  virtual ~Bundle();

  template<typename T2, typename APPROX1, typename APPROX2>
  friend void swap(Bundle<T2, APPROX1> &A, Bundle<T2, APPROX2> &B);

  template<typename T2, typename APPROX2>
  friend Bundle<T2, APPROX2>
  over_approximate_union(const Bundle<T2, APPROX2> &b1,
                         const Bundle<T2, APPROX2> &b2);

  template<typename T2, typename APPROX2>
  friend SetsUnion<Bundle<T2, APPROX2>>
  subtract_and_close(const Bundle<T2, APPROX2> &b1,
                     const Bundle<T2, APPROX2> &b2);

  template<typename T2>
  friend class Evolver;

  template<typename T2>
  friend class CachedEvolver;
};

/**
 * @brief Test whether two bundles are disjoint
 *
 * @tparam T is the numeric type used to represent bundle
 * @tparam APPROX_TYPE is the type of approximation used in place of T
 *      during the computation
 * @param A is a bundle
 * @param B is a bundle
 * @return `true` if and only if `A` and `B` are disjoint
 * @return `true` when the method can establish that `A` and
 *         `B` are disjoint. `false` when it can establish
 *         that `b1` and `b2` are not disjoint. `uncertain`
 *         in the remaining cases
 */

template<typename T, typename APPROX_TYPE>
TriBool are_disjoint(const Bundle<T, APPROX_TYPE> &A,
                     const Bundle<T, APPROX_TYPE> &B);

/**
 * @brief Swap the content of two bundles
 *
 * @tparam T is the numeric type used to represent bundle
 * @tparam APPROX1 is the type of approximation used in place of T
 *      during the computation on the first bundle
 * @tparam APPROX2 is the type of approximation used in place of T
 *      during the computation on the second bundle
 * @param A is the first bundle to be swapped
 * @param B is the second bundle to be swapped
 */
template<typename T, typename APPROX1, typename APPROX2>
void swap(Bundle<T, APPROX1> &A, Bundle<T, APPROX2> &B);

/**
 * @brief Equality between bundles
 *
 * @tparam T is the numeric type used to represent bundle
 * @tparam APPROX_TYPE is the type of approximation used in place of T
 *      during the computation
 * This method tests whether the space represented by two bundles
 * is the same. Please, notice that this method does not consider
 * bundle templates; thus, event if two bundles are considered
 * equivalent by it, their evolutions according to the very same
 * dynamic laws may differ.
 *
 * @param b1 is a bundle
 * @param b2 is a bundle
 * @return `true` when the method can establish that `b1` and
 *         `b2` are the same. `false` when it can  establish
 *         that `b1` and `b2` differ. `uncertain` in the
 *         remaining cases
 */
template<typename T, typename APPROX_TYPE>
TriBool operator==(const Bundle<T, APPROX_TYPE> &b1,
                   const Bundle<T, APPROX_TYPE> &b2);

/**
 * @brief Test whether two sets are the same
 *
 * @tparam T is the numeric type used to represent bundle
 * @tparam APPROX_TYPE is the type of approximation used in place of T
 *      during the computation
 * @param A is a polytope
 * @param B is a bundle
 * @return `true` when the method can establish that `A` and `B`
 *         represent the same set. `false` when it can
 *         establish `A` and `B` differ. `uncertain` in the
 *         remaining cases
 */
template<typename T, typename APPROX_TYPE>
TriBool operator==(const Polytope<T, APPROX_TYPE> &A,
                   const Bundle<T, APPROX_TYPE> &B);

/**
 * @brief Test whether two sets are the same
 *
 * @tparam T is the numeric type used to represent bundle
 * @tparam APPROX_TYPE is the type of approximation used in place of T
 *      during the computation
 * @param A is a bundle
 * @param B is a polytope
 * @return `true` when the method can establish that `A` and `B`
 *         represent the same set. `false` when it can
 *         establish `A` and `B` differ. `uncertain` in the
 *         remaining cases
 */
template<typename T, typename APPROX_TYPE>
TriBool operator==(const Bundle<T, APPROX_TYPE> &A,
                   const Polytope<T, APPROX_TYPE> &B);

/**
 * @brief Test whether two sets differ
 *
 * @tparam T is the numeric type used to represent bundle
 * @tparam APPROX_TYPE is the type of approximation used in place of T
 *      during the computation
 * @param A is a polytope
 * @param B is a parallelotope
 * @return `true` when the method can establish that `A` and `B`
 *         differ. `false` when it can establish `A` and `B`
 *         represent the same set. `uncertain` in the remaining
 *         cases
 */
template<typename T, typename APPROX_TYPE>
TriBool operator!=(const Polytope<T, APPROX_TYPE> &A,
                   const Bundle<T, APPROX_TYPE> &B);

/**
 * @brief Test whether two sets differ
 *
 * @tparam T is the numeric type used to represent bundle
 * @tparam APPROX_TYPE is the type of approximation used in place of T
 *      during the computation
 * @param A is a bundle
 * @param B is a polytope
 * @return `true` when the method can establish that `A` and `B`
 *         differ. `false` when it can establish `A` and `B`
 *         represent the same set. `uncertain` in the remaining
 *         cases
 */
template<typename T, typename APPROX_TYPE>
TriBool operator!=(const Bundle<T, APPROX_TYPE> &A,
                   const Polytope<T, APPROX_TYPE> &B);

/**
 * @brief Test whether two bundles differ
 *
 * @tparam T is the numeric type used to represent bundle
 * @tparam APPROX_TYPE is the type of approximation used in place of T
 *      during the computation
 * @param A is a bundle
 * @param B is a bundle
 * @return `true` when the method can establish that `A` and `B`
 *         differ. `false` when it can establish `A` and `B`
 *         represent the same set. `uncertain` in the remaining
 *         cases
 */
template<typename T, typename APPROX_TYPE>
TriBool operator!=(const Bundle<T, APPROX_TYPE> &A,
                   const Bundle<T, APPROX_TYPE> &B);

/**
 * @brief Get the intersection between two bundles
 *
 * @tparam T is the numeric type used to represent bundle
 * @tparam APPROX_TYPE is the type of approximation used in place of T
 *      during the computation
 * @param b1 is a bundle
 * @param b2 is a bundle
 * @return A bundle representing the intersection between
 *         the bundles `b1` and `b2`
 */
template<typename T, typename APPROX_TYPE>
Bundle<T, APPROX_TYPE> intersect(const Bundle<T, APPROX_TYPE> &b1,
                                 const Bundle<T, APPROX_TYPE> &b2);

/**
 * @brief Get the intersection between a bundle and a linear system
 *
 * @tparam T is the numeric type used to represent bundle
 * @tparam APPROX_TYPE is the type of approximation used in place of T
 *      during the computation
 * @param bundle is a bundle
 * @param linear_system is a linear system
 * @return A bundle representing the intersection between
 *         `bundle` and `linear_system`
 */
template<typename T, typename APPROX_TYPE>
Bundle<T, APPROX_TYPE>
intersect(const Bundle<T, APPROX_TYPE> &bundle,
          const LinearSystem<T, APPROX_TYPE> &linear_system);

/**
 * @brief Get the intersection between a linear system and a bundle
 *
 * @tparam T is the numeric type used to represent bundle
 * @tparam APPROX_TYPE is the type of approximation used in place of T
 *      during the computation
 * @param linear_system is a linear system
 * @param bundle is a bundle
 * @return A bundle representing the intersection between
 *         `linear_system` and `bundle`
 */
template<typename T, typename APPROX_TYPE>
Bundle<T, APPROX_TYPE>
intersect(const LinearSystem<T, APPROX_TYPE> &linear_system,
          const Bundle<T, APPROX_TYPE> &bundle);

/**
 * @brief Compute the over-approximation of the union of two bundles
 *
 * @tparam T is the numeric type used to represent bundle
 * @tparam APPROX_TYPE is the type of approximation used in place of T
 *      during the computation
 * @param b1 is a bundle
 * @param b2 is a bundle
 * @return A bundle over-approximating the union of `b1` and `b2`
 */
template<typename T, typename APPROX_TYPE>
Bundle<T, APPROX_TYPE>
over_approximate_union(const Bundle<T, APPROX_TYPE> &b1,
                       const Bundle<T, APPROX_TYPE> &b2);

// Inplementation

template<typename T, typename APPROX_TYPE>
Bundle<T, APPROX_TYPE>::Bundle()
{
}

template<typename T, typename APPROX_TYPE>
template<typename APPROX2>
Bundle<T, APPROX_TYPE>::Bundle(const Bundle<T, APPROX2> &orig):
    _directions(orig._directions),
    _adaptive_directions(orig._adaptive_directions),
    _lower_bounds(orig._lower_bounds), _upper_bounds(orig._upper_bounds),
    _templates(orig._templates)
{
}

template<typename T, typename APPROX_TYPE>
template<typename APPROX2>
Bundle<T, APPROX_TYPE>::Bundle(Bundle<T, APPROX2> &&orig)
{
  swap(*this, orig);
}

template<typename T, typename APPROX_TYPE>
Bundle<T, APPROX_TYPE>::Bundle(
    const std::set<size_t> &adaptive_directions,
    const std::vector<LinearAlgebra::Vector<T>> &directions,
    const LinearAlgebra::Vector<T> &lower_bounds,
    const LinearAlgebra::Vector<T> &upper_bounds,
    const std::set<BundleTemplate> &templates):
    _directions(directions),
    _adaptive_directions(adaptive_directions), _lower_bounds(lower_bounds),
    _upper_bounds(upper_bounds), _templates(templates)
{
}

template<typename T, typename APPROX_TYPE>
Bundle<T, APPROX_TYPE>::Bundle(
    const std::set<size_t> &adaptive_directions,
    std::vector<LinearAlgebra::Vector<T>> &&directions,
    LinearAlgebra::Vector<T> &&lower_bounds,
    LinearAlgebra::Vector<T> &&upper_bounds,
    const std::set<BundleTemplate> &templates):
    _directions(std::move(directions)),
    _adaptive_directions(adaptive_directions),
    _lower_bounds(std::move(lower_bounds)),
    _upper_bounds(std::move(upper_bounds)), _templates(templates)
{
}

template<typename T, typename APPROX_TYPE>
std::set<std::vector<unsigned int>>
Bundle<T, APPROX_TYPE>::duplicate_adaptive_directions(
    std::set<std::vector<unsigned int>> &templates)
{
  std::vector<bool> require_duplicate(_directions.size(), false);

  // build a new dynamic template set
  std::set<std::vector<unsigned int>> new_templates;

  // for each template in the dynamic template set
  for (auto bundle_template: templates) {
    bool changed = false;
    for (auto &idx: bundle_template) {

      // if one of template directions has been already mentioned
      if (_adaptive_directions.count(idx) > 0 && require_duplicate[idx]) {

        // the template must be changed
        changed = true;

        // add the index of the direction that we are going to
        // insert at the end of the direction vector among
        // the adaptive direction indices
        _adaptive_directions.insert(_directions.size());

        // copy the mentioned direction as a new direction
        _directions.emplace_back(_directions[idx]);
        _lower_bounds.push_back(_lower_bounds[idx]);
        _upper_bounds.push_back(_upper_bounds[idx]);

        // change the direction id in the template
        idx = _directions.size() - 1;
      } else { // if not yet mentioned, mark as to-be-duplicated
        require_duplicate[idx] = true;
      }
    }

    // if the bundle changed
    if (changed) {
      BundleTemplate::canonize(bundle_template);
    }
    new_templates.insert(bundle_template);
  }

  return new_templates;
}

template<typename T, typename APPROX_TYPE>
void Bundle<T, APPROX_TYPE>::validate_directions_and_bounds() const
{
  if (_directions.size() == 0) {
    SAPO_ERROR("direction vector must be non empty", std::domain_error);
  }

  const size_t dim = _directions[0].size();
  for (const auto &dir: _directions) {
    if (dir.size() != dim) {
      SAPO_ERROR("all the directions must have the same dimension",
                 std::domain_error);
    }

    if (LinearAlgebra::norm_infinity(dir) == 0) {
      SAPO_ERROR("all the directions must be non-null", std::domain_error);
    }
  }

  for (const auto &dir_index: _adaptive_directions) {
    if (dir_index >= _directions.size()) {
      SAPO_ERROR("the adaptive directions must be valid "
                 "indices for the direction vector",
                 std::domain_error);
    }
  }

  if (_directions.size() != _upper_bounds.size()) {
    SAPO_ERROR("directions and upper_bounds must have the same size",
               std::domain_error);
  }
  if (_directions.size() != _lower_bounds.size()) {
    SAPO_ERROR("directions and lower_bounds must have the same size",
               std::domain_error);
  }
}

template<typename T, typename APPROX_TYPE>
bool Bundle<T, APPROX_TYPE>::is_a_complete_basis(
    const std::vector<unsigned int> &bundle_template) const
{
  using namespace LinearAlgebra;

  Dense::Matrix<T> A;
  for (const auto &dir_idx: bundle_template) {
    if (dir_idx >= _directions.size()) {
      SAPO_ERROR("templates must contains valid "
                 "indices for the directions vectors",
                 std::domain_error);
    }

    A.emplace_back(_directions[dir_idx]);
  }

  return A.size() == dim() && LinearAlgebra::rank(A) == A.size();
}

template<typename T, typename APPROX_TYPE>
void Bundle<T, APPROX_TYPE>::validate_templates(
    const std::set<std::vector<unsigned int>> &templates) const
{
  for (const auto &bundle_template: templates) {
    if (bundle_template.size() != dim()) {
      SAPO_ERROR("templates must have as many columns "
                     << "as directions",
                 std::domain_error);
    }

    if (!is_a_complete_basis(bundle_template)) {
      SAPO_ERROR(bundle_template << " does not contain linearly "
                                 << "independent directions",
                 std::domain_error);
    }
  }
}

template<typename T, typename APPROX_TYPE>
void Bundle<T, APPROX_TYPE>::mark_linearly_dependent_directions(
    std::vector<size_t> &indices, const size_t &i) const
{
  using namespace LinearAlgebra;

  const bool dir_is_dyn = (_adaptive_directions.count(i) > 0);

  // compare the i-th direction all the directions in the interval
  // [i+1, |_directions|]
  for (size_t j = i + 1; j < _directions.size(); ++j) {
    if (indices[j] == j) { // the j-th direction is not dependent to
                           // any direction coming before j
      if ((_adaptive_directions.count(j) > 0) == dir_is_dyn) {
        // if the j-th direction is adaptive as well as the i-th direction
        // and they are linearly dependent
        if (are_linearly_dependent(_directions[i], _directions[j])) {

          // update the index of the j-th direction to be i
          indices[j] = i;
        }
      }
    }
  }
}

template<typename T, typename APPROX_TYPE>
std::vector<size_t>
Bundle<T, APPROX_TYPE>::filter_duplicated_directions()
{
  std::vector<LinearAlgebra::Vector<T>> new_directions;
  std::set<size_t> new_adaptive_dirs;
  LinearAlgebra::Vector<T> new_lower_bounds, new_upper_bounds;

  std::vector<size_t> new_pos(_directions.size());
  std::iota(std::begin(new_pos), std::end(new_pos), 0);

  for (size_t i = 0; i < _directions.size(); ++i) {

    // if directions[i] is the first instance of this direction
    // in directions
    if (new_pos[i] == i) {

      // get the new index for positions[i] in new_directions
      // and set in new_pos[i]
      const size_t new_idx = new_directions.size();
      new_pos[i] = new_idx;

      // add directions[i], lower_bounds[i], and
      // upper_bounds[i] in new_directions, new_lower_bounds,
      // and new_upper_bounds, respectively
      new_directions.push_back(_directions[i]);
      new_lower_bounds.push_back(_lower_bounds[i]);
      new_upper_bounds.push_back(_upper_bounds[i]);

      // if the i-th direction is dynamic add among the
      // adaptive directions of the new direction vector
      if (_adaptive_directions.count(i) > 0) {
        new_adaptive_dirs.insert(new_idx);
      }

      // mark in new_pos directions all direction vectors that are
      // linearly dependent to `_directions[i]` and are among dynamic
      // directions as `_direction[i]` is/is not.
      mark_linearly_dependent_directions(new_pos, i);
    } else {
      // directions[i] has been already inserted in new_directions
      // in position new_pos[i]
      using namespace LinearAlgebra;

      // get the index of directions[i] in new_directions
      const unsigned int new_idx = new_pos[i];

      const APPROX_TYPE coeff = get_dependency_coefficient<T, APPROX_TYPE>(
          new_directions[new_idx], _directions[i]);

      // if necessary, increase the lower bound
      new_lower_bounds[new_idx]
          = std::max(new_lower_bounds[new_idx],
                     ::get_lower_bound(_lower_bounds[i] * coeff));

      // if necessary, decrease the upper bound
      new_upper_bounds[new_idx]
          = std::min(new_upper_bounds[new_idx],
                     ::get_upper_bound(_upper_bounds[i] * coeff));
    }
  }

  // replace old directions, adaptive_directions, lower_bounds, and
  // upper_bounds
  std::swap(_directions, new_directions);
  std::swap(_adaptive_directions, new_adaptive_dirs);
  std::swap(_lower_bounds, new_lower_bounds);
  std::swap(_upper_bounds, new_upper_bounds);

  return new_pos;
}

template<typename T, typename APPROX_TYPE>
void Bundle<T, APPROX_TYPE>::resort_directions(
    const std::map<unsigned int, unsigned int> &new_positions)
{
  std::vector<LinearAlgebra::Vector<T>> new_directions(new_positions.size());
  std::set<size_t> new_adaptive_dirs;
  LinearAlgebra::Vector<T> new_lower_bounds(new_positions.size()),
      new_upper_bounds(new_positions.size());

  for (const auto &new_pos: new_positions) {
    std::swap(_directions[new_pos.first], new_directions[new_pos.second]);
    std::swap(_lower_bounds[new_pos.first], new_lower_bounds[new_pos.second]);
    std::swap(_upper_bounds[new_pos.first], new_upper_bounds[new_pos.second]);
  }

  for (const auto &adaptive_direction: _adaptive_directions) {
    if (new_positions.find(adaptive_direction) != new_positions.end()) {
      new_adaptive_dirs.insert(new_positions.at(adaptive_direction));
    }
  }

  std::swap(_directions, new_directions);
  std::swap(_adaptive_directions, new_adaptive_dirs);
  std::swap(_lower_bounds, new_lower_bounds);
  std::swap(_upper_bounds, new_upper_bounds);
}

template<typename T, typename APPROX_TYPE>
std::map<unsigned int, unsigned int>
Bundle<T, APPROX_TYPE>::find_new_position_for(
    const std::set<unsigned int> &positions)
{
  std::map<unsigned int, unsigned int> new_positions;
  unsigned int new_id = 0;

  // for each of the positions in the set of the positions to
  // be reassigned
  for (const auto &dir_id: positions) {

    // assign a new position
    new_positions[dir_id] = new_id++;
  }

  return new_positions;
}

template<typename T, typename APPROX_TYPE>
Bundle<T, APPROX_TYPE>::Bundle(
    const std::vector<LinearAlgebra::Vector<T>> &directions,
    const LinearAlgebra::Vector<T> &lower_bounds,
    const LinearAlgebra::Vector<T> &upper_bounds,
    std::set<std::vector<unsigned int>> templates,
    const std::set<size_t> &adaptive_directions,
    const bool remove_unused_directions):
    _directions(directions),
    _adaptive_directions(adaptive_directions), _lower_bounds(lower_bounds),
    _upper_bounds(upper_bounds), _templates()
{
  validate_directions_and_bounds();
  validate_templates(templates);

  // filter duplicated direction, update templates, and bring them in canonical
  // form
  auto new_pos = filter_duplicated_directions();

  templates = BundleTemplate::update_directions(templates, new_pos);

  auto used_dirs = BundleTemplate::collect_directions(templates);

  if (!remove_unused_directions) {
    // add missing directions in templates
    std::set<size_t> missing_dirs;
    for (size_t idx = 0; idx < _directions.size(); ++idx) {
      if (used_dirs.count(idx) == 0) {
        missing_dirs.insert(idx);
      }
    }
    BundleTemplate::add_missing_templates(templates, _directions,
                                          missing_dirs);
  } else {
    if (templates.size() == 0) {
      SAPO_ERROR("template vector must be non empty", std::domain_error);
    }

    // if the number of directions appearing in the templates is smaller
    // than the number of directions
    if (used_dirs.size() != directions.size()) {

      // find a new position for the used directions
      auto new_pos = find_new_position_for(used_dirs);

      resort_directions(new_pos);
      templates = BundleTemplate::update_directions(templates, new_pos);
    }
  }
  templates = duplicate_adaptive_directions(templates);

  for (auto &bundle_template: templates) {
    _templates.emplace(bundle_template, _adaptive_directions);
  }
}

template<typename T, typename APPROX_TYPE>
Bundle<T, APPROX_TYPE>::Bundle(
    const std::vector<LinearAlgebra::Vector<T>> &directions,
    const LinearAlgebra::Vector<T> &lower_bounds,
    const LinearAlgebra::Vector<T> &upper_bounds,
    const std::set<std::vector<unsigned int>> &templates,
    const bool remove_unused_directions):
    Bundle<T, APPROX_TYPE>(directions, lower_bounds, upper_bounds, templates, {},
           remove_unused_directions)
{
}

template<typename T, typename APPROX_TYPE>
Bundle<T, APPROX_TYPE>::Bundle(
    const std::vector<LinearAlgebra::Vector<T>> &directions,
    const LinearAlgebra::Vector<T> &lower_bounds,
    const LinearAlgebra::Vector<T> &upper_bounds,
    const bool remove_unused_directions):
    Bundle<T, APPROX_TYPE>(directions, lower_bounds, upper_bounds, {}, {},
           remove_unused_directions)
{
}

template<typename T, typename APPROX_TYPE>
template<typename APPROX2>
Bundle<T, APPROX_TYPE>::Bundle(const Polytope<T, APPROX2> &P)
{
  using namespace LinearAlgebra;

  std::vector<size_t> double_row(P.size(), false);

  for (size_t i = 0; i < P.size(); ++i) {
    if (!double_row[i]) {
      const auto &Ai = P.A(i);
      const auto nAi = -P.A(i);
      for (size_t j = i + 1; j < P.size(); ++j) {
        if (!double_row[j] && (P.A(j) == Ai || P.A(j) == nAi)) {
          double_row[j] = true;
        }
      }
      _directions.push_back(P.A(i));
      _upper_bounds.push_back(P.maximize(P.A(i)).objective_value());
      _lower_bounds.push_back(P.minimize(P.A(i)).objective_value());
    }
  }

  std::set<size_t> missing_dirs;
  for (size_t i = 0; i < _directions.size(); ++i) {
    missing_dirs.insert(i);
  }

  BundleTemplate::add_missing_templates(_templates, _directions, {},
                                        missing_dirs);
}

template<typename T, typename APPROX_TYPE>
TriBool
Bundle<T, APPROX_TYPE>::satisfies(const LinearSystem<T, APPROX_TYPE> &ls) const
{
  // if this object is empty
  if (is_true(this->is_empty())) {
    return true;
  }

  // if the parameter has no solutions and this object is not,
  // return false
  if (is_false(ls.has_solutions())) {
    return false;
  }

  TriBool sat{true};

  Polytope<double> P_this = *this;
  // for each direction in the bundle
  for (unsigned int dir_idx = 0; dir_idx < ls.size(); ++dir_idx) {

    auto max_res = P_this.maximize(ls.A(dir_idx));

    sat = sat && (max_res.objective_value() <= ls.b(dir_idx));
    if (is_false(sat)) {
      return false;
    }
  }

  return sat;
}

template<typename T, typename APPROX_TYPE>
std::list<Bundle<T, APPROX_TYPE>> &
Bundle<T, APPROX_TYPE>::split_bundle(std::list<Bundle<T, APPROX_TYPE>> &res,
                                     LinearAlgebra::Vector<T> &lower_bounds,
                                     LinearAlgebra::Vector<T> &upper_bounds,
                                     const size_t idx, const T &max_magnitude,
                                     const T &split_ratio) const
{
  if (idx == this->size()) {
    Bundle new_bundle(this->_adaptive_directions, this->_directions,
                      this->_lower_bounds, this->_upper_bounds,
                      this->_templates);
    res.push_back(std::move(new_bundle));

    return res;
  }

  if (std::abs(_upper_bounds[idx] - _lower_bounds[idx]) > max_magnitude) {
    T lower_bound = _lower_bounds[idx];

    do {
      const T upper_bound = std::min(lower_bound + split_ratio * max_magnitude,
                                     _upper_bounds[idx]);

      lower_bounds[idx] = lower_bound;
      upper_bounds[idx] = upper_bound;
      this->split_bundle(res, lower_bounds, upper_bounds, idx + 1,
                         max_magnitude, split_ratio);

      lower_bound = upper_bound;
    } while (_upper_bounds[idx] != lower_bound);
  } else {
    lower_bounds[idx] = _lower_bounds[idx];
    upper_bounds[idx] = _upper_bounds[idx];
    this->split_bundle(res, lower_bounds, upper_bounds, idx + 1, max_magnitude,
                       split_ratio);
  }
  return res;
}

template<typename T, typename APPROX_TYPE>
template<typename APPROX2>
Bundle<T, APPROX_TYPE> &
Bundle<T, APPROX_TYPE>::operator=(const Bundle<T, APPROX2> &orig)
{
  this->_directions = orig._directions;
  this->_adaptive_directions = orig._adaptive_directions;
  this->_templates = orig._templates;

  this->_upper_bounds = orig._upper_bounds;
  this->_lower_bounds = orig._lower_bounds;

  return *this;
}

template<typename T, typename APPROX_TYPE>
template<typename APPROX2>
Bundle<T, APPROX_TYPE> &
Bundle<T, APPROX_TYPE>::operator=(Bundle<T, APPROX2> &&orig)
{
  swap(*this, orig);

  return *this;
}

template<typename T, typename APPROX_TYPE>
Bundle<T, APPROX_TYPE>::operator Polytope<T, APPROX_TYPE>() const
{
  using namespace std;
  using namespace LinearAlgebra;

  std::vector<Vector<T>> A;
  Vector<T> b;
  for (unsigned int i = 0; i < this->size(); i++) {
    A.push_back(this->_directions[i]);
    b.push_back(this->_upper_bounds[i]);
  }
  for (unsigned int i = 0; i < this->size(); i++) {
    A.push_back(-this->_directions[i]);
    b.push_back(AVOID_NEG_ZERO(-this->_lower_bounds[i]));
  }

  return Polytope<T, APPROX_TYPE>(std::move(A), std::move(b));
}

template<typename T, typename APPROX_TYPE>
Parallelotope<T, APPROX_TYPE> Bundle<T, APPROX_TYPE>::get_parallelotope(
    const BundleTemplate &bundle_template) const
{
  using namespace std;

  vector<T> lbound, ubound;
  vector<LinearAlgebra::Vector<T>> Lambda;

  auto it = std::begin(bundle_template.direction_indices());
  // upper facets
  for (size_t j = 0; j < this->dim(); j++) {
    const unsigned int idx = *it;
    if (idx >= this->_directions.size()) {
      SAPO_ERROR("the parameter is not a template for the bundle",
                 std::domain_error);
    }
    Lambda.push_back(this->_directions[idx]);
    ubound.push_back(this->_upper_bounds[idx]);
    lbound.push_back(this->_lower_bounds[idx]);

    ++it;
  }

  return Parallelotope<T, APPROX_TYPE>(Lambda, lbound, ubound);
}

template<typename T, typename APPROX_TYPE>
Bundle<T, APPROX_TYPE> Bundle<T, APPROX_TYPE>::get_canonical() const
{
  Bundle res(*this);

  res.canonize();

  return res;
}

template<typename T, typename APPROX_TYPE>
Bundle<T, APPROX_TYPE> &Bundle<T, APPROX_TYPE>::canonize()
{
  if (this->size() == 0) {
    return *this;
  }

  // if the bundle is empty
  if (is_true(is_empty())) {
    return *this;
  }

  // get current polytope
  Polytope<T, APPROX_TYPE> bund = *this;
  for (unsigned int i = 0; i < this->size(); ++i) {
    _lower_bounds[i] = bund.minimize(this->_directions[i]).objective_value();
    _upper_bounds[i] = bund.maximize(this->_directions[i]).objective_value();
  }
  return *this;
}

template<typename T, typename APPROX_TYPE>
std::list<Bundle<T, APPROX_TYPE>>
Bundle<T, APPROX_TYPE>::split(const T max_magnitude, const T split_ratio) const
{
  std::list<Bundle<T, APPROX_TYPE>> split_list;

  LinearAlgebra::Vector<T> upper_bounds(this->size());
  LinearAlgebra::Vector<T> lower_bounds(this->size());

  this->split_bundle(split_list, lower_bounds, upper_bounds, 0, max_magnitude,
                     split_ratio);

  return split_list;
}

template<typename T, typename APPROX_TYPE>
size_t Bundle<T, APPROX_TYPE>::get_a_linearly_dependent_direction(
    const LinearAlgebra::Vector<T> &v) const
{
  for (size_t i = 0; i < _directions.size(); ++i) {
    using namespace LinearAlgebra;
    if (are_linearly_dependent(v, _directions[i])) {
      return i;
    }
  }
  return _directions.size();
}

template<typename T, typename APPROX_TYPE>
void Bundle<T, APPROX_TYPE>::add_templates(
    const std::set<BundleTemplate> &new_templates,
    const std::vector<unsigned int> &new_direction_indices)
{
  // check whether some of the templates must be copied
  for (const BundleTemplate &temp: new_templates) {
    // if this is the case, copy and update the template indices
    std::vector<unsigned int> t_copy(temp.direction_indices());
    for (unsigned int j = 0; j < t_copy.size(); ++j) {
      t_copy[j] = new_direction_indices[t_copy[j]];
    }

    // add the new template to the intersected bundle
    _templates.emplace(t_copy, std::set<size_t>{});
  }
}

template<typename T, typename APPROX_TYPE>
LinearAlgebra::Vector<T> Bundle<T, APPROX_TYPE>::edge_lengths() const
{
  using namespace LinearAlgebra;

  Vector<double> dist(this->size());
  for (unsigned int i = 0; i < this->size(); i++) {
    dist[i] = std::abs(this->_upper_bounds[i] - this->_lower_bounds[i])
              / norm_2(this->_directions[i]);
  }
  return dist;
}

template<typename T, typename APPROX_TYPE>
Bundle<T, APPROX_TYPE> &
Bundle<T, APPROX_TYPE>::intersect_with(const Bundle<T, APPROX_TYPE> &A)
{
  using namespace LinearAlgebra;

  if (dim() != A.dim()) {
    SAPO_ERROR("the two bundles differ in dimensions", std::domain_error);
  }

  std::vector<unsigned int> new_ids(A.size());

  // for each direction in A
  for (size_t i = 0; i < A.size(); ++i) {
    const Vector<T> &A_dir = A._directions[i];
    new_ids[i] = get_a_linearly_dependent_direction(A_dir);

    // if the direction is not present in this object
    if (new_ids[i] == this->size()) {

      // add the direction and the corresponding boundaries
      this->_directions.push_back(A_dir);
      this->_lower_bounds.push_back(A._lower_bounds[i]);
      this->_upper_bounds.push_back(A._upper_bounds[i]);
    } else { // if the direction is already included in this object

      // compute the dependency coefficient
      const auto &new_i = new_ids[i];
      const APPROX_TYPE dep_coeff = get_dependency_coefficient<T, APPROX_TYPE>(
          this->_directions[new_i], A_dir);

      T scaled_A_lb, scaled_A_ub;
      if (is_true(dep_coeff > 0)) {
        scaled_A_lb = ::get_lower_bound(A._lower_bounds[i] * dep_coeff);
        scaled_A_ub = ::get_upper_bound(A._upper_bounds[i] * dep_coeff);
      } else {
        scaled_A_lb = ::get_lower_bound(A._upper_bounds[i] * dep_coeff);
        scaled_A_ub = ::get_upper_bound(A._lower_bounds[i] * dep_coeff);
      }

      // increase the lower bound if necessary
      if (scaled_A_lb > this->_lower_bounds[new_i]) {
        this->_lower_bounds[new_i] = scaled_A_lb;
      }

      // decrease the upper bound if necessary
      if (scaled_A_ub < this->_upper_bounds[new_i]) {
        this->_upper_bounds[new_i] = scaled_A_ub;
      }
    }
  }

  add_templates(A.templates(), new_ids);

  return *this;
}

template<typename T, typename APPROX_TYPE>
Bundle<T, APPROX_TYPE> &
Bundle<T, APPROX_TYPE>::intersect_with(const LinearSystem<T, APPROX_TYPE> &ls)
{
  using namespace LinearAlgebra;

  // If the linear system has no constraint or involves no variables
  if (ls.dim() == 0) {
    return *this;
  }

  if (dim() != ls.dim()) {
    SAPO_ERROR("the bundle and the linear system differ in "
               "dimensions",
               std::domain_error);
  }

  std::set<size_t> outside_templates;
  const Dense::Matrix<T> &A = ls.A();

  std::vector<unsigned int> new_ids(A.size());
  // for each row in the linear system
  for (size_t i = 0; i < A.size(); ++i) {
    const LinearAlgebra::Vector<double> &A_row = A[i];

    new_ids[i] = get_a_linearly_dependent_direction(A_row);

    // if the direction is not present in this object
    if (new_ids[i] == this->size()) {

      // add the direction and the corresponding boundaries
      outside_templates.insert(new_ids[i]);
      this->_directions.push_back(A_row);
      this->_lower_bounds.push_back(-std::numeric_limits<T>::infinity());
      this->_upper_bounds.push_back(ls.b(i));
    } else {
      // compute the dependency coefficient
      const unsigned int &new_i = new_ids[i];
      const APPROX_TYPE dep_coeff = get_dependency_coefficient<T, APPROX_TYPE>(
          this->_directions[new_i], A_row);

      auto b_rescaled = ls.b(i) * dep_coeff;
      if (dep_coeff > 0) {
        // if necessary, decrease the upper bound
        if (this->_upper_bounds[new_i] > ::get_upper_bound(b_rescaled)) {
          this->_upper_bounds[new_i] = ::get_upper_bound(b_rescaled);
        }
      } else {
        // or increase the lower bound
        if (this->_lower_bounds[new_i] < ::get_lower_bound(b_rescaled)) {
          this->_lower_bounds[new_i] = ::get_lower_bound(b_rescaled);
        }
      }
    }
  }

  // add missing templates
  BundleTemplate::add_missing_templates(
      _templates, _directions, _adaptive_directions, outside_templates);

  return this->canonize();
}

template<typename T, typename APPROX_TYPE>
TriBool Bundle<T, APPROX_TYPE>::is_subset_of(
    const Bundle<T, APPROX_TYPE> &bundle) const
{
  // if this object is empty
  if (is_true(this->is_empty())) {
    return true;
  }

  // if the parameter is empty and this object is not,
  // return false
  if (is_true(bundle.is_empty())) {
    return false;
  }

  TriBool is_sub{true};

  Polytope<double> P_this = *this;
  // for each direction in the bundle
  for (unsigned int dir_idx = 0; dir_idx < bundle.size(); ++dir_idx) {

    auto min_value = P_this.minimize(bundle.get_direction(dir_idx));

    // if the minimum of this object on that direction is lesser than
    // the bundle minimum, this object is not a subset of the bundle
    is_sub
        = is_sub
          && (min_value.objective_value() >= bundle.get_lower_bound(dir_idx));
    if (is_false(is_sub)) {
      return is_sub;
    }

    auto max_value = P_this.maximize(bundle.get_direction(dir_idx));

    // if the maximum of this object on that direction is greater than
    // the bundle maximum, this object is not a subset of the bundle
    is_sub
        = is_sub
          && (max_value.objective_value() <= bundle.get_upper_bound(dir_idx));
    if (is_false(is_sub)) {
      return is_sub;
    }
  }

  return is_sub;
}

template<typename T, typename APPROX_TYPE>
Bundle<T, APPROX_TYPE> &Bundle<T, APPROX_TYPE>::expand_by(const T delta)
{
  if (is_true(is_empty())) {
    return *this;
  }

  for (size_t i = 0; i < _directions.size(); ++i) {
    using namespace LinearAlgebra;

    const auto dir_delta = norm_2(_directions[i]) * delta;
    _lower_bounds[i] = ::get_lower_bound(_lower_bounds[i] - dir_delta);
    _upper_bounds[i] = ::get_upper_bound(_upper_bounds[i] + dir_delta);
  }

  return *this;
}

template<typename T, typename APPROX_TYPE>
Bundle<T, APPROX_TYPE>::~Bundle()
{
}

template<typename T, typename APPROX_TYPE>
TriBool are_disjoint(const Bundle<T, APPROX_TYPE> &A,
                     const Bundle<T, APPROX_TYPE> &B)
{
  if (A.dim() != B.dim()) {
    SAPO_ERROR("the two bundles must have the same dimension",
               std::domain_error);
  }

  if (A.dim() == 0) {
    return false;
  }

  using namespace LinearAlgebra;

  Dense::Matrix<T> ls_A;
  Vector<T> ls_b;

  auto add_bundle_constraints
      = [&ls_A, &ls_b](const Bundle<T, APPROX_TYPE> &bundle) {
          for (size_t i = 0; i < bundle.size(); ++i) {
            ls_A.push_back(bundle.get_direction(i));
            ls_b.push_back(bundle.get_upper_bound(i));
            ls_A.push_back(-bundle.get_direction(i));
            ls_b.push_back(-bundle.get_lower_bound(i));
          }
        };

  add_bundle_constraints(A);
  add_bundle_constraints(B);

  SimplexMethodOptimizer<T, APPROX_TYPE> optimizer;

  return !optimizer.is_feasible(ls_A, ls_b);
}

template<typename T, typename APPROX1, typename APPROX2>
void swap(Bundle<T, APPROX1> &A, Bundle<T, APPROX2> &B)
{
  std::swap(A._directions, B._directions);
  std::swap(A._adaptive_directions, B._adaptive_directions);
  std::swap(A._upper_bounds, B._upper_bounds);
  std::swap(A._lower_bounds, B._lower_bounds);
  std::swap(A._templates, B._templates);
}

template<typename T, typename APPROX_TYPE>
inline TriBool operator==(const Bundle<T, APPROX_TYPE> &b1,
                          const Bundle<T, APPROX_TYPE> &b2)
{
  return b1.is_subset_of(b2) && b2.is_subset_of(b1);
}

template<typename T, typename APPROX_TYPE>
inline TriBool operator==(const Polytope<T, APPROX_TYPE> &A,
                          const Bundle<T, APPROX_TYPE> &B)
{
  return static_cast<Polytope<double>>(B) == A;
}

template<typename T, typename APPROX_TYPE>
inline TriBool operator==(const Bundle<T, APPROX_TYPE> &A,
                          const Polytope<T, APPROX_TYPE> &B)
{
  return B == A;
}

template<typename T, typename APPROX_TYPE>
inline TriBool operator!=(const Polytope<T, APPROX_TYPE> &A,
                          const Bundle<T, APPROX_TYPE> &B)
{
  return !(A == B);
}

template<typename T, typename APPROX_TYPE>
inline TriBool operator!=(const Bundle<T, APPROX_TYPE> &A,
                          const Polytope<T, APPROX_TYPE> &B)
{
  return !(A == B);
}

template<typename T, typename APPROX_TYPE>
inline TriBool operator!=(const Bundle<T, APPROX_TYPE> &A,
                          const Bundle<T, APPROX_TYPE> &B)
{
  return !(A == B);
}

template<typename T, typename APPROX_TYPE>
inline Bundle<T, APPROX_TYPE> intersect(const Bundle<T, APPROX_TYPE> &b1,
                                        const Bundle<T, APPROX_TYPE> &b2)
{
  Bundle res(b1);

  res.intersect_with(b2);

  return res;
}

template<typename T, typename APPROX_TYPE>
inline Bundle<T, APPROX_TYPE>
intersect(const Bundle<T, APPROX_TYPE> &bundle,
          const LinearSystem<T, APPROX_TYPE> &linear_system)
{
  Bundle res(bundle);

  res.intersect_with(linear_system);

  return res;
}

template<typename T, typename APPROX_TYPE>
inline Bundle<T, APPROX_TYPE>
intersect(const LinearSystem<T, APPROX_TYPE> &linear_system,
          const Bundle<T, APPROX_TYPE> &bundle)
{
  return intersect(bundle, linear_system);
}

template<typename T, typename APPROX_TYPE>
Bundle<T, APPROX_TYPE> over_approximate_union(const Bundle<T, APPROX_TYPE> &b1,
                                              const Bundle<T, APPROX_TYPE> &b2)
{
  using namespace LinearAlgebra;
  using namespace LinearAlgebra::Dense;

  if (b1.dim() != b2.dim()) {
    SAPO_ERROR("the two bundles differ in dimensions", std::domain_error);
  }

  if (is_true(b1.is_empty())) {
    return b2;
  }

  if (is_true(b2.is_empty())) {
    return b1;
  }

  Polytope<T, APPROX_TYPE> p2(b2);
  Bundle<T, APPROX_TYPE> res(b1);
  const Matrix<T> &res_dirs = res.directions();

  // Updates res boundaries to include p2
  for (unsigned int i = 0; i < res_dirs.size(); ++i) {
    const LinearAlgebra::Vector<double> &res_dir = res_dirs[i];

    const auto res_min = p2.minimize(res_dir).objective_value();
    res._lower_bounds[i]
        = std::min(::get_lower_bound(res_min), res._lower_bounds[i]);
    const auto res_max = p2.maximize(res_dir).objective_value();
    res._upper_bounds[i]
        = std::max(::get_upper_bound(res_max), res._upper_bounds[i]);
  }

  const Matrix<T> &b2_dirs = b2.directions();

  Polytope<T, APPROX_TYPE> p1(b1);
  std::vector<unsigned int> new_ids(b2.size());
  // for each row in the linear system
  for (size_t i = 0; i < b2_dirs.size(); ++i) {
    const LinearAlgebra::Vector<T> &b2_dir = b2_dirs[i];

    new_ids[i] = res.get_a_linearly_dependent_direction(b2_dir);

    // if the direction is not present in this object
    if (new_ids[i] == res.size()) {
      const auto dir_min = p1.minimize(b2_dir).objective_value();
      const T lower_bound
          = std::min(get_lower_bound(dir_min), b2._lower_bounds[i]);
      const auto dir_max = p1.maximize(b2_dir).objective_value();
      const T upper_bound
          = std::max(get_lower_bound(dir_max), b2._upper_bounds[i]);

      // add the direction and the corresponding boundaries
      res._directions.push_back(b2_dir);
      res._lower_bounds.push_back(lower_bound);
      res._upper_bounds.push_back(upper_bound);
    } else {
      // compute the dependency coefficient
      const unsigned int &new_i = new_ids[i];
      const APPROX_TYPE dep_coeff = get_dependency_coefficient<T, APPROX_TYPE>(
          res._directions[new_i], b2_dir);

      T lb_rescaled = ::get_lower_bound(b2._lower_bounds[i] * dep_coeff);
      T ub_rescaled = ::get_upper_bound(b2._upper_bounds[i] * dep_coeff);
      if (dep_coeff < 0) {
        std::swap(lb_rescaled, ub_rescaled);
      }

      // if necessary, update the boundaries
      if (res._upper_bounds[new_i] < ub_rescaled) {
        res._upper_bounds[new_i] = ub_rescaled;
      }
      if (res._lower_bounds[new_i] > lb_rescaled) {
        res._lower_bounds[new_i] = lb_rescaled;
      }
    }
  }

  res.add_templates(b2.templates(), new_ids);

  return res;
}

template<typename T, typename APPROX_TYPE>
SetsUnion<Bundle<T, APPROX_TYPE>>
subtract_and_close(const Bundle<T, APPROX_TYPE> &b1,
                   const Bundle<T, APPROX_TYPE> &b2)
{
  SetsUnion<Bundle<T, APPROX_TYPE>> su;
  if (is_true(b2.includes(b1))) {
    return su;
  }

  if (is_true(are_disjoint(b1, b2))) {
    su.add(b1);

    return su;
  }

  Polytope<T, APPROX_TYPE> p1 = b1;

  for (size_t i = 0; i < b2.size(); ++i) {
    const auto b2_max = p1.maximize(b2._directions[i]).objective_value();
    if (::get_upper_bound(b2_max) > b2._upper_bounds[i]) {
      Bundle<T, APPROX_TYPE> new_b1{b1};
      new_b1._directions.push_back(b2._directions[i]);
      new_b1._lower_bounds.push_back(b2._upper_bounds[i]);
      new_b1._upper_bounds.push_back(::get_upper_bound(b2_max));

      su.add(std::move(new_b1.canonize()));
    }

    const auto b2_min = p1.minimize(b2._directions[i]).objective_value();
    if (::get_lower_bound(b2_min) < b2._lower_bounds[i]) {
      Bundle<T, APPROX_TYPE> new_b1{b1};
      new_b1._directions.push_back(b2._directions[i]);
      new_b1._lower_bounds.push_back(::get_lower_bound(b2_min));
      new_b1._upper_bounds.push_back(b2._lower_bounds[i]);

      su.add(std::move(new_b1.canonize()));
    }
  }

  return su;
}

#endif /* BUNDLE_H_ */
