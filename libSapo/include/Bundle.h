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
#include "Bernstein.h"
#include "Polytope.h"
#include "Parallelotope.h"

#include "STL/Atom.h"

#define SPLIT_MAGNITUDE_RATIO                                                 \
  0.75 //!< define the default versor magnitude multiplier for bundle splits

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
 * @brief A class for parallelotope bundles
 *
 * A parallelotope bundle represents the intersection between different
 * non-singular parallelotopes. This class stores all of the parallelotope
 * directions/axes in one single array and constraints each of them between an
 * upper and a lower bound. The parallelotopes are represented by means of
 * bundle templates. A template is the vector of the indices of the bundle
 * directions involved in the corresponding parallelotope.
 */
class Bundle;

/**
 * @brief Bundle template
 *
 * This class represents bundle templates. A template of a $n$-dimensional
 * bundle is a vector of the indices of $n$ linearly independent directions of
 * the bundle itself. A bundle consists in the intersection of parallelotopes.
 * A template symbolically represents one of these parallelotopes.
 */
class BundleTemplate
{
  std::vector<unsigned int>
      _dir_indices; //!< Indices of the template directions

public:
  /**
   * @brief A copy constructor
   *
   * @param[in] orig is the copied bundle template
   */
  BundleTemplate(const BundleTemplate &orig);

  /**
   * @brief A move constructor
   *
   * @param[in] orig is the moved bundle template
   */
  BundleTemplate(BundleTemplate &&orig);

  /**
   * @brief A constructor for bundle templates
   *
   * @param[in] dir_indices is the vector of linearly independent directions of
   * the template
   */
  BundleTemplate(const std::vector<unsigned int> &dir_indices);

  /**
   * @brief A move constructor for bundle templates
   *
   * @param[in] dir_indices is the vector of linearly independent directions of
   * the template
   */
  BundleTemplate(std::vector<unsigned int> &&dir_indices);

  /**
   * @brief A move operator
   *
   * @param orig is the original template
   * @return a reference to the updated template
   */
  inline BundleTemplate &operator=(BundleTemplate &&orig)
  {
    std::swap(_dir_indices, orig._dir_indices);

    return *this;
  }

  /**
   * @brief A copy operator
   *
   * @param orig is the original template
   * @return a reference to the updated template
   */
  inline BundleTemplate &operator=(const BundleTemplate &orig)
  {
    _dir_indices = orig._dir_indices;

    return *this;
  }

  /**
   * @brief Get the vector of the template direction indices
   *
   * @return the vector of the template direction indices
   */
  inline const std::vector<unsigned int> &direction_indices() const
  {
    return _dir_indices;
  }

  /**
   * @brief Get the number of direction indices in the template
   *
   * @return the number of direction indices in the template
   */
  inline size_t dim() const
  {
    return _dir_indices.size();
  }

  /**
   * @brief Get the i-th direction index in the template
   *
   * @param i is the position of the aimed direction index
   * @return a reference to the i-th direction index in the
   *        template
   */
  inline const unsigned int &operator[](const size_t &i) const
  {
    return _dir_indices[i];
  }
};

/**
 * @brief Print a bundle template
 *
 * @param os is the output stream
 * @param bundle_template is the bundle template to be printed
 * @return a reference to the output stream
 */
std::ostream &operator<<(std::ostream &os,
                         const BundleTemplate &bundle_template);

template<>
struct std::less<BundleTemplate> {
  bool operator()(const BundleTemplate &a, const BundleTemplate &b) const;
};

class Bundle
{
  std::vector<LinearAlgebra::Vector<double>>
      _directions; //!< the vector of directions
  std::set<size_t>
      _dynamic_directions; //!< the set of dynamic direction indices
  LinearAlgebra::Vector<double> _lower_bounds; //!< direction upper bounds
  LinearAlgebra::Vector<double> _upper_bounds; //!< direction lower bounds
  std::set<BundleTemplate> _templates;         //!< bundle templates

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
  std::list<Bundle> &split_bundle(std::list<Bundle> &res,
                                  LinearAlgebra::Vector<double> &lower_bounds,
                                  LinearAlgebra::Vector<double> &upper_bounds,
                                  const size_t idx,
                                  const double &max_magnitude,
                                  const double &split_ratio) const;

  /**
   * @brief A constructor
   *
   * This constructor does not perform any check on the parameter
   * consistency.
   *
   * @param[in] dynamic_directions is the set of dynamic direction indices
   * @param[in] directions is the direction vector
   * @param[in] lower_bounds is the vector of direction lower bounds
   * @param[in] upper_bounds is the vector of direction upper bounds
   * @param[in] templates is the set of the bundle templates
   */
  Bundle(const std::set<size_t> &dynamic_directions,
         const std::vector<LinearAlgebra::Vector<double>> &directions,
         const LinearAlgebra::Vector<double> &lower_bounds,
         const LinearAlgebra::Vector<double> &upper_bounds,
         const std::set<BundleTemplate> &templates);

  /**
   * @brief A constructor
   *
   * This constructor does not perform any check on the parameter
   * consistency.
   *
   * @param[in] dynamic_directions is the set of dynamic direction indices
   * @param[in] directions is the direction vector
   * @param[in] lower_bounds is the vector of direction lower bounds
   * @param[in] upper_bounds is the vector of direction upper bounds
   * @param[in] templates is the set of the bundle templates
   */
  Bundle(const std::set<size_t> &dynamic_directions,
         std::vector<LinearAlgebra::Vector<double>> &&directions,
         LinearAlgebra::Vector<double> &&lower_bounds,
         LinearAlgebra::Vector<double> &&upper_bounds,
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
   * @param[in] templates is the set of the bundle templates
   * @param[in] dynamic_direction is the set of dynamic directions
   * @param[in] remove_unused_directions is a flag to remove
   *      directions not belonging to any template
   */
  Bundle(const std::vector<LinearAlgebra::Vector<double>> &directions,
         const LinearAlgebra::Vector<double> &lower_bounds,
         const LinearAlgebra::Vector<double> &upper_bounds,
         std::set<std::vector<unsigned int>> templates,
         const std::set<size_t> &dynamic_directions,
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
  Bundle(const std::vector<LinearAlgebra::Vector<double>> &directions,
         const LinearAlgebra::Vector<double> &lower_bounds,
         const LinearAlgebra::Vector<double> &upper_bounds,
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
  Bundle(const std::vector<LinearAlgebra::Vector<double>> &directions,
         const LinearAlgebra::Vector<double> &lower_bounds,
         const LinearAlgebra::Vector<double> &upper_bounds,
         const bool remove_unused_directions = false);

  /**
   * @brief A constructor
   *
   * Build the bundle representing a polytope
   *
   * @param P is the polytope whose bundle representation is aimed
   */
  explicit Bundle(const Polytope &P);

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
  inline const std::vector<LinearAlgebra::Vector<double>> &directions() const
  {
    return _directions;
  }

  /**
   * @brief Get the i-th direction in the bundle
   *
   * @param i is the index of the aimed direction
   * @return  a reference to the i-th direction in the bundle
   */
  inline const LinearAlgebra::Vector<double> &
  get_direction(const size_t &i) const
  {
    return _directions[i];
  }

  /**
   * @brief Check whether the i-th direction is dynamic
   *
   * @param i is the index of the considered direction
   * @return `true` if and only if the i-th direction
   *      is dynamic
   */
  inline bool is_direction_dynamic(const size_t &i) const
  {
    return _dynamic_directions.count(i) > 0;
  }

  /**
   * @brief Get the set of the bundle dynamic direction indices
   *
   * @return the set of bundle dynamic direction indices
   */
  inline const std::set<size_t> &dynamic_directions() const
  {
    return _dynamic_directions;
  }

  /**
   * @brief Get the upper bound of the i-th direction
   *
   * @param i is the index of the direction whose
   *        upper bound must be returned
   * @return a reference to the upper bound of the i-th direction
   */
  const double &get_upper_bound(const unsigned int &i) const
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
  const double &get_lower_bound(const unsigned int &i) const
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
  inline bool is_empty() const
  {
    return static_cast<Polytope>(*this).is_empty();
  }

  /**
   * @brief Test whether the bundle interior is empty
   *
   * @return `true` if and only if the bundle interior is empty
   */
  inline bool is_interior_empty() const
  {
    return static_cast<Polytope>(*this).is_interior_empty();
  }

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
   * @param[in] bundle_template is the bundle template
   * @returns the parallelotope of the bundle associated to `bundle_template`
   */
  Parallelotope get_parallelotope(const BundleTemplate &bundle_template) const;

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
   * is moved by a value `delta`.
   *
   * @param delta is the aimed expansion
   * @return a reference to the updated bundle
   */
  Bundle &expand_by(const double delta);

  /**
   * Compute the edge lengths
   *
   * This method compute the edge lengths, i.e., the distances between the
   * half-spaces bounding the bundle.
   *
   * @returns the vector of the edge lengths
   */
  LinearAlgebra::Vector<double> edge_lengths() const;

  virtual ~Bundle();

  friend void swap(Bundle &A, Bundle &B);

  friend Bundle over_approximate_union(const Bundle &b1, const Bundle &b2);

  friend SetsUnion<Bundle> subtract_and_close(const Bundle &b1,
                                              const Bundle &b2);

  template<typename T>
  friend class Evolver;
  template<typename T>
  friend class CachedEvolver;
};

/**
 * @brief Test whether two bundles are disjoint
 *
 * @param A is a bundle
 * @param B is a bundle
 * @return `true` if and only if `A` and `B` are disjoint
 */
bool are_disjoint(const Bundle &A, const Bundle &B);

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

/**
 * @brief Subtract two bundle and close the result
 *
 * @param b1 is a bundle
 * @param b2 is a bundle
 * @return a union of bundles obtained by closing the set
 *         \f$b1\setminus b2\f$
 */
SetsUnion<Bundle> subtract_and_close(const Bundle &b1, const Bundle &b2);

#endif /* BUNDLE_H_ */
