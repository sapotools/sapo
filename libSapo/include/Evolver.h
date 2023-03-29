/**
 * @file Evolver.h
 * @author Alberto Casagrande <acasagrande@units.it>
 * @brief Handle evolution of geometric sets subject to a discrete system
 * @version 0.1
 * @date 2022-07-02
 *
 * @copyright Copyright (c) 2022
 */

#ifndef EVOLVER_H_
#define EVOLVER_H_

#ifdef WITH_THREADS
#include <shared_mutex>
#endif // WITH_THREADS

#include "DiscreteSystem.h"

#include "Bundle.h"
#include "Polytope.h"
#include "SetsUnion.h"

#include "STL/Atom.h"

// define the maximum admissible length for a parallelotope edge
#define EDGE_MAX_LENGTH 1e18

template<typename T>
class BernsteinCache
{
  /**
   * @brief Direction type
   */
  using direction_type = LinearAlgebra::Vector<T>;

  /**
   * @brief Symbolic Bernstein coefficient vector type
   */
  using coefficients_type = std::vector<SymbolicAlgebra::Expression<T>>;

  /**
   * @brief Parallelotope generators type
   */
  using generators_type = std::vector<LinearAlgebra::Vector<T>>;

  /**
   * @brief Maps that associate directions to Bernstein coefficients
   */
  using direction2coefficient_type
      = std::map<direction_type, coefficients_type>;

  /**
   * @brief Maps that associate generators to direction to Bernstein
   * coefficient map
   */
  using cache_type = std::map<generators_type, direction2coefficient_type>;

  cache_type _cache;

#ifdef WITH_THREADS

  mutable std::shared_timed_mutex _mutex; //!< Cache mutex

#endif // WITH_THREADS

public:
  /**
   * @brief The empty constructor
   */
  BernsteinCache(): _cache() {}

  /**
   * @brief The copy constructor
   *
   * @param orig is the original instance of the object
   */
  BernsteinCache(const BernsteinCache &orig): _cache(orig._cache) {}

  /**
   * @brief Check whether the Bernstein coefficients are cached
   *
   * This method checks whether the symbolic Bernstein coefficients
   * for a parallelotope and a direction are stored in the cache.
   *
   * @param P is a parallelotope
   * @param direction is a direction
   * @return `true` if and only if the coefficients for `P` and
   *         `direction` are cached
   */
  bool coefficients_in_cache(const Parallelotope<T> &P,
                             const LinearAlgebra::Vector<T> &direction) const
  {
#ifdef WITH_THREADS
    std::shared_lock<std::shared_timed_mutex> readlock(_mutex);
#endif // WITH_THREADS

    auto P_found = _cache.find(P.generators());

    if (P_found == std::end(_cache)) {
      return false;
    }

    auto dir_found = P_found->second.find(direction);

    return dir_found != std::end(P_found->second);
  }

  /**
   * @brief Get the cached Bernstein coefficients
   *
   * This method returns a reference to the cached symbolic
   * Bernstein coefficients for a parallelotope and a
   * direction.
   *
   * @param P is a parallelotope
   * @param direction is a direction
   * @return the cached symbolic Bernstein coefficients for `P`,
   *         `direction`, and the current evolver
   */
  const std::vector<SymbolicAlgebra::Expression<T>> &
  get_coefficients(const Parallelotope<T> &P,
                   const LinearAlgebra::Vector<T> &direction) const
  {
#ifdef WITH_THREADS
    std::shared_lock<std::shared_timed_mutex> readlock(_mutex);
#endif // WITH_THREADS

    return _cache.at(P.generators()).at(direction);
  }

public:
  /**
   * @brief Store the Bernstein coefficients
   *
   * This method stores the Bernstein coefficients
   * of a parallelotope and a direction in the cache.
   *
   * @param[in] P is a parallelotope
   * @param[in] direction is a direction
   * @param[in] coefficients is the vector of Bernstein coefficients
   */
  const std::vector<SymbolicAlgebra::Expression<T>> &save_coefficients(
      const Parallelotope<T> &P, const LinearAlgebra::Vector<T> &direction,
      const std::vector<SymbolicAlgebra::Expression<T>> &coefficients)
  {
#ifdef WITH_THREADS
    std::unique_lock<std::shared_timed_mutex> writelock(_mutex);
#endif // WITH_THREADS

    auto &cached_coefficients = _cache[P.generators()][direction];

    cached_coefficients = coefficients;

    return cached_coefficients;
  }

  /**
   * @brief Store the Bernstein coefficients
   *
   * This method stores the Bernstein coefficients
   * of a parallelotope and a direction in the cache.
   *
   * @param[in] P is a parallelotope
   * @param[in] direction is a direction
   * @param[in] coefficients is the vector of Bernstein coefficients
   */
  const std::vector<SymbolicAlgebra::Expression<T>> &
  save_coefficients(const Parallelotope<T> &P,
                    const LinearAlgebra::Vector<T> &direction,
                    std::vector<SymbolicAlgebra::Expression<T>> &&coefficients)
  {
#ifdef WITH_THREADS
    std::unique_lock<std::shared_timed_mutex> writelock(_mutex);
#endif // WITH_THREADS

    auto &cached_coefficients = _cache[P.generators()][direction];

    cached_coefficients = std::move(coefficients);

    return cached_coefficients;
  }
};

/**
 * @brief Handle evolution of geometric sets subject to a discrete system
 *
 * The objects of this class let geometric objects, such as
 * `Parallelotope` or `Bundle`, evolve according to a discrete system.
 * This class computes the Bernstein coefficients of the associated
 * discrete system and the evolving geometric set and use them
 * to over-approximate the set reachable from the source.
 *
 * @tparam T is the type of expression value domain
 */
template<typename T>
class Evolver
{
protected:
  DiscreteSystem<T> _ds; //!< the dynamic system

  BernsteinCache<T> *_cache; //!< the symbolic Bernstein coefficient cache
public:
  /**
   * @brief Approach to evaluate the image of a bundle
   *
   * The are two different approaches to evaluate the
   * image of a bundle through a polynomial function:
   * 1. One-For-One: the boundaries of any parallelotope
   *                 image in the bundle are evaluated by
   *                 exclusively considering the original
   *                 parallelotope itself
   * 2. All-For-One: the boundaries of any parallelotope
   *                 image in the bundle are evaluated by
   *                 exploiting all the bundle templates
   */
  typedef enum {
    ONE_FOR_ONE, /* One-For-One */
    ALL_FOR_ONE  /* All-For-One */
  } evolver_mode;

  evolver_mode mode; //!< the mode used to compute the evolution

  /**
   * @brief A constructor
   *
   * @param discrete_system is the investigated discrete system
   * @param mode is the mode used to compute the evolution
   */
  Evolver(const DiscreteSystem<T> &discrete_system,
          const bool cache_Bernstein_coefficients = true,
          const evolver_mode mode = ALL_FOR_ONE):
      _ds(discrete_system),
      _cache(nullptr), mode(mode)
  {
    if (cache_Bernstein_coefficients) {
      _cache = new BernsteinCache<T>();
    }
  }

  /**
   * @brief A constructor
   *
   * @param discrete_system is the investigated discrete system
   * @param mode is the mode used to compute the evolution
   */
  Evolver(const DiscreteSystem<T> &&discrete_system,
          const bool cache_Bernstein_coefficients = true,
          const evolver_mode mode = ALL_FOR_ONE):
      _ds(std::move(discrete_system)),
      _cache(nullptr), mode(mode)
  {
    if (cache_Bernstein_coefficients) {
      _cache = new BernsteinCache<T>();
    }
  }

  /**
   * @brief Return the dynamical system
   *
   * @return the dynamical system associated to the evolution
   */
  inline const DynamicalSystem<T> &dynamical_system() const
  {
    return _ds;
  }

  /**
   * @brief Transform a bundle according with the system dynamics
   *
   * @param bundle is the bundle to evolve
   * @return an over-approximation of the bundle transformed
   *         by the dynamic laws
   */
  inline Bundle<double> operator()(const Bundle<double> &bundle)
  {
    if (_ds.parameters().size() != 0) {
      SAPO_ERROR("the parameter set has not been specified",
                 std::domain_error);
    }

    return this->operator()(bundle, Polytope<double>());
  }

  /**
   * @brief Compute the image of a bundle
   *
   * This method over-approximates the image of a
   * bundle according to the dynamical system associated
   * to the current evolver.
   *
   * @param bundle is the bundle to evolve
   * @param parameter_set is the parameter set
   * @return an over-approximation of the bundle transformed
   *         by the dynamic laws
   */
  Bundle<double> operator()(const Bundle<double> &bundle,
                            const Polytope<double> &parameter_set);

  /**
   * @brief Transform a bundles union according with the system dynamics
   *
   * @param bundles_union is the set of bundles to evolve
   * @return an over-approximation of set reached from `bundles_union`
   *         by the dynamic laws
   */
  SetsUnion<Bundle<double>>
  operator()(const SetsUnion<Bundle<double>> &bundles_union)
  {
    if (_ds.parameters().size() != 0) {
      SAPO_ERROR("the parameter set has not been specified",
                 std::domain_error);
    }

    SetsUnion<Bundle<double>> result;

    for (auto b_it = std::begin(bundles_union);
         b_it != std::end(bundles_union); ++b_it) {

      result.add(this->operator()(*b_it));
    }

    return result;
  }

  /**
   * @brief Transform a bundles union according with the system dynamics
   *
   * @param bundles_union is the set of bundles to evolve
   * @param parameter_set is the parameter set
   * @return an over-approximation of set reached from `bundles_union`
   *         by the dynamic laws
   */
  SetsUnion<Bundle<double>>
  operator()(const SetsUnion<Bundle<double>> &bundles_union,
             const Polytope<double> &parameter_set)
  {
    SetsUnion<Bundle<double>> result;

    for (auto it = std::begin(bundles_union); it != std::end(bundles_union);
         ++it) {
      result.add(this->operator()(*it, parameter_set));
    }

    return result;
  }

  /**
   * @brief Transform a bundles union according with the system dynamics
   *
   * @param bundles_union is the set of bundles to evolve
   * @param parameter_set is the parameter set
   * @return an over-approximation of set reached from `bundles_union`
   *         by the dynamic laws
   */
  SetsUnion<Bundle<double>>
  operator()(const SetsUnion<Bundle<double>> &bundles_union,
             const SetsUnion<Polytope<double>> &parameter_set)
  {
    SetsUnion<Bundle<double>> result;

    for (auto p_it = std::begin(parameter_set);
         p_it != std::end(parameter_set); ++p_it) {
      for (auto b_it = std::begin(bundles_union);
           b_it != std::end(bundles_union); ++b_it) {
        result.add(this->operator()(*b_it, *p_it));
      }
    }

    return result;
  }

  /**
   * @brief Destroyer
   */
  ~Evolver()
  {
    if (_cache != nullptr) {
      delete _cache;
    }
  }
};

/**
 * @brief Perform the parametric synthesis for an atom
 *
 * This method computes a set of parameters such that the
 * evolution of a bundle satisfies the provided atom.
 *
 * @param[in] evolver is the discrete system evolver
 * @param[in] bundle is the set from which the evolution occurs
 * @param[in] parameter_set is the initial parameter set
 * @param[in] atom is the specification to be satisfied
 * @return a subset of `parameter_set` such that the
 *         evolution of a `bundle` through the dynamical
 *         system when the parameters are in the returned set
 *         satisfies the atomic formula
 */
SetsUnion<Polytope<double>>
synthesize(const Evolver<double> &evolver, const Bundle<double> &bundle,
           const SetsUnion<Polytope<double>> &parameter_set,
           const std::shared_ptr<STL::Atom> atom);

#endif // EVOLVER_H_
