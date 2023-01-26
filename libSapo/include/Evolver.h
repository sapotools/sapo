/**
 * @file Evolver.h
 * @author Alberto Casagrande <acasagrande@units.it>
 * @brief Handle evolution of geometric sets subject to a dynamical system
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

#include "DynamicalSystem.h"

#include "Bundle.h"
#include "Polytope.h"
#include "SetsUnion.h"

#include "STL/Atom.h"

/**
 * @brief Handle evolution of geometric sets subject to a dynamical system
 *
 * The objects of this class let geometric objects, such as
 * `Parallelotope` or `Bundle`, subject to a dynamical system evolve.
 * This class computes the Bernstein coefficients of the associated
 * dynamical system and the evolving geometric set and use them
 * to over-approximate the set reachable from the source.
 *
 * @tparam T is the type of expression value domain
 */
template<typename T>
class Evolver
{
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

protected:
  DynamicalSystem<T> _ds; //!< The dynamical system

  const T _edge_threshold{static_cast<T>(1e18)};

public:
  evolver_mode mode; //!< the mode used to compute the evolution

  /**
   * @brief A constructor
   *
   * @param dynamical_system is the investigated dynamical system
   * @param mode is the mode used to compute the evolution
   */
  Evolver(const DynamicalSystem<T> &dynamical_system,
          const evolver_mode mode = ALL_FOR_ONE):
      _ds(dynamical_system),
      mode(mode)
  {
  }

  /**
   * @brief A constructor
   *
   * @param dynamical_system is the investigated dynamical system
   * @param mode is the mode used to compute the evolution
   */
  Evolver(const DynamicalSystem<T> &&dynamical_system,
          const evolver_mode mode = ALL_FOR_ONE):
      _ds(std::move(dynamical_system)),
      mode(mode)
  {
  }

  /**
   * @brief A copy constructor
   *
   * @param orig is the template object
   */
  Evolver(const Evolver<T> &orig): _ds(orig._ds), mode(orig.mode) {}

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
  inline Bundle operator()(const Bundle &bundle)
  {
    if (_ds.parameters().size() != 0) {
      throw std::domain_error("The parameter set has not been specified");
    }

    return this->operator()(bundle, Polytope());
  }

  /**
   * @brief Transform a bundle according with the system dynamics
   *
   * @param bundle is the bundle to evolve
   * @param parameter_set is the parameter set
   * @return an over-approximation of the bundle transformed
   *         by the dynamic laws
   */
  virtual Bundle operator()(const Bundle &bundle,
                            const Polytope &parameter_set);

  /**
   * @brief Transform a bundles union according with the system dynamics
   *
   * @param bundles_union is the set of bundles to evolve
   * @return an over-approximation of set reached from `bundles_union`
   *         by the dynamic laws
   */
  SetsUnion<Bundle> operator()(const SetsUnion<Bundle> &bundles_union)
  {
    if (_ds.parameters().size() != 0) {
      throw std::domain_error("The parameter set has not been specified");
    }

    SetsUnion<Bundle> result;

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
  SetsUnion<Bundle> operator()(const SetsUnion<Bundle> &bundles_union,
                               const Polytope &parameter_set)
  {
    SetsUnion<Bundle> result;

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
  SetsUnion<Bundle> operator()(const SetsUnion<Bundle> &bundles_union,
                               const SetsUnion<Polytope> &parameter_set)
  {
    SetsUnion<Bundle> result;

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
  virtual ~Evolver() {}
};

/**
 * @brief Perform the parametric synthesis for an atom
 *
 * This method computes a set of parameters such that the
 * evolution of a bundle satisfies the provided atom.
 *
 * @param[in] evolver is the dynamical system evolver
 * @param[in] bundle is the set from which the evolution occurs
 * @param[in] parameter_set is the initial parameter set
 * @param[in] atom is the specification to be satisfied
 * @return a subset of `parameter_set` such that the
 *         evolution of a `bundle` through the dynamical
 *         system when the parameters are in the returned set
 *         satisfies the atomic formula
 */
SetsUnion<Polytope> synthesize(const Evolver<double> &evolver,
                               const Bundle &bundle,
                               const SetsUnion<Polytope> &parameter_set,
                               const std::shared_ptr<STL::Atom> atom);

/**
 * @brief Handle evolution of geometric sets subject to a dynamical system
 *
 * The objects of this class let geometric objects, such as
 * `Parallelotope` or `Bundle`, subject to a dynamical system evolve.
 * This class computes the symbolic Bernstein coefficients of the
 * associated dynamical system and the evolving geometric set and
 * store them in a cache. When needed the computed symbolic
 * coefficient are recovered and instantiated for the specific
 * geometric set.
 *
 * @tparam T is the type of expression value domain
 */
template<typename T>
class CachedEvolver : public Evolver<T>
{
protected:
  /**
   * @brief Direction type
   */
  typedef LinearAlgebra::Vector<T> dir_type;

  /**
   * @brief Maps that associate directions to Bernstein coefficients
   */
  typedef std::map<dir_type, std::vector<SymbolicAlgebra::Expression<T>>>
      dir2coeffs_type;

  /**
   * @brief Maps that associate generators to direction to Bernstein
   * coefficient map
   */
  typedef std::map<LinearAlgebra::Dense::Matrix<T>, dir2coeffs_type>
      cache_type;

  cache_type _cache; //!< The evolver cache for symbolic Bernstein coefficients

#ifdef WITH_THREADS
  mutable std::shared_timed_mutex _cache_mutex; //!< Cache mutex
#endif                                          // WITH_THREADS

public:
  /**
   * @brief Check whether the Bernstein coefficients are cached
   *
   * @param P is a parallelotope
   * @param direction is a direction
   * @return `true` if and only if the coefficients for `P` and
   *         `direction` are cached
   */
  bool coefficients_in_cache(const Parallelotope &P,
                             const LinearAlgebra::Vector<T> &direction) const
  {
#ifdef WITH_THREADS
    std::shared_lock<std::shared_timed_mutex> readlock(_cache_mutex);
#endif // WITH_THREADS

    auto P_found = _cache.find(P.generators());

    if (P_found == std::end(_cache)) {
      return false;
    }

    auto dir_found = P_found->second.find(direction);

    return dir_found != std::end(P_found->second);
  }

  /**
   * @brief Get the cached Bernstein coefficients of a parallelotope
   *
   * @param P is a parallelotope
   * @param direction is a direction
   * @return the cached symbolic Bernstein coefficients for `P`,
   *         `direction`, and the current evolver
   */
  const std::vector<SymbolicAlgebra::Expression<T>> &
  get_cached_coefficients(const Parallelotope &P,
                          const LinearAlgebra::Vector<T> &direction) const
  {
#ifdef WITH_THREADS
    std::shared_lock<std::shared_timed_mutex> readlock(_cache_mutex);
#endif // WITH_THREADS

    return _cache.at(P.generators()).at(direction);
  }

  /**
   * @brief Cache the Bernstein coefficients
   *
   * This method saves the Bernstein coefficients
   * of a parallelotope, a direction, and the
   * dynamics of the current evaluator in the cache.
   *
   * @param[in] P is a parallelotope
   * @param[in] direction is a direction
   * @param[in] coeffs is the vector of Bernstein coefficients
   */
  const std::vector<SymbolicAlgebra::Expression<T>> &set_cache_coefficients(
      const Parallelotope &P, const LinearAlgebra::Vector<T> &direction,
      const std::vector<SymbolicAlgebra::Expression<T>> &coeffs)
  {
#ifdef WITH_THREADS
    std::unique_lock<std::shared_timed_mutex> writelock(_cache_mutex);
#endif // WITH_THREADS

    auto &cache = _cache[P.generators()][direction];

    cache = coeffs;

    return cache;
  }

  /**
   * @brief Cache the Bernstein coefficients
   *
   * This method saves the Bernstein coefficients
   * of a parallelotope, a direction, and the
   * dynamics of the current evaluator in the cache.
   *
   * @param P is a parallelotope
   * @param direction is a direction
   * @param coeffs is the vector of Bernstein coefficients
   */
  const std::vector<SymbolicAlgebra::Expression<T>> &
  set_cache_coefficients(const Parallelotope &P,
                         const LinearAlgebra::Vector<T> &direction,
                         std::vector<SymbolicAlgebra::Expression<T>> &&coeffs)
  {
#ifdef WITH_THREADS
    std::unique_lock<std::shared_timed_mutex> writelock(_cache_mutex);
#endif // WITH_THREADS

    // auto & cache = _cache[P.generators()][direction];

    _cache[P.generators()][direction] = std::move(coeffs);

    return _cache[P.generators()][direction];
  }

  /**
   * @brief A constructor
   *
   * @param dynamical_system is the investigated dynamical system
   * @param mode is the mode used to compute the evolution
   */
  CachedEvolver(const DynamicalSystem<T> &dynamical_system,
                const typename Evolver<T>::evolver_mode mode
                = Evolver<T>::ALL_FOR_ONE):
      Evolver<T>(dynamical_system, mode),
      _cache()
  {
  }

  /**
   * @brief A constructor
   *
   * @param dynamical_system is the investigated dynamical system
   * @param mode is the mode used to compute the evolution
   */
  CachedEvolver(DynamicalSystem<T> &&dynamical_system,
                const typename Evolver<T>::evolver_mode mode
                = Evolver<T>::ALL_FOR_ONE):
      Evolver<T>(std::move(dynamical_system), mode),
      _cache()
  {
  }

  /**
   * @brief A constructor
   *
   * @param evolver is the original evolver
   */
  CachedEvolver(const Evolver<T> &evolver): Evolver<T>(evolver), _cache() {}

  /**
   * @brief A copy constructor
   *
   * @param orig is the template object
   */
  CachedEvolver(const CachedEvolver<T> &orig):
      Evolver<T>(orig), _cache(orig._cache)
  {
  }

  /**
   * @brief Let a bundle evolve according with the dynamical system
   *
   * @param bundle is the bundle to evolve
   * @param parameter_set is the parameter set
   * @return an over-approximation of the bundle transformed
   *         by the dynamic laws
   */
  Bundle operator()(const Bundle &bundle, const Polytope &parameter_set);
};

#endif // EVOLVER_H_
