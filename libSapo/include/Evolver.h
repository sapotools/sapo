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
protected:
  DynamicalSystem<T> _ds; //!< The dynamical system

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
  Evolver(DynamicalSystem<T> &&dynamical_system,
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
  inline Bundle operator()(const Bundle &bundle) const
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
  Bundle operator()(const Bundle &bundle, const Polytope &parameter_set) const;

  /**
   * @brief Transform a bundles union according with the system dynamics
   *
   * @param bundles_union is the set of bundles to evolve
   * @return an over-approximation of set reached from `bundles_union`
   *         by the dynamic laws
   */
  SetsUnion<Bundle> operator()(const SetsUnion<Bundle> &bundles_union) const
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
                               const Polytope &parameter_set) const
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
                               const SetsUnion<Polytope> &parameter_set) const
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

#endif // EVOLVER_H_
