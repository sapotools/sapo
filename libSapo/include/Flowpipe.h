/**
 * @file Flowpipe.h
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @author Alberto Casagrande <acasagrande@units.it>
 * @brief Represent and manipulate reachability flowpipes
 * @version 0.2
 * @date 2022-05-04
 *
 * @copyright Copyright (c) 2016-2022
 */

#ifndef FLOWPIPE_H_
#define FLOWPIPE_H_

#include <vector>

#include "Bundle.h"
#include "SetsUnion.h"

/**
 * @brief A representation for reachability flowpipe
 *
 * This class represents reachability flowpipe as
 * a sequence of polytope union.
 * @todo reimplement this class as a list of `SetsUnion`
 */
class Flowpipe : public std::vector<SetsUnion<Polytope>>
{

public:
  /**
   * An empty flowpipe constructor
   */
  Flowpipe();

  const SetsUnion<Polytope> &
  get(const unsigned int i) const; // get i-th polytopes union

  /**
   * Append a polytopes union to the flowpipe
   *
   * @param[in] Pu is the polytopes union to be appended
   * @return a reference to the new flowpipe
   */
  inline Flowpipe &push_back(const SetsUnion<Polytope> &Pu)
  {
    std::vector<SetsUnion<Polytope>>::push_back(Pu);

    return *this;
  }

  /**
   * Append a polytopes union to the flowpipe
   *
   * @param[in] Pu is the polytopes union to be appended
   * @return a reference to the new flowpipe
   */
  inline Flowpipe &push_back(SetsUnion<Polytope> &&Pu)
  {
    std::vector<SetsUnion<Polytope>>::push_back(std::move(Pu));

    return *this;
  }

  /**
   * Append a bundles union to the flowpipe
   *
   * @param[in] Bu is the bundles union to be appended
   * @return a reference to the new flowpipe
   */
  inline Flowpipe &push_back(const SetsUnion<Bundle> &Bu)
  {
    SetsUnion<Polytope> Pu;

    for (const Bundle &b: Bu) {
      Pu.add(b);
    }

    push_back(Pu);

    return *this;
  }

  /**
   * Append a bundle to the flowpipe
   *
   * @param[in] bundle bundle to be appended
   * @return a reference to the new flowpipe
   */
  inline Flowpipe &push_back(const Bundle &bundle)
  {
    this->emplace(this->end(), bundle);

    return *this;
  }

  /**
   * Get the number of variables
   *
   * @returns number of variables stored in the flowpipe
   */
  unsigned int dim() const;
};

#endif /* FLOWPIPE_H_ */
