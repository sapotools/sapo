/**
 * @file Flowpipe.h
 * Represent and manipulate flowpipes of bundle
 * Used to represent the reachable set of a dynamical system
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#ifndef FLOWPIPE_H_
#define FLOWPIPE_H_

#include <vector>

#include "Bundle.h"
#include "PolytopesUnion.h"
#include "OutputFormatter.h"

/**
 * @brief A representation for reachability flowpipe
 * 
 * This class represents reachability flowpipe as 
 * a sequence of polytope union.
 */
class Flowpipe: public std::vector<PolytopesUnion>
{

public:
  /**
   * An empty flowpipe constructor
   */
  Flowpipe();

  const PolytopesUnion &
  get(const unsigned int i) const; // get i-th PolytopesUnion

  /**
   * Append a polytopes union to the flowpipe
   *
   * @param[in] Pu is the polytopes union to be appended
   * @return a reference to the new flowpipe
   */
  inline Flowpipe &push_back(const PolytopesUnion &Pu)
  {
    std::vector<PolytopesUnion>::push_back(Pu);

    return *this;
  }

  /**
   * Append a polytopes union to the flowpipe
   *
   * @param[in] Pu is the polytopes union to be appended
   * @return a reference to the new flowpipe
   */
  inline Flowpipe &push_back(PolytopesUnion &&Pu)
  {
    std::vector<PolytopesUnion>::push_back(std::move(Pu));

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

  /**
   * Stream a flowpipe
   *
   * @param[in] os is the output stream
   * @param[in] fs is the flowpipe to be streamed
   * @return the output stream
   */
  template<typename OSTREAM>
  friend OSTREAM &operator<<(OSTREAM &os, const Flowpipe &fp)
  {
    using OF = OutputFormatter<OSTREAM>;

    os << OF::sequence_begin();
    for (auto it = std::cbegin(fp); it != std::cend(fp); ++it) {
      if (it != std::cbegin(fp)) {
        os << OF::sequence_separator();
      }
      os << *it;
    }
    os << OF::sequence_end();

    return os;
  }
};

#endif /* BUNDLE_H_ */
