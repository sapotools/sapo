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

class Flowpipe
{

private:
  std::vector<std::vector<double>> v_templates;
  std::vector<PolytopesUnion> flowpipe; // flowpipe

public:
  /**
   * An empty flowpipe constructor
   */
  Flowpipe();

  /**
   * An empty flowpipe with variable templates constructor
   */
  Flowpipe(const std::vector<std::vector<double>> &variable_templates);

  const PolytopesUnion &
  get(const unsigned int i) const; // get i-th PolytopesUnion

  /**
   * Append a polytope to the flowpipe
   *
   * @param[in] P is the polytopes union to be appended
   * @return a reference to the new flowpipe
   */
  Flowpipe &append(const Polytope &P);

  /**
   * Append a polytopes union to the flowpipe
   *
   * @param[in] Pu is the polytopes union to be appended
   * @return a reference to the new flowpipe
   */
  Flowpipe &append(const PolytopesUnion &Pu);

  /**
   * Append a bundle to the flowpipe
   *
   * @param[in] bundle bundle to be appended
   * @return a reference to the new flowpipe
   */
  Flowpipe &append(const Bundle &bundle);

  std::size_t size() const
  {
    return this->flowpipe.size();
  }

  /**
   * Get the number of variables
   *
   * @returns number of variables stored in the flowpipe
   */
  unsigned int dim() const;

  /**
   * Print the flowpipe in Matlab format (for plotregion script)
   *
   * @param[in] os is the output stream
   * @param[in] color color of the polytope to plot
   */
  void plotRegion(std::ostream &os = std::cout, const char color = ' ') const;

  void plotProj(std::ostream &os, const unsigned int var,
                const double time_step, const char color) const;

  friend void swap(Flowpipe &A, Flowpipe &B)
  {
    swap(A.flowpipe, B.flowpipe);
    swap(A.v_templates, B.v_templates);
  }

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
    for (auto it = std::cbegin(fp.flowpipe); it != std::cend(fp.flowpipe);
         ++it) {
      if (it != std::cbegin(fp.flowpipe)) {
        os << OF::sequence_separator();
      }
      os << *it;
    }
    os << OF::sequence_end();

    return os;
  }

  ~Flowpipe() {}
};

#endif /* BUNDLE_H_ */
