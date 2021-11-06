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
#include "LinearSystemSet.h"

class Flowpipe
{

private:
  std::vector<std::vector<double>> v_templates;
  std::vector<LinearSystemSet> flowpipe; // flowpipe

public:
  /**
   * An empty flowpipe constructor
   */
  Flowpipe();

  /**
   * An empty flowpipe with variable templates constructor
   */
  Flowpipe(const std::vector<std::vector<double>> &variable_templates);

  const LinearSystemSet &
  get(const unsigned int i) const; // get i-th LinearSystemSet

  /**
   * Append a linear system to the flowpipe
   *
   * @param[in] ls is the linear system set to be appended
   * @return a reference to the new flowpipe
   */
  Flowpipe &append(const LinearSystem &ls);

  /**
   * Append a linear system set to the flowpipe
   *
   * @param[in] ls is the linear system set to be appended
   * @return a reference to the new flowpipe
   */
  Flowpipe &append(const LinearSystemSet &ls);

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
  friend std::ostream &operator<<(std::ostream &out, const Flowpipe &fp);

  ~Flowpipe()
  {}
};

#endif /* BUNDLE_H_ */
