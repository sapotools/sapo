/**
 * @file Flowpipe.cpp
 * Represent and manipulate flowpipes of bundle
 * Used to represent the reachable set of a dynamical system
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#include <fstream>
#include <string>

#include "Flowpipe.h"

/**
 * Constructor that instantiates Flowpipe
 */
Flowpipe::Flowpipe(): v_templates(), flowpipe() {}

/**
 * Constructor that instantiates Flowpipe
 */
Flowpipe::Flowpipe(const std::vector<std::vector<double>> &variable_templates):
    v_templates(variable_templates), flowpipe()
{
}

/**
 * Return the i-th polytopes union
 *
 * @param[in] i index
 * @return i-th polytopes union
 */
const PolytopesUnion &Flowpipe::get(const unsigned int i) const
{
  if (i >= this->size()) {
    std::domain_error("Flowpipe::get: i must be between 0 and "
                      "the flowpipe size");
  }

  return this->flowpipe[i];
}

/**
 * Append a bundle to the flowpipe
 *
 * @param[in] bundle bundle to be appended
 * @return a reference to the new flowpipe
 */
Flowpipe &Flowpipe::append(const Bundle &bundle)
{
  this->flowpipe.push_back(PolytopesUnion(bundle));

  return *this;
}

/**
 * Append a polytopes union to the flowpipe
 *
 * @param[in] P is the polytopes union to be appended
 * @return a reference to the new flowpipe
 */
Flowpipe &Flowpipe::append(const PolytopesUnion &Ps)
{
  this->flowpipe.push_back(Ps);

  return *this;
}

/**
 * Append a polytope to the flowpipe
 *
 * @param[in] P is the polytopes union to be appended
 * @return a reference to the new flowpipe
 */
Flowpipe &Flowpipe::append(const Polytope &P)
{
  this->flowpipe.push_back(PolytopesUnion(P));

  return *this;
}

unsigned int Flowpipe::dim() const
{
  if (this->flowpipe.empty()) {
    return 0;
  }

  return (this->flowpipe[0]).dim();
}

/**
 * Print the polytope in Matlab format (for plotregion script)
 *
 * @param[in] os is the output stream
 * @param[in] color color of the polytope to plot
 */
void Flowpipe::plotRegion(std::ostream &os, const char color) const
{
  for (auto it = std::begin(flowpipe); it != std::end(flowpipe); ++it) {
    it->plotRegion(os, color);
  }
}

/**
 * Print the projection of the variable in time using Matlab format into a file
 *
 * @param[in] os is the output stream
 * @param[in] var variable to be plotted
 * @param[in] time_step time step for time projection
 * @param[in] color color of the polytope to plot
 */
void Flowpipe::plotProj(std::ostream &os, const unsigned int var,
                        const double time_step, const char color) const
{

  if (size() == 0 || this->flowpipe[0].is_empty()) {
    std::domain_error("Flowpipe::plotProj: The flowpipe must be "
                      "non-empty");
  }

  if (var >= this->flowpipe[0].dim()) {
    std::domain_error("Flowpipe::plotProj: var must be between "
                      "0 and the system dimension");
  }

  // select figure
  os << "figure(" << var + 1 << ")" << std::endl;

  // print time
  os << "t = [ ";
  for (unsigned int i = 0; i < this->size(); i++) {
    os << i * time_step << " ";
  }
  os << " ];" << std::endl;

  // print lower offsets
  os << "varm = [";
  for (auto it = std::begin(this->flowpipe); it != std::end(this->flowpipe);
       ++it) {
    PolytopesUnion::const_iterator ls_it(it->begin());

    double min_value = ls_it->minimize(v_templates[var]).optimum();

    for (; ls_it != it->end(); ++ls_it) {
      double min_var_value = ls_it->minimize(v_templates[var]).optimum();

      min_value = std::min(min_var_value, min_value);
    }
    os << " " << min_value;
  }
  os << " ];" << std::endl;

  // print upper offsets
  os << "varp = [";
  for (auto it = std::begin(this->flowpipe); it != std::end(this->flowpipe);
       ++it) {
    std::vector<double> obj_funct(dim(), 0.0);
    obj_funct[var] = 1.0;

    PolytopesUnion::const_iterator ls_it(it->begin());

    double max_value = ls_it->maximize(v_templates[var]).optimum();

    for (; ls_it != it->end(); ++ls_it) {
      double max_var_value = ls_it->maximize(v_templates[var]).optimum();

      max_value = std::min(max_var_value, max_value);
    }
    os << " " << max_value;
  }
  os << " ];" << std::endl;

  os << "T = [t,fliplr(t)];" << std::endl;
  os << "X = [varm,fliplr(varp)];" << std::endl;
  os << "fill(T,X,'" << color << "');" << std::endl;
}
