/**
 * @file Model.h
 * Represent a discrete-time (eventually parameteric) dynamical system
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#ifndef MODEL_H_
#define MODEL_H_

#include <ginac/ginac.h>

#include "Bundle.h"
#include "LinearSystemSet.h"
#include "STL.h"

class Model
{

protected:
  GiNaC::lst vars;   // variables
  GiNaC::lst params; // parameters
  GiNaC::lst dyns;   // dynamics

  Bundle *reachSet; // Initial reach set
  LinearSystemSet *paraSet;

  std::shared_ptr<STL> spec;

  std::string name;

public:
  const std::string &getName() const
  {
    return this->name;
  }

  const GiNaC::lst &getVars() const
  {
    return this->vars;
  }

  const GiNaC::lst &getParams() const
  {
    return this->params;
  }

  const GiNaC::lst &getDyns() const
  {
    return this->dyns;
  }

  const Bundle *getReachSet() const
  {
    return this->reachSet;
  }
  const LinearSystemSet *getParaSet() const
  {
    return this->paraSet;
  }
  const std::shared_ptr<STL> getSpec() const
  {
    return this->spec;
  }

  ~Model();
};

#endif /* MODEL_H_ */
