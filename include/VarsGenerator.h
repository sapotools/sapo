/**
 * @file VarsGenerator.h
 * Automatically generate variables for paralleltope generator functions.
 * For high dimensions declaring manually the variables can be tedious...
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#ifndef VARSGENERATOR_H_
#define VARSGENERATOR_H_

#include "Common.h"
#include "LinearSystem.h"

class VarsGenerator
{

private:
  unsigned int dim;

  GiNaC::lst qs;              // base vertex
  GiNaC::lst as;              // alphas (free variables)
  GiNaC::lst bs;              // betas (lenghts)
  GiNaC::lst ls;              // directions (lambdas)
  std::vector<GiNaC::lst> us; // genertor versors

public:
  VarsGenerator(const unsigned int dim);

  const GiNaC::lst &getBaseVertex() const
  {
    return this->qs;
  };
  const GiNaC::lst &getFreeVars() const
  {
    return this->as;
  };
  const GiNaC::lst &getLenghts() const
  {
    return this->bs;
  };
  const GiNaC::lst &getDirections() const
  {
    return this->ls;
  };
  const std::vector<GiNaC::lst> &getVersors() const
  {
    return this->us;
  };

  LinearSystem genBox(const std::vector<double> &b) const;

  virtual ~VarsGenerator();
};

#endif /* VARSGENERATOR_H_*/
