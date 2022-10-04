#ifndef __INPUTDATA_H__
#define __INPUTDATA_H__

#include <vector>
#include <algorithm>

#include <STL/STL.h>
#include <STL/Conjunction.h>

#include <LinearSystem.h>
#include <LinearAlgebra.h>

#include "AbsSynIO.h"

#include "Sapo.h"
#include "types.h"
#include "Variable.h"
#include "Parameter.h"
#include "Constant.h"
#include "Definition.h"
#include "Direction.h"

namespace AbsSyn
{

class InputData
{
  friend std::ostream &operator<<(std::ostream &os, const InputData &id);

public:
  InputData();

  ~InputData();

  bool isProblemDefined() const
  {
    return problem != problemType::P_UNDEF;
  }
  const problemType &getProblem() const
  {
    return problem;
  }
  void setProblem(problemType t)
  {
    problem = t;
  }

  bool isVarModeDefined() const
  {
    return varMode != modeType::M_UNDEF;
  }
  const modeType &getVarMode() const
  {
    return varMode;
  }
  void setVarMode(modeType t)
  {
    varMode = t;
  }

  bool isParamModeDefined() const
  {
    return paramMode != modeType::M_UNDEF;
  }
  const modeType &getParamMode() const
  {
    return paramMode;
  }
  void setParamMode(modeType t)
  {
    paramMode = t;
  }

  bool isIterationSet() const
  {
    return iter_set;
  }

  const unsigned int &getIterations() const
  {
    return iterations;
  }

  void setIterations(unsigned int n)
  {
    iter_set = true;
    iterations = n;
  }

  const unsigned int &getMaxKInduction() const
  {
    return max_k_induction;
  }

  void setMaxKInduction(unsigned int n)
  {
    max_k_induction = n;
  }

  const double &getMaxVersorMagnitude() const
  {
    return max_bundle_magnitude;
  }

  const unsigned int &getMaxParameterSplits() const
  {
    return max_param_splits;
  }

  void setMaxParameterSplits(unsigned int n)
  {
    max_param_splits = n;
  }

  void setMaxVersorMagnitude(const double magnitude)
  {
    max_bundle_magnitude = magnitude;
  }

  const bool &isPreSplitsSet() const
  {
    return presplits;
  }

  void setPreSplits(const bool presplits)
  {
    this->presplits = presplits;
  }

  unsigned getVarNum() const
  {
    return vars.size();
  }
  unsigned getParamNum() const
  {
    return params.size();
  }
  unsigned getConstNum() const
  {
    return consts.size();
  }
  unsigned getDefNumber() const
  {
    return defs.size();
  }

  bool isVarDefined(const std::string &name)
      const; // checks if a variable named 'name' already exists
  bool isParamDefined(const std::string &name)
      const; // checks if a parameter named 'name' already exists
  bool isConstDefined(const std::string &name)
      const; // checks if a constant named 'name' already exists
  bool isDefDefined(const std::string &name)
      const; // checks if a definition named 'name' already exists
  bool isDirectionDefined(const std::string &name)
      const; // checks if a direction named "name" already exists
  bool isSymbolDefined(
      const std::string &name) const; // checks if a symbol (var, param, const
                                      // or def) named 'name' already exists

  SymbolicAlgebra::Symbol<> getSymbol(
      std::string name) const; // returns the symbol named "name" (must exist)

  const Variable *getVar(int i) const
  {
    return vars[i];
  }
  const Variable *getVar(const std::string &name)
      const; // return the variable named 'name', which must exist
  Variable *getVar(int i)
  {
    return vars[i];
  }
  Variable *getVar(const std::string &name);
  int getVarPos(const std::string &name)
      const; // return an index such that vars[i] has name 'name'
  std::vector<SymbolicAlgebra::Symbol<>> getVarSymbols() const;

  const Parameter *getParam(int i) const
  {
    return params[i];
  }
  const Parameter *getParam(const std::string &name)
      const; // return the parameter named 'name', which must exist
  int getParamPos(const std::string &name)
      const; // return an index i such that params[i] has name 'name'
  std::vector<SymbolicAlgebra::Symbol<>> getParamSymbols() const;

  const Constant *getConst(int i) const
  {
    return consts[i];
  }
  const Constant *getConst(const std::string &name)
      const; // return the constant named 'name', which must exist
  const Definition *getDef(int i) const
  {
    return defs[i];
  }
  const Definition *getDef(const std::string &name)
      const; // return the definition named 'name', which must exist
  int getDefPos(const std::string &name)
      const; // return an index i such that defs[i] has name 'name'

  const std::list<Direction *> &getAssumptions() const
  {
    return assumptions;
  }

  void addVariable(
      Variable *v) // adds a new variable, which name is not already used
  {
    vars.push_back(v);
  }
  void addParameter(
      Parameter *p) // adds a new parameter, which name is not already used
  {
    params.push_back(p);
  }
  void addConstant(Constant *c)
  {
    consts.push_back(c);
  } // adds a new constant, which name is not already used
  void addDefinition(Definition *d)
  {
    defs.push_back(d);
  } // adds a new definition, which name is not already used
  void addAssumption(Direction *d)
  {
    assumptions.push_back(d);
  } // adds a new assumption

  bool isSpecDefined() const
  {
    return spec != NULL;
  }

  void addSpec(std::shared_ptr<STL::STL> f)
  {
    if (spec == NULL) {
      spec = f;
    } else {
      spec = std::make_shared<STL::Conjunction>(spec, f);
    }
  }

  const std::shared_ptr<STL::STL> specification() const
  {
    return spec;
  }

  // add direction with specified name
  void addVarDirectionConstraint(Direction *d);

  unsigned int getDirectionsNum() const
  {
    return directions.size();
  }
  const std::vector<Direction *> &getDirections() const
  {
    return directions;
  }
  const Direction *getDirection(const unsigned int &i) const
  {
    return directions[i];
  }
  bool isBounded(const unsigned int &d) const
  {
    return directions[d]->hasLB() && directions[d]->hasUB();
  }

  void addInvariantConstraint(Direction *d);

  const std::list<Direction *> &getInvariant() const
  {
    return invariant;
  }

  inline void setApproxType(const Sapo::joinApproxType type)
  {
    approx_type = type;
  }

  inline const Sapo::joinApproxType &getApproxType() const
  {
    return approx_type;
  }

  unsigned int findDirectionPos(const std::string &name) const;

  unsigned templateRows() const
  {
    return templateMatrix.size();
  }
  unsigned templateCols() const
  {
    return templateMatrix[0].size();
  }
  void setTemplate(const std::vector<std::vector<unsigned int>> &m)
  {
    templateMatrix = m;
  }
  const std::vector<std::vector<unsigned int>> &getTemplate() const
  {
    return templateMatrix;
  }

  void addParamDirectionConstraint(Direction *d);
  unsigned paramDirectionsNum() const
  {
    return paramDirections.size();
  }

  const std::vector<Direction *> &getParameterDirections() const
  {
    return paramDirections;
  }

  Direction *getParamDirection(int i) const
  {
    return paramDirections[i];
  }

  bool isTransModeDefined() const
  {
    return trans != transType::T_UNDEF;
  }

  void setTransMode(transType t)
  {
    trans = t;
  }

  const transType &getTransMode() const
  {
    return trans;
  }

  transType getTransValue() const
  {
    return trans;
  }

  const bool &isDecompositionDefined() const
  {
    return decomp_defined;
  }
  void setDecomposition()
  {
    decomp = true;
  }
  const bool &getDecomposition() const
  {
    return decomp;
  }

  bool isAlphaDefined() const
  {
    return alphaDefined;
  }
  void setAlpha(double a)
  {
    alpha = a;
  }
  const double &getAlpha() const
  {
    return alpha;
  }

  bool isDynamicCompositionEnabled() const
  {
    return compose_dynamic;
  }

  unsigned getDynamicDegree() const
  {
    return dynamic_degree;
  }

  void setDynamicDegree(unsigned d)
  {
    compose_dynamic = true;
    dynamic_degree = d;
  }

  void setBernsteinCaching(bool flag)
  {
    bern_caching = flag;
  }

  bool useBernsteinCaching() const
  {
    return bern_caching;
  }

  void setUseInvariantDirections(bool flag)
  {
    use_invariant_directions = flag;
  }

  bool getUseInvariantDirections() const
  {
    return use_invariant_directions;
  }

  /**
   * @brief Fix boundaries according to the whole input
   *
   * There are cases in which a set of constraints is bounded
   * even though not all constrains are explicitly provided
   * with an upper and a lower bound. For instance, neither
   * \f$x\f$ nor \f$y\f$ are apparently upper bounded by the
   * system \f$\{x + y <= 10, x >= 0, y >= 0\}\f$, however,
   * the constrains \f$x + y <= 10\f$ forces both of them
   * to be $10$ at most.
   * This method considers both the initial set and the
   * parameter set, discovers the upper and the lower
   * boundaries over each constraint, and updates their
   * values in the current object.
   */
  void optimize_boundaries();

  bool check(); // checks for errors in model

protected:
  problemType problem;

  modeType varMode;
  modeType paramMode;

  unsigned int iterations;
  bool iter_set;

  unsigned int max_k_induction;

  unsigned int max_param_splits;
  bool presplits;

  double max_bundle_magnitude;

  std::vector<Variable *> vars;
  std::vector<Parameter *> params;
  std::vector<Constant *> consts;
  std::vector<Definition *> defs;

  std::list<Direction *> assumptions;
  std::list<Direction *> invariant;

  std::shared_ptr<STL::STL> spec;

  std::vector<Direction *> directions;

  std::vector<std::vector<unsigned int>> templateMatrix;

  std::vector<Direction *> paramDirections;

  // SAPO options
  transType trans;
  bool decomp, decomp_defined;
  double alpha;
  bool alphaDefined;
  bool compose_dynamic;
  unsigned dynamic_degree;
  Sapo::joinApproxType approx_type;
  bool bern_caching;
  bool use_invariant_directions;

  // addition of direction to params or vars
  void addDirectionConstraint(Direction *d, bool isVar);
};

/**
 * @brief Build a constraints linear system from a set of constraints
 *
 * @tparam T is the type of the symbol domain
 * @param constraints the set of constraints
 * @param symbols the order of the constraint symbols in the linear system
 * @return A linear system representing the provided linear systems
 */
template<typename T, template<class, class> class CONTAINER>
LinearSystem getConstraintsSystem(
    const CONTAINER<Direction *, std::allocator<Direction *>> &constraints,
    const std::vector<SymbolicAlgebra::Symbol<T>> &symbols)
{
  using namespace LinearAlgebra;

  std::vector<std::vector<T>> A;
  std::vector<T> b;

  for (auto dir_it = std::cbegin(constraints);
       dir_it != std::cend(constraints); ++dir_it) {
    const AbsSyn::Direction &dir = **dir_it;

    auto systemRow = dir.getConstraintVector(symbols);
    // add only bounded constraints
    if (dir.hasUB()) {
      A.push_back(systemRow);
      b.push_back(dir.getUB());
    }
    if (dir.hasLB()) {
      A.push_back(-systemRow);
      b.push_back(-dir.getLB());
    }
  }
  return LinearSystem(A, b);
}

}

#endif
