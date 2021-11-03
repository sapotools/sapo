/**
 * @file Sapo.cpp
 * Core of Sapo tool.
 * Here the reachable set and the parameter synthesis are done.
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#include "Sapo.h"

/**
 * Constructor that instantiates Sapo
 *
 * @param[in] model model to analyize
 * @param[in] sapo_opt options to tune sapo
 */
Sapo::Sapo(Model *model, sapo_opt options):
    dyns(model->getDyns()), vars(model->getVars()), params(model->getParams()),
    options(options)
{
}

/**
 * Reachable set computation
 *
 * @param[in] initSet bundle with the initial set
 * @param[in] k time horizon
 * @returns flowpipe of bundles
 */
Flowpipe Sapo::reach(const Bundle &initSet, unsigned int k)
{
  using namespace std;
  using namespace GiNaC;

  map<vector<int>, pair<lst, lst>> controlPts;

  Flowpipe flowpipe(initSet.getDirectionMatrix());

  if (this->options.verbose) {
    cout << initSet.getLinearSystem() << endl << endl;
  }

  cout << "Computing reach set..." << flush;

  Bundle X = initSet;
  LinearSystem Xls = initSet.getLinearSystem();
  flowpipe.append(Xls);

  unsigned int i = 0;
  while (i < k && !Xls.isEmpty()) {
    i++;

    X = X.transform(this->vars, this->dyns, controlPts,
                    this->options.trans); // transform it

    if (this->options.decomp > 0) { // eventually decompose it
      X = X.decompose(this->options.alpha, this->options.decomp);
    }

    flowpipe.append(X.getLinearSystem()); // store result

    if (this->options.verbose) {
      cout << flowpipe.get(i) << endl << endl;
    }
  }

  cout << "done" << endl;

  return flowpipe;
}

/**
 * Reachable set computation for parameteric dynamical systems
 *
 * @param[in] initSet bundle with the initial set
 * @param[in] paraSet set of parameters
 * @param[in] k time horizon
 * @returns flowpipe of bundles
 */
Flowpipe Sapo::reach(const Bundle &initSet, LinearSystemSet &paraSet,
                     unsigned int k)
{
  using namespace std;
  using namespace GiNaC;

  map<vector<int>, pair<lst, lst>> controlPts;

  cout << "Parameter set" << endl
       << paraSet << endl
       << endl
       << endl
       << "Computing parametric reach set..." << flush;

  if (this->options.verbose) {
    cout << initSet.getLinearSystem() << endl << endl;
  }

  std::vector<Bundle> cbundles(paraSet.size(), initSet);
  LinearSystemSet Xls(
      std::make_shared<LinearSystem>(initSet.getLinearSystem()));

  Flowpipe flowpipe(initSet.getDirectionMatrix());
  flowpipe.append(initSet);

  unsigned int i = 0;
  while (i < k && !Xls.isEmpty()) {
    i++;

    LinearSystemSet Xls;

    LinearSystemSet::iterator pset_it(paraSet.begin());
    for (unsigned int b_idx = 0; b_idx < cbundles.size(); b_idx++) {
      cbundles[b_idx] = cbundles[b_idx].transform(
          this->vars, this->params, this->dyns, *(pset_it++), controlPts,
          this->options.trans); // transform it

      if (this->options.decomp > 0) { // eventually decompose it
        cbundles[b_idx] = cbundles[b_idx].decompose(this->options.alpha,
                                                    this->options.decomp);
      }
      Xls.add(cbundles[b_idx].getLinearSystem());
    }

    flowpipe.append(Xls); // store result

    if (this->options.verbose) {
      cout << flowpipe.get(i) << endl << endl;
    }
  }

  cout << "done" << endl;

  return flowpipe;
}

/**
 * Parameter synthesis procedure
 *
 * @param[in] reachSet bundle with the initial set
 * @param[in] parameterSet set of sets of parameters
 * @param[in] formula STL contraint to impose over the model
 * @param[in] max_splits maximum number of splits of the original
 *                       parameter set to identify a non-null solution
 * @returns refined sets of parameters
 */
LinearSystemSet Sapo::synthesize(Bundle &reachSet,
                                 LinearSystemSet &parameterSet,
                                 const std::shared_ptr<STL> formula,
                                 const unsigned int max_splits)
{
  using namespace std;

  cout << "Initial parameter set" << endl;
  parameterSet.print();

  cout << endl << "Specification: ";
  formula->print();

  cout << endl << endl << "Synthesizing parameters..." << flush;

  unsigned int num_of_splits = 0;
  LinearSystemSet splitParamSet = parameterSet;
  LinearSystemSet res = this->synthesizeSTL(reachSet, parameterSet, formula);

  while (res.isEmpty() && num_of_splits++ < max_splits) {
    parameterSet = splitParamSet.get_a_finer_covering();
    splitParamSet = parameterSet;

    res = this->synthesizeSTL(reachSet, parameterSet, formula);
  }

  res.simplify();

  cout << "done" << endl;

  return res;
}

/**
 * Internal parameter synthesis procedure
 *
 * @param[in] reachSet bundle with the initial set
 * @param[in] parameterSet set of sets of parameters
 * @param[in] formula STL contraint to impose over the model
 * @returns refined sets of parameters
 */
LinearSystemSet Sapo::synthesizeSTL(Bundle &reachSet,
                                    LinearSystemSet &parameterSet,
                                    const std::shared_ptr<STL> formula)
{
  switch (formula->getType()) {

  // Atomic predicate
  case ATOM:
    return this->refineParameters(reachSet, parameterSet,
                                  std::dynamic_pointer_cast<Atom>(formula));

  // Conjunction
  case CONJUNCTION: {
    const std::shared_ptr<Conjunction> conj
        = std::dynamic_pointer_cast<Conjunction>(formula);
    LinearSystemSet LS1 = this->synthesizeSTL(reachSet, parameterSet,
                                              conj->getLeftSubFormula());
    LinearSystemSet LS2 = this->synthesizeSTL(reachSet, parameterSet,
                                              conj->getRightSubFormula());
    return intersection(LS1, LS2);
  }

  // Disjunction
  case DISJUNCTION: {
    const std::shared_ptr<Disjunction> disj
        = std::dynamic_pointer_cast<Disjunction>(formula);
    LinearSystemSet LS1 = this->synthesizeSTL(reachSet, parameterSet,
                                              disj->getLeftSubFormula());
    LS1.unionWith(this->synthesizeSTL(reachSet, parameterSet,
                                      disj->getRightSubFormula()));

    return LS1;
  }

  // Until
  case UNTIL:
    return this->synthesizeUntil(reachSet, parameterSet,
                                 std::dynamic_pointer_cast<Until>(formula), 0);

  // Always
  case ALWAYS:
    return this->synthesizeAlways(
        reachSet, parameterSet, std::dynamic_pointer_cast<Always>(formula), 0);

  // Eventually
  case EVENTUALLY: {
    std::shared_ptr<Atom> a = std::make_shared<Atom>(-1);
    const std::shared_ptr<Eventually> ev
        = std::dynamic_pointer_cast<Eventually>(formula);

    std::shared_ptr<Until> u = std::make_shared<Until>(
        a, ev->time_bounds().begin(), ev->time_bounds().end(),
        ev->getSubFormula());
    return this->synthesizeUntil(reachSet, parameterSet, u, 0);
  }
  default:
    throw std::logic_error("Unsupported formula");
  }
}

/**
 * Parameter synthesis w.r.t. an atomic formula
 *
 * @param[in] reachSet bundle with the initial set
 * @param[in] parameterSet set of sets of parameters
 * @param[in] sigma STL atomic formula
 * @returns refined sets of parameters
 */
LinearSystemSet Sapo::refineParameters(Bundle &reachSet,
                                       LinearSystemSet &parameterSet,
                                       const std::shared_ptr<Atom> atom)
{
  using namespace std;
  using namespace GiNaC;

  LinearSystemSet result;

  for (unsigned int i = 0; i < reachSet.getCard();
       i++) { // for each parallelotope

    // complete the key
    vector<int> key = reachSet.getTemplate(i);
    key.push_back(atom->getID());

    Parallelotope P = reachSet.getParallelotope(i);
    lst genFun = P.getGeneratorFunction();
    lst controlPts;

    if (this->synthControlPts.count(key) == 0
        || (!this->synthControlPts[key].first.is_equal(genFun))) {
      // compose f(gamma(x))
      lst sub, fog;
      for (unsigned int j = 0; j < this->vars.nops(); j++) {
        sub.append(vars[j] == genFun[j]);
      }
      for (unsigned int j = 0; j < vars.nops(); j++) {
        fog.append(this->dyns[j].subs(sub));
      }

      // compose sigma(f(gamma(x)))
      lst sub_sigma;
      for (unsigned int j = 0; j < this->vars.nops(); j++) {
        sub_sigma.append(vars[j] == fog[j]);
      }
      ex sofog;
      sofog = atom->getPredicate().subs(sub_sigma);

      // compute the Bernstein control points
      controlPts = BaseConverter(P.getAlpha(), sofog).getBernCoeffsMatrix();
      this->synthControlPts[key].first = genFun;
      this->synthControlPts[key].second = controlPts;

    } else {
      controlPts = this->synthControlPts[key].second;
    }

    // substitute numerical values in sofog
    vector<double> base_vertex = P.getBaseVertex();
    vector<double> lengths = P.getLenghts();

    lst qvars(P.getQ());
    lst bvars(P.getBeta());
    lst para_sub;
    for (unsigned int j = 0; j < this->vars.nops(); j++) {
      para_sub.append(qvars[j] == base_vertex[j]);
      para_sub.append(bvars[j] == lengths[j]);
    }
    ex num_sofog;
    lst synth_controlPts;
    // for (int j=0; j<controlPts.nops(); j++) {
    for (lst::const_iterator j = controlPts.begin(); j != controlPts.end();
         ++j) {
      synth_controlPts.append((*j).subs(para_sub));
    }

    // cout<<synth_controlPts;

    LinearSystem num_constraintLS(this->params, synth_controlPts);
    result.unionWith(intersection(parameterSet, num_constraintLS));
  }

  return result;
}

/**
 * Parameter synthesis w.r.t. an until formula
 *
 * @param[in] reachSet bundle with the initial set
 * @param[in] parameterSet set of sets of parameters
 * @param[in] sigma STL until formula
 * @param[in] time is the time of the current evaluation
 * @returns refined sets of parameters
 */
LinearSystemSet Sapo::synthesizeUntil(Bundle &reachSet,
                                      LinearSystemSet &parameterSet,
                                      const std::shared_ptr<Until> formula,
                                      const int time)
{
  LinearSystemSet result;
  const TimeInterval &t_itvl = formula->time_bounds();

  // Base case
  if (t_itvl.isEmpty())
    return result;

  // Until interval far
  if (t_itvl > time) {
    // Synthesize wrt phi1
    LinearSystemSet P1 = this->synthesizeSTL(reachSet, parameterSet,
                                             formula->getLeftSubFormula());
    if (P1.isEmpty()) {
      return P1; // false until
    } else {
      // Reach step wrt to the i-th linear system of P1
      for (LinearSystemSet::iterator P1_it = P1.begin(); P1_it != P1.end();
           ++P1_it) {
        // TODO : add the decomposition
        Bundle newReachSet
            = reachSet.transform(this->vars, this->params, this->dyns, *P1_it,
                                 this->reachControlPts, this->options.trans);

        // TODO: Check whether the object tmpLSset can be removed
        LinearSystemSet tmpLSset(*P1_it);
        result.unionWith(
            synthesizeUntil(newReachSet, tmpLSset, formula, time + 1));
      }
      return result;
    }
  }

  // Inside until interval
  if (t_itvl.end() > time) {
    // Refine wrt phi1 and phi2
    LinearSystemSet P1 = this->synthesizeSTL(reachSet, parameterSet,
                                             formula->getLeftSubFormula());

    if (P1.isEmpty()) {
      return this->synthesizeSTL(reachSet, parameterSet,
                                 formula->getRightSubFormula());
    }

    // shift until interval
    for (LinearSystemSet::iterator P1_it = P1.begin(); P1_it != P1.end();
         ++P1_it) {
      // 	TODO : add decomposition
      Bundle newReachSet
          = reachSet.transform(this->vars, this->params, this->dyns, *P1_it,
                               this->reachControlPts, this->options.trans);
      LinearSystemSet tmpLSset(*P1_it);

      /* TODO: Add a constructor that copy the parameter set */
      result.unionWith(
          synthesizeUntil(newReachSet, tmpLSset, formula, time + 1));
    }
    result.unionWith(this->synthesizeSTL(reachSet, parameterSet,
                                         formula->getRightSubFormula()));

    return result;
  }

  // If none of the above condition holds, then it must holds that :
  // 			t_itvl.begin()<=time and t_itvl.end()==time
  return this->synthesizeSTL(reachSet, parameterSet,
                             formula->getRightSubFormula());
}

/**
 * Parameter synthesis w.r.t. an always formula
 *
 * @param[in] reachSet bundle with the initial set
 * @param[in] parameterSet set of sets of parameters
 * @param[in] sigma STL always formula
 * @returns refined sets of parameters
 */
LinearSystemSet Sapo::synthesizeAlways(Bundle &reachSet,
                                       LinearSystemSet &parameterSet,
                                       const std::shared_ptr<Always> formula,
                                       const int time)
{
  LinearSystemSet result;
  const TimeInterval &t_itvl = formula->time_bounds();

  // Base case
  if (t_itvl.isEmpty()) {
    return result;
  }

  // Always interval far
  if (t_itvl > time) {
    // Reach step wrt to the i-th linear system of parameterSet
    for (LinearSystemSet::iterator pset_it = parameterSet.begin();
         pset_it != parameterSet.end(); ++pset_it) {

      Bundle newReachSet
          = reachSet.transform(this->vars, this->params, this->dyns, *pset_it,
                               this->reachControlPts, options.trans);
      LinearSystemSet tmpLSset(*pset_it);

      /* TODO: Add a constructor that copy the parameter set */
      result.unionWith(
          synthesizeAlways(newReachSet, tmpLSset, formula, time + 1));
    }
    return result;
  }

  // Inside Always interval
  if (t_itvl.end() > time) {

    // Refine wrt phi
    LinearSystemSet P = this->synthesizeSTL(reachSet, parameterSet,
                                            formula->getSubFormula());

    if (!P.isEmpty()) {
      // Reach step wrt to the i-th linear system of P
      for (LinearSystemSet::iterator P_it = P.begin(); P_it != P.end();
           ++P_it) {
        Bundle newReachSet
            = reachSet.transform(this->vars, this->params, this->dyns, *P_it,
                                 this->reachControlPts, options.trans);

        /* TODO: Add a constructor that copy a const parameter set */
        LinearSystemSet tmpLSset(*P_it);
        result.unionWith(
            synthesizeAlways(newReachSet, tmpLSset, formula, time + 1));
      }

      return result;
    }

    return P;
  }

  // If none of the above condition holds, then it must holds that :
  // 			t_itvl.begin()<=time and t_itvl.end()==time
  return this->synthesizeSTL(reachSet, parameterSet, formula->getSubFormula());
}

Sapo::~Sapo()
{
  // TODO Auto-generated destructor stub
}
