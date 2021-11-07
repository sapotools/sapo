/**
 * @file Sapo.cpp
 * Core of Sapo tool.
 * Here the reachable set and the parameter synthesis are done.
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#include "Sapo.h"

#include <thread>

/**
 * Constructor that instantiates Sapo
 *
 * @param[in] model model to analyize
 * @param[in] sapo_opt options to tune sapo
 */
Sapo::Sapo(Model *model):
    dyns(model->getDyns()), vars(model->getVars()), params(model->getParams())
{
}

/**
 * Reachable set computation
 *
 * @param[in] initSet bundle representing the current reached set
 * @param[in] k time horizon
 * @returns the reached flowpipe
 */
Flowpipe Sapo::reach(const Bundle &initSet, unsigned int k)
{
  using namespace std;
  using namespace GiNaC;

  ControlPointStorage controlPts;

  Flowpipe flowpipe(initSet.getDirectionMatrix());

  if (this->verbose) {
    cout << "Initial Set" << endl
         << initSet.getLinearSystem() << endl
         << endl
         << "Computing reach set..." << flush;
  }

  Bundle X = initSet;
  LinearSystem Xls = initSet.getLinearSystem();
  flowpipe.append(Xls);

  unsigned int i = 0;
  while (i < k && !Xls.isEmpty()) {
    i++;

    X = X.transform(this->vars, this->dyns, controlPts,
                    this->trans); // transform it

    if (this->decomp > 0) { // if requested, decompose it
      X = X.decompose(this->alpha, this->decomp);
    }

    flowpipe.append(X.getLinearSystem()); // store result

    if (this->verbose) {
      cout << flowpipe.get(i) << endl << endl;
    }
  }

  if (this->verbose) {
    cout << "done" << endl;
  }

  return flowpipe;
}

/**
 * Reachable set computation for parameteric dynamical systems
 *
 * @param[in] initSet bundle representing the current reached set
 * @param[in] pSet the set of parameters
 * @param[in] k time horizon
 * @returns the reached flowpipe
 */
Flowpipe Sapo::reach(const Bundle &initSet, const LinearSystemSet &pSet,
                     unsigned int k)
{
  using namespace std;
  using namespace GiNaC;

  if (this->verbose) {
    cout << "Initial Set" << endl
         << initSet.getLinearSystem() << endl
         << endl
         << "Parameter set" << endl
         << pSet << endl
         << endl
         << "Computing parametric reach set...";
  }

  std::list<Bundle> cbundles{initSet};
  ControlPointStorage ctrlPts;
  LinearSystemSet last_step;

  last_step.add(initSet.getLinearSystem());

  Flowpipe flowpipe(initSet.getDirectionMatrix());
  flowpipe.append(initSet);

  unsigned int i = 0;

  // while time horizon has not been reached and last step reached is not empty
  // TODO: check whether there exists any chance for the last_step to be empty
  while (i < k && last_step.size() != 0) {

    // create a new last step reach set
    last_step = LinearSystemSet();
    i++;

    // create a list for the new reached bundles
    std::list<Bundle> nbundles;

    // for all the old bundles
    for (auto b_it = std::cbegin(cbundles); b_it != std::cend(cbundles);
         ++b_it) {
      // for all the parameter sets
      for (auto p_it = pSet.begin(); p_it != pSet.end(); ++p_it) {

        // get the transformed bundle
        Bundle bundle
            = b_it->transform(this->vars, this->params, this->dyns, *p_it,
                              ctrlPts, this->trans); // transform it

        if (this->decomp > 0) { // if requested, decompose it
          bundle = bundle.decompose(this->alpha, this->decomp);
        }

        LinearSystem bls = bundle.getLinearSystem();

        // TODO: check whether there is any chance for a transformed bundle to
        // be empty
        if (!bls.isEmpty()) {
          // add to the new reached bundle
          nbundles.push_back(bundle);

          last_step.add(bundle.getLinearSystem());
        }
      }
    }

    // swap current bundles and new bundles
    swap(cbundles, nbundles);

    // add the last step to the flow pipe
    flowpipe.append(last_step); // store result

    if (this->verbose) {
      cout << flowpipe.get(i) << endl << endl;
    }
  }

  if (this->verbose) {
    cout << "done" << endl;
  }

  return flowpipe;
}

std::list<LinearSystemSet>
get_a_finer_covering(const std::list<LinearSystemSet> &orig)
{
  std::list<LinearSystemSet> result;
  for (auto ps_it = std::cbegin(orig); ps_it != std::cend(orig); ++ps_it) {

    // if the linear system set contains more than one linear system
    switch (ps_it->size()) {
    case 0: // the linear system set does not contain any linear system

      // nothing to add to the resulting list
      break;
    case 1: // the linear system set contains exacly one linear system
    {       // then, split it by using LinearSystem::get_a_finer_covering();

      std::list<LinearSystem> f_cov = (ps_it->begin())->get_a_finer_covering();
      for (auto ls_it = std::begin(f_cov); ls_it != std::end(f_cov); ++ls_it) {
        result.push_back(*ls_it);
      }
    } break;
    case 2: // the linear system set contains more than one linear system
    {       // then, unpack them
      for (auto ls_it = ps_it->begin(); ls_it != ps_it->end(); ++ls_it) {
        result.push_back(*ls_it);
      }
    } break;
    }
  }

  return result;
}

/**
 * Parameter synthesis
 *
 * @param[in] sapo the Sapo object that performs the computation
 * @param[in] reachSet bundle representing the current reached set
 * @param[in] pSetList a list of parameter sets
 * @param[in] formula is an STL formula providing the specification
 * @returns the list of refined parameter sets
 */
std::list<LinearSystemSet>
synthesize_list(Sapo &sapo, Bundle reachSet,
                const std::list<LinearSystemSet> &pSetList,
                const std::shared_ptr<STL> &formula)
{
  std::list<LinearSystemSet> results;

  for (auto ps_it = std::begin(pSetList); ps_it != std::end(pSetList);
       ++ps_it) {
    results.push_back(sapo.synthesize(reachSet, *ps_it, formula));
  }

  return results;
}

/**
 * Parameter synthesis with splits
 *
 * @param[in] reachSet bundle representing the current reached set
 * @param[in] pSet is the current parameter sets
 * @param[in] formula is an STL formula providing the specification
 * @param[in] max_splits maximum number of splits of the original
 *                       parameter set to identify a non-null solution
 * @returns the list of refined parameter sets
 */
std::list<LinearSystemSet> Sapo::synthesize(const Bundle &reachSet,
                                            const LinearSystemSet &pSet,
                                            const std::shared_ptr<STL> formula,
                                            const unsigned int max_splits)
{
  using namespace std;

  std::list<LinearSystemSet> pSetList{pSet};

  unsigned int num_of_splits = 0;
  std::list<LinearSystemSet> res
      = synthesize_list(*this, reachSet, pSetList, formula);

  while (every_set_is_empty(res) && num_of_splits++ < max_splits) {
    res.clear();
    pSetList = get_a_finer_covering(pSetList);

    res = synthesize_list(*this, reachSet, pSetList, formula);
  }

  for (auto lss_it = std::begin(res); lss_it != std::end(res); ++lss_it) {
    lss_it->simplify();
  }

  if (this->verbose) {
    cout << "done" << endl;
  }

  return res;
}

/**
 * Parmeter synthesis for conjunctions
 *
 * @param[in] reachSet bundle representing the current reached set
 * @param[in] pSet the current parameter set
 * @param[in] conj is an STL conjunction providing the specification
 * @returns refined parameter set
 */
LinearSystemSet Sapo::synthesize(const Bundle &reachSet,
                                 const LinearSystemSet &pSet,
                                 const std::shared_ptr<Conjunction> conj)
{
  LinearSystemSet LS1
      = this->synthesize(reachSet, pSet, conj->getLeftSubFormula());
  LinearSystemSet LS2
      = this->synthesize(reachSet, pSet, conj->getRightSubFormula());
  return intersection(LS1, LS2);
}

/**
 * Parmeter synthesis for disjunctions
 *
 * @param[in] reachSet bundle representing the current reached set
 * @param[in] pSet the current parameter set
 * @param[in] conj is an STL disjunction providing the specification
 * @returns refined parameter set
 */
LinearSystemSet Sapo::synthesize(const Bundle &reachSet,
                                 const LinearSystemSet &pSet,
                                 const std::shared_ptr<Disjunction> disj)
{
  LinearSystemSet LS1
      = this->synthesize(reachSet, pSet, disj->getLeftSubFormula());
  LS1.unionWith(this->synthesize(reachSet, pSet, disj->getRightSubFormula()));

  return LS1;
}

/**
 * Parmeter synthesis for the eventually fomulas
 *
 * @param[in] reachSet bundle representing the current reached set
 * @param[in] pSet the current parameter set
 * @param[in] ev is an STL eventually formula providing the specification
 * @returns refined parameter set
 */
LinearSystemSet Sapo::synthesize(const Bundle &reachSet,
                                 const LinearSystemSet &pSet,
                                 const std::shared_ptr<Eventually> ev)
{
  std::shared_ptr<Atom> true_atom = std::make_shared<Atom>(-1);

  std::shared_ptr<Until> u
      = std::make_shared<Until>(true_atom, ev->time_bounds().begin(),
                                ev->time_bounds().end(), ev->getSubFormula());

  return this->synthesize(reachSet, pSet, u, 0);
}

/**
 * Parameter synthesis method
 *
 * @param[in] reachSet bundle representing the current reached set
 * @param[in] pSet the current parameter set
 * @param[in] formula is an STL specification for the model
 * @returns refined parameter set
 */
LinearSystemSet Sapo::synthesize(const Bundle &reachSet,
                                 const LinearSystemSet &pSet,
                                 const std::shared_ptr<STL> formula)
{
  switch (formula->getType()) {

  // Atomic predicate
  case ATOM:
    return synthesize(reachSet, pSet,
                      std::dynamic_pointer_cast<Atom>(formula));

  // Conjunction
  case CONJUNCTION:
    return synthesize(reachSet, pSet,
                      std::dynamic_pointer_cast<Conjunction>(formula));

  // Disjunction
  case DISJUNCTION:
    return synthesize(reachSet, pSet,
                      std::dynamic_pointer_cast<Disjunction>(formula));

  // Until
  case UNTIL:
    return this->synthesize(reachSet, pSet,
                            std::dynamic_pointer_cast<Until>(formula), 0);

  // Always
  case ALWAYS:
    return synthesize(reachSet, pSet,
                      std::dynamic_pointer_cast<Always>(formula), 0);

  // Eventually
  case EVENTUALLY:
    return synthesize(reachSet, pSet,
                      std::dynamic_pointer_cast<Eventually>(formula));

  default:
    throw std::logic_error("Unsupported formula");
  }
}

// TODO: the following method is too long and must be split.
/**
 * Parameter synthesis for atomic formulas
 *
 * @param[in] reachSet bundle representing the current reached set
 * @param[in] pSet the current set of parameters
 * @param[in] sigma STL atomic formula
 * @returns refined parameter set
 */
LinearSystemSet Sapo::synthesize(const Bundle &reachSet,
                                 const LinearSystemSet &pSet,
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

    if (!(this->synthControlPts.contains(key)
          && this->synthControlPts.gen_fun_is_equal_to(key, genFun))) {
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
      this->synthControlPts.set(key, genFun, controlPts);

    } else {
      controlPts = this->synthControlPts.get_ctrl_pts(key);
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
    result.unionWith(intersection(pSet, num_constraintLS));
  }

  return result;
}

/**
 * Parameter synthesis for until formulas
 *
 * @param[in] reachSet bundle representing the current reached set
 * @param[in] pSet the current set of parameters
 * @param[in] formula STL until formula
 * @param[in] time is the time of the current evaluation
 * @returns refined parameter set
 */
LinearSystemSet Sapo::synthesize(const Bundle &reachSet,
                                 const LinearSystemSet &pSet,
                                 const std::shared_ptr<Until> formula,
                                 const int time)
{
  const TimeInterval &t_itvl = formula->time_bounds();

  // Base case
  if (t_itvl.isEmpty())
    return LinearSystemSet();

  // Until interval far
  if (t_itvl > time) {
    // Synthesize wrt phi1
    LinearSystemSet P1
        = this->synthesize(reachSet, pSet, formula->getLeftSubFormula());
    if (P1.isEmpty()) {
      return P1; // false until
    } else {
      return transition_and_synthesis(reachSet, P1, formula, time);
    }
  }

  // Inside until interval
  if (t_itvl.end() > time) {
    // Refine wrt phi1 and phi2
    LinearSystemSet P1
        = this->synthesize(reachSet, pSet, formula->getLeftSubFormula());

    if (P1.isEmpty()) {
      return this->synthesize(reachSet, pSet, formula->getRightSubFormula());
    }

    LinearSystemSet result
        = transition_and_synthesis(reachSet, P1, formula, time);
    result.unionWith(
        this->synthesize(reachSet, pSet, formula->getRightSubFormula()));

    return result;
  }

  // If none of the above condition holds, then it must holds that :
  // 			t_itvl.begin()<=time and t_itvl.end()==time
  return this->synthesize(reachSet, pSet, formula->getRightSubFormula());
}

/**
 * Parameter synthesis for always formulas
 *
 * @param[in] reachSet bundle representing the current reached set
 * @param[in] pSet the current set of parameters
 * @param[in] formula STL always formula
 * @param[in] time is the time of the current evaluation
 * @returns refined parameter set
 */
LinearSystemSet Sapo::synthesize(const Bundle &reachSet,
                                 const LinearSystemSet &pSet,
                                 const std::shared_ptr<Always> formula,
                                 const int time)
{
  const TimeInterval &t_itvl = formula->time_bounds();

  // Base case
  if (t_itvl.isEmpty()) {
    return LinearSystemSet();
  }

  // Always interval far
  if (t_itvl > time) {
    // transition and synthesis
    return transition_and_synthesis(reachSet, pSet, formula, time);
  }

  // Inside Always interval
  if (t_itvl.end() > time) {

    // Refine wrt phi
    LinearSystemSet P
        = this->synthesize(reachSet, pSet, formula->getSubFormula());

    if (P.isEmpty()) {
      return P;
    }

    return transition_and_synthesis(reachSet, P, formula, time);
  }

  // If none of the above condition holds, then it must holds that :
  // 			t_itvl.begin()<=time and t_itvl.end()==time
  return this->synthesize(reachSet, pSet, formula->getSubFormula());
}

Sapo::~Sapo()
{
  // TODO Auto-generated destructor stub
}
