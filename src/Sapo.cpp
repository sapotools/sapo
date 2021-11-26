/**
 * @file Sapo.cpp
 * Core of Sapo tool.
 * Here the reachable set and the parameter synthesis are done.
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#include "Sapo.h"

#if WITH_THREADS
#include <thread>
#include <shared_mutex>

#endif // WITH_THREADS

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

  Flowpipe flowpipe(initSet.get_directions());

  Polytope Xls = initSet;

  if (this->verbose) {
    cout << "Initial Set" << endl
         << Xls << endl
         << endl
         << "Computing reach set..." << flush;
  }

  Bundle X = initSet;
  flowpipe.append(Xls);

  unsigned int i = 0;
  while (i < k && !Xls.is_empty()) {
    i++;

    X = X.transform(this->vars, this->dyns,
                    this->trans); // transform it

    if (this->decomp > 0) { // if requested, decompose it
      X = X.decompose(this->decomp_weight, this->decomp);
    }

    flowpipe.append(X); // store result

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
Flowpipe Sapo::reach(const Bundle &initSet, const PolytopesUnion &pSet,
                     unsigned int k)
{
  using namespace std;

  if (this->verbose) {
    cout << "Initial Set" << endl
         << (Polytope)initSet << endl
         << endl
         << "Parameter set" << endl
         << pSet << endl
         << endl
         << "Computing parametric reach set...";
  }

  std::list<Bundle> cbundles{initSet};
  PolytopesUnion last_step;

  last_step.add(initSet);

  Flowpipe flowpipe(initSet.get_directions());
  flowpipe.append(initSet);

  unsigned int i = 0;

  // while time horizon has not been reached and last step reached is not empty
  // TODO: check whether there exists any chance for the last_step to be empty
  while (i < k && last_step.size() != 0) {

    // create a new last step reach set
    last_step = PolytopesUnion();
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
                              this->trans); // transform it

        if (this->decomp > 0) { // if requested, decompose it
          bundle = bundle.decompose(this->decomp_weight, this->decomp);
        }

        Polytope bls = bundle;

        // TODO: check whether there is any chance for a transformed bundle to
        // be empty
        if (!bls.is_empty()) {
          // add to the new reached bundle
          nbundles.push_back(bundle);

          last_step.add(bls);
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

std::list<PolytopesUnion>
get_a_finer_covering(const std::list<PolytopesUnion> &orig)
{
  std::list<PolytopesUnion> result;
  for (auto ps_it = std::cbegin(orig); ps_it != std::cend(orig); ++ps_it) {

    // if the polytopes union contains more than one polytope
    switch (ps_it->size()) {
    case 0: // the polytopes union does not contain any polytope

      // nothing to add to the resulting list
      break;
    case 1: // the polytopes union contains exacly one polytope
    {       // then, split it by using Polytope::get_a_finer_covering();

      std::list<Polytope> f_cov = (ps_it->begin())->split();
      for (auto ls_it = std::begin(f_cov); ls_it != std::end(f_cov); ++ls_it) {
        result.push_back(*ls_it);
      }
    } break;
    case 2: // the polytopes union contains more than one polytope
    {       // then, unpack them
      for (auto ls_it = ps_it->begin(); ls_it != ps_it->end(); ++ls_it) {
        result.push_back(*ls_it);
      }
    } break;
    }
  }

  return result;
}

#if WITH_THREADS
template<typename T>
class ThreadSafeList
{
  std::list<T> list;
  mutable std::shared_timed_mutex mutex;

public:
  ThreadSafeList(): list() {}

  ThreadSafeList(const std::list<T> &list): list(list) {}

  ThreadSafeList<T> &push_back(T &&obj)
  {
    std::shared_lock<std::shared_timed_mutex> writelock(mutex,
                                                        std::defer_lock);

    list.push_back(obj);

    return *this;
  }

  ThreadSafeList<T> &push_back(const T &obj)
  {
    std::shared_lock<std::shared_timed_mutex> writelock(mutex,
                                                        std::defer_lock);

    list.push_back(obj);

    return *this;
  }

  const std::list<T> &get_list() const
  {
    return list;
  }
};
#endif // WITH_THREADS

/**
 * Parameter synthesis
 *
 * @param[in] sapo the Sapo object that performs the computation
 * @param[in] reachSet bundle representing the current reached set
 * @param[in] pSetList a list of parameter sets
 * @param[in] formula is an STL formula providing the specification
 * @returns the list of refined parameter sets
 */
std::list<PolytopesUnion>
synthesize_list(Sapo &sapo, const Bundle &reachSet,
                const std::list<PolytopesUnion> &pSetList,
                const std::shared_ptr<STL> &formula)
{

#if WITH_THREADS && false

  ThreadSafeList<PolytopesUnion> results;
  auto synthesize_funct
      = [&results, &sapo, &reachSet, &formula](const PolytopesUnion &pSet) {
          results.push_back(sapo.synthesize(reachSet, pSet, formula));
        };

  std::vector<std::thread> threads;
  for (auto ps_it = std::begin(pSetList); ps_it != std::end(pSetList);
       ++ps_it) {
    threads.push_back(std::thread(synthesize_funct, std::ref(*ps_it)));
  }

  for (std::thread &th: threads) {
    if (th.joinable())
      th.join();
  }

  return results.get_list();

#else // WITH_THREADS
  std::list<PolytopesUnion> results;

  std::vector<std::thread> threads;
  for (auto ps_it = std::begin(pSetList); ps_it != std::end(pSetList);
       ++ps_it) {
    results.push_back(sapo.synthesize(reachSet, *ps_it, formula));
  }

  return results;

#endif // WITH_THREADS
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
std::list<PolytopesUnion> Sapo::synthesize(const Bundle &reachSet,
                                           const PolytopesUnion &pSet,
                                           const std::shared_ptr<STL> formula,
                                           const unsigned int max_splits)
{
  using namespace std;

  std::list<PolytopesUnion> pSetList{pSet};

  unsigned int num_of_splits = 0;
  std::list<PolytopesUnion> res
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
PolytopesUnion Sapo::synthesize(const Bundle &reachSet,
                                const PolytopesUnion &pSet,
                                const std::shared_ptr<Conjunction> conj)
{
  PolytopesUnion Pu1
      = this->synthesize(reachSet, pSet, conj->getLeftSubFormula());
  PolytopesUnion Pu2
      = this->synthesize(reachSet, pSet, conj->getRightSubFormula());
  return intersect(Pu1, Pu2);
}

/**
 * Parmeter synthesis for disjunctions
 *
 * @param[in] reachSet bundle representing the current reached set
 * @param[in] pSet the current parameter set
 * @param[in] conj is an STL disjunction providing the specification
 * @returns refined parameter set
 */
PolytopesUnion Sapo::synthesize(const Bundle &reachSet,
                                const PolytopesUnion &pSet,
                                const std::shared_ptr<Disjunction> disj)
{
  PolytopesUnion Pu
      = this->synthesize(reachSet, pSet, disj->getLeftSubFormula());
  Pu.add(this->synthesize(reachSet, pSet, disj->getRightSubFormula()));

  return Pu;
}

/**
 * Parmeter synthesis for the eventually fomulas
 *
 * @param[in] reachSet bundle representing the current reached set
 * @param[in] pSet the current parameter set
 * @param[in] ev is an STL eventually formula providing the specification
 * @returns refined parameter set
 */
PolytopesUnion Sapo::synthesize(const Bundle &reachSet,
                                const PolytopesUnion &pSet,
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
PolytopesUnion Sapo::synthesize(const Bundle &reachSet,
                                const PolytopesUnion &pSet,
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
PolytopesUnion Sapo::synthesize(const Bundle &reachSet,
                                const PolytopesUnion &pSet,
                                const std::shared_ptr<Atom> atom)
{
  using namespace std;
  using namespace GiNaC;

  PolytopesUnion result;

  for (unsigned int i = 0; i < reachSet.num_of_templates();
       i++) { // for each parallelotope

    Parallelotope P = reachSet.getParallelotope(i);
    lst genFun = build_instanciated_generator_functs(reachSet.get_alpha(), P);

    const lst fog = sub_vars(this->dyns, vars, genFun);

    // compose sigma(f(gamma(x)))
    lst sub_sigma;
    for (unsigned int j = 0; j < this->vars.nops(); j++) {
      sub_sigma.append(vars[j] == fog[j]);
    }

    const ex sofog = atom->getPredicate().subs(sub_sigma);

    // compute the Bernstein control points
    lst controlPts
        = BaseConverter(reachSet.get_alpha(), sofog).getBernCoeffsMatrix();

    Polytope constraints(this->params, controlPts);
    result.add(intersect(pSet, constraints));
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
PolytopesUnion Sapo::synthesize(const Bundle &reachSet,
                                const PolytopesUnion &pSet,
                                const std::shared_ptr<Until> formula,
                                const int time)
{
  const TimeInterval &t_itvl = formula->time_bounds();

  // Base case
  if (t_itvl.is_empty())
    return PolytopesUnion();

  // Until interval far
  if (t_itvl > time) {
    // Synthesize wrt phi1
    PolytopesUnion P1
        = this->synthesize(reachSet, pSet, formula->getLeftSubFormula());
    if (P1.is_empty()) {
      return P1; // false until
    } else {
      return transition_and_synthesis(reachSet, P1, formula, time);
    }
  }

  // Inside until interval
  if (t_itvl.end() > time) {
    // Refine wrt phi1 and phi2
    PolytopesUnion P1
        = this->synthesize(reachSet, pSet, formula->getLeftSubFormula());

    if (P1.is_empty()) {
      return this->synthesize(reachSet, pSet, formula->getRightSubFormula());
    }

    PolytopesUnion result
        = transition_and_synthesis(reachSet, P1, formula, time);
    result.add(
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
PolytopesUnion Sapo::synthesize(const Bundle &reachSet,
                                const PolytopesUnion &pSet,
                                const std::shared_ptr<Always> formula,
                                const int time)
{
  const TimeInterval &t_itvl = formula->time_bounds();

  // Base case
  if (t_itvl.is_empty()) {
    return PolytopesUnion();
  }

  // Always interval far
  if (t_itvl > time) {
    // transition and synthesis
    return transition_and_synthesis(reachSet, pSet, formula, time);
  }

  // Inside Always interval
  if (t_itvl.end() > time) {

    // Refine wrt phi
    PolytopesUnion P
        = this->synthesize(reachSet, pSet, formula->getSubFormula());

    if (P.is_empty()) {
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
