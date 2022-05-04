/**
 * @file Sapo.cpp
 * Core of Sapo tool.
 * Here the reachable set and the parameter synthesis are done.
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#include "Sapo.h"

#include <limits>

#ifdef WITH_THREADS
#include <shared_mutex>

#include "SapoThreads.h"
#endif // WITH_THREADS

/**
 * Constructor that instantiates Sapo
 *
 * @param[in] model model to analyize
 */
Sapo::Sapo(Model *model):
    tmode(Bundle::AFO), decomp(0), max_param_splits(0), num_of_presplits(0),
    max_bundle_magnitude(std::numeric_limits<double>::max()),
    dyns(model->getDyns()), vars(model->getVars()), params(model->getParams())
{
}

/**
 * Reachable set computation
 *
 * @param[in] initSet bundle representing the current reached set
 * @param[in] k time horizon
 * @param[in,out] accounter acccounts for the computation progress
 * @returns the reached flowpipe
 */
Flowpipe Sapo::reach(const Bundle &initSet, unsigned int k,
                     ProgressAccounter *accounter) const
{
  // create current bundles list
  std::list<Bundle> cbundles = initSet.split(max_bundle_magnitude, 1.0);

  // create next bundles list
  std::list<Bundle> nbundles;

  // last polytope union in flowpipe
  PolytopesUnion last_step;
  last_step.add(initSet);

  // create flowpipe
  Flowpipe flowpipe(initSet.get_directions());
  flowpipe.append(initSet);

#ifdef WITH_THREADS
  std::mutex mutex;

  auto compute_next_bundles_and_add_to_last
      = [&nbundles, &last_step, &mutex](const Sapo *sapo, const Bundle &bundle)
#else
  auto compute_next_bundles_and_add_to_last
      = [&nbundles, &last_step](const Sapo *sapo, const Bundle &bundle)
#endif

  {
    // get the transformed bundle
    Bundle nbundle = bundle.transform(sapo->vars, sapo->dyns,
                                      sapo->tmode); // transform it

    if (sapo->decomp > 0) { // if requested, decompose it
      nbundle = nbundle.decompose(sapo->decomp_weight, sapo->decomp);
    }

    Polytope bls = nbundle;

    // TODO: check whether there is any chance for a transformed bundle to
    // be empty
    if (!bls.is_empty()) {
      // split if necessary the new reached bundle and add the resulting
      // bundles to the nbundles list
      nbundles.splice(nbundles.end(),
                      nbundle.split(sapo->max_bundle_magnitude));

      {
#ifdef WITH_THREADS
        std::unique_lock<std::mutex> lock(mutex);
#endif // WITH_THREADS
        last_step.add(bls);
      }
    }
  };

  unsigned int i = 0;

  // while time horizon has not been reached and last step reached is not empty
  // TODO: check whether there exists any chance for the last_step to be empty
  while (i < k && last_step.size() != 0) {

    // create a new last step reach set
    last_step = PolytopesUnion();
    i++;

#ifdef WITH_THREADS
    ThreadPool::BatchId batch_id = thread_pool.create_batch();

    // for all the old bundles
    for (auto b_it = std::cbegin(cbundles); b_it != std::cend(cbundles);
         ++b_it) {
      // submit the task to the thread pool
      thread_pool.submit_to_batch(batch_id,
                                  compute_next_bundles_and_add_to_last, this,
                                  std::ref(*b_it));
    }

    // join to the pool threads
    thread_pool.join_threads(batch_id);

    // close the batch
    thread_pool.close_batch(batch_id);
#else  // WITH_THREADS

    // for all the old bundles
    for (auto b_it = std::cbegin(cbundles); b_it != std::cend(cbundles);
         ++b_it) {

      compute_next_bundles_and_add_to_last(this, std::ref(*b_it));
    }
#endif // WITH_THREADS

    // swap current bundles and new bundles
    std::swap(cbundles, nbundles);

    nbundles = std::list<Bundle>();

    // add the last step to the flow pipe
    flowpipe.append(last_step); // store result

    if (accounter != NULL) {
      accounter->increase_performed();
    }
  }

  return flowpipe;
}

/**
 * Reachable set computation for parameteric dynamical systems
 *
 * @param[in] initSet bundle representing the current reached set
 * @param[in] pSet the set of parameters
 * @param[in] k time horizon
 * @param[in,out] accounter acccounts for the computation progress
 * @returns the reached flowpipe
 */
Flowpipe Sapo::reach(const Bundle &initSet, const PolytopesUnion &pSet,
                     unsigned int k, ProgressAccounter *accounter) const
{
  using namespace std;
  const unsigned int num_p_poly = pSet.size();

  // Each polytope in pSet corresponds to the list of bundles
  // reachable by using the parameters in that polytope.
  // This list is stored in an element of the vector cbundles

  // create current bundle list vector
  std::list<Bundle> cbundle_list = initSet.split(max_bundle_magnitude, 1.0);
  std::vector<std::list<Bundle>> cbundles(num_p_poly, cbundle_list);

  // create next bundles list
  std::vector<std::list<Bundle>> nbundles;

  // last polytope union in flowpipe
  PolytopesUnion last_step;
  last_step.add(initSet);

  // create flowpipe
  Flowpipe flowpipe(initSet.get_directions());
  flowpipe.append(initSet);

#ifdef WITH_THREADS_TEMP_DISABLED
  std::mutex add_mtx;

  auto compute_next_bundles_and_add_to_last
      = [&nbundles, &cbundles, &last_step, &add_mtx](
            const Sapo *sapo, const Polytope &pSet, const unsigned int pos)
#else
  auto compute_next_bundles_and_add_to_last
      = [&nbundles, &cbundles, &last_step](
            const Sapo *sapo, const Polytope &pSet, const unsigned int pos)
#endif

  {
    // for all the parameter sets
    for (auto b_it = std::cbegin(cbundles[pos]);
         b_it != std::cend(cbundles[pos]); ++b_it) {
      // get the transformed bundle
      Bundle nbundle
          = b_it->transform(sapo->vars, sapo->params, sapo->dyns, pSet,
                            sapo->tmode); // transform it

      if (sapo->decomp > 0) { // if requested, decompose it
        nbundle = nbundle.decompose(sapo->decomp_weight, sapo->decomp);
      }

      Polytope bls = nbundle;

      // TODO: check whether there is any chance for a transformed bundle to
      // be empty
      if (!bls.is_empty()) {
        // split if necessary the new reached bundle and add the resulting
        // bundles to the nbundles list
        nbundles[pos].splice(nbundles[pos].end(),
                             nbundle.split(sapo->max_bundle_magnitude));

        {
#ifdef WITH_THREADS_TEMP_DISABLED
          std::unique_lock<std::mutex> lock(add_mtx);
#endif
          last_step.add(bls);
        }
      }
    }
  };

  unsigned int i = 0;

  // while time horizon has not been reached and last step reached is not empty
  // TODO: check whether there exists any chance for the last_step to be empty
  while (i < k && last_step.size() != 0) {

    nbundles = std::vector<std::list<Bundle>>(num_p_poly);

    // create a new last step reach set
    last_step = PolytopesUnion();
    i++;

    unsigned int pSet_idx = 0;
#ifdef WITH_THREADS_TEMP_DISABLED
    ThreadPool::BatchId batch_id = thread_pool.create_batch();

    // for all the old bundles
    for (auto p_it = std::cbegin(pSet); p_it != std::cend(pSet); ++p_it) {
      // submit the task to the thread pool
      thread_pool.submit_to_batch(batch_id,
                                  compute_next_bundles_and_add_to_last, this,
                                  std::ref(*p_it), pSet_idx++);
    }

    // join to the pool threads
    thread_pool.join_threads(batch_id);

    // close the batch
    thread_pool.close_batch(batch_id);
#else  // WITH_THREADS

    // for all the old bundles
    for (auto p_it = std::cbegin(pSet); p_it != std::cend(pSet); ++p_it) {

      compute_next_bundles_and_add_to_last(this, std::ref(*p_it), pSet_idx++);
    }
#endif // WITH_THREADS

    // move the new bundles content in the current bundles
    cbundles = std::move(nbundles);

    // add the last step to the flow pipe
    flowpipe.append(last_step); // store result

    if (accounter != NULL) {
      accounter->increase_performed();
    }
  }

  if (accounter != NULL) {
    accounter->increase_performed_to(k);
  }

  return flowpipe;
}

std::list<PolytopesUnion>
get_a_finer_covering(const std::list<PolytopesUnion> &orig,
                     const unsigned int num_of_polytope_splits
                     = std::numeric_limits<unsigned>::max())
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

      std::list<Polytope> f_cov
          = (ps_it->begin())->split(num_of_polytope_splits);
      for (auto ls_it = std::begin(f_cov); ls_it != std::end(f_cov); ++ls_it) {
        result.push_back(*ls_it);
      }
      break;
    }
    case 2:  // the polytopes union contains more than one polytope
    default: // then, unpack them
    {
      for (auto ls_it = ps_it->begin(); ls_it != ps_it->end(); ++ls_it) {
        result.push_back(*ls_it);
      }
      break;
    }
    }
  }

  return result;
}

#ifdef WITH_THREADS
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
    std::unique_lock<std::shared_timed_mutex> writelock(mutex);

    list.push_back(obj);

    return *this;
  }

  ThreadSafeList<T> &push_back(const T &obj)
  {
    std::unique_lock<std::shared_timed_mutex> writelock(mutex);

    list.push_back(obj);

    return *this;
  }

  const std::list<T> &get_list() const
  {
    std::shared_lock<std::shared_timed_mutex> readlock(mutex);
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
synthesize_list(const Sapo &sapo, const Bundle &reachSet,
                const std::list<PolytopesUnion> &pSetList,
                const std::shared_ptr<STL::STL> &formula,
                ProgressAccounter *accounter)
{

#ifdef WITH_THREADS
  std::vector<PolytopesUnion> vect_res(pSetList.size());

  auto synthesize_funct
      = [&vect_res, &sapo, &reachSet, &formula,
         &accounter](const PolytopesUnion pSet, const unsigned int idx) {
          vect_res[idx] = sapo.synthesize(reachSet, pSet, formula);

          if (accounter != NULL) {
            accounter->increase_performed(formula->time_bounds().end());
          }
        };

  ThreadPool::BatchId batch_id = thread_pool.create_batch();

  unsigned int res_idx = 0;
  for (auto ps_it = std::begin(pSetList); ps_it != std::end(pSetList);
       ++ps_it) {
    // submit the task to the thread pool
    thread_pool.submit_to_batch(batch_id, synthesize_funct, *ps_it, res_idx++);
  }

  // join to the pool threads
  thread_pool.join_threads(batch_id);

  // close the batch
  thread_pool.close_batch(batch_id);

  return std::list<PolytopesUnion>(std::make_move_iterator(vect_res.begin()),
                                   std::make_move_iterator(vect_res.end()));

#else  // WITH_THREADS
  std::list<PolytopesUnion> results;

  for (auto ps_it = std::begin(pSetList); ps_it != std::end(pSetList);
       ++ps_it) {
    results.push_back(sapo.synthesize(reachSet, *ps_it, formula));
    if (accounter != NULL) {
      accounter->increase_performed(formula->time_bounds().end());
    }
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
 * @param[in] num_of_presplits is number of splits to be performed before
 *                             the computation
 * @param[in,out] accounter acccounts for the computation progress
 * @returns the list of refined parameter sets
 */
std::list<PolytopesUnion> Sapo::synthesize(const Bundle &reachSet,
                                           const PolytopesUnion &pSet,
                                           const std::shared_ptr<STL::STL> formula,
                                           const unsigned int max_splits,
                                           const unsigned int num_of_presplits,
                                           ProgressAccounter *accounter) const
{
  std::list<PolytopesUnion> pSetList{pSet};

  if (num_of_presplits > 1) {
    pSetList = get_a_finer_covering(pSetList, num_of_presplits);
  }

  const unsigned int max_time = formula->time_bounds().end();
  unsigned int already_performed_steps = 0;

  unsigned int num_of_splits = 0;
  std::list<PolytopesUnion> res
      = synthesize_list(*this, reachSet, pSetList, formula, accounter);

  if (accounter) {
    already_performed_steps = max_time * pSetList.size();
    accounter->increase_performed_to(already_performed_steps);
  }

  while (every_set_is_empty(res) && num_of_splits++ < max_splits) {
    pSetList = get_a_finer_covering(pSetList);

    res = synthesize_list(*this, reachSet, pSetList, formula, accounter);

    if (accounter != NULL) {
      already_performed_steps += max_time * pSetList.size();
      accounter->increase_performed_to(already_performed_steps);
    }
  }

  for (auto lss_it = std::begin(res); lss_it != std::end(res); ++lss_it) {
    lss_it->simplify();
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
                                const std::shared_ptr<STL::Conjunction> conj) const
{
  PolytopesUnion Pu1
      = this->synthesize(reachSet, pSet, conj->get_left_subformula());
  PolytopesUnion Pu2
      = this->synthesize(reachSet, pSet, conj->get_right_subformula());
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
                                const std::shared_ptr<STL::Disjunction> disj) const
{
  PolytopesUnion Pu
      = this->synthesize(reachSet, pSet, disj->get_left_subformula());
  Pu.add(this->synthesize(reachSet, pSet, disj->get_right_subformula()));

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
                                const std::shared_ptr<STL::Eventually> ev) const
{
  std::shared_ptr<STL::Atom> true_atom = std::make_shared<STL::Atom>(-1);

  std::shared_ptr<STL::Until> u
      = std::make_shared<STL::Until>(true_atom, ev->time_bounds().begin(),
                                     ev->time_bounds().end(),
                                     ev->get_subformula());

  return this->synthesize(reachSet, pSet, u, 0);
}

/**
 * Parameter synthesis method
 *
 * @param[in] reachSet bundle representing the current reached set
 * @param[in] pSet the current parameter set
 * @param[in] formula is an STL specification for the model
 * @param[in,out] accounter acccounts for the computation progress
 * @returns a parameter set refined according with `formula`
 */
PolytopesUnion Sapo::synthesize(const Bundle &reachSet,
                                const PolytopesUnion &pSet,
                                const std::shared_ptr<STL::STL> formula,
                                ProgressAccounter *accounter) const
{
  (void)accounter;

  switch (formula->get_type()) {

  // Atomic predicate
  case STL::ATOM:
    return synthesize(reachSet, pSet,
                      std::dynamic_pointer_cast<STL::Atom>(formula));

  // Conjunction
  case STL::CONJUNCTION:
    return synthesize(reachSet, pSet,
                      std::dynamic_pointer_cast<STL::Conjunction>(formula));

  // Disjunction
  case STL::DISJUNCTION:
    return synthesize(reachSet, pSet,
                      std::dynamic_pointer_cast<STL::Disjunction>(formula));

  // Until
  case STL::UNTIL:
    return synthesize(reachSet, pSet,
                      std::dynamic_pointer_cast<STL::Until>(formula), 0);

  // Always
  case STL::ALWAYS:
    return synthesize(reachSet, pSet,
                      std::dynamic_pointer_cast<STL::Always>(formula), 0);

  // Eventually
  case STL::EVENTUALLY:
    return synthesize(reachSet, pSet,
                      std::dynamic_pointer_cast<STL::Eventually>(formula));

  default:
    throw std::logic_error("Unsupported formula");
  }
}

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
                                const std::shared_ptr<STL::Atom> atom) const
{
  return reachSet.synthesize(vars, params, dyns, pSet, atom);
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
                                const std::shared_ptr<STL::Until> formula,
                                const int time) const
{
  const TimeInterval &t_itvl = formula->time_bounds();

  // Base case
  if (t_itvl.is_empty())
    return PolytopesUnion();

  // Until interval far
  if (t_itvl > time) {
    // Synthesize wrt phi1
    PolytopesUnion P1
        = this->synthesize(reachSet, pSet, formula->get_left_subformula());
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
        = this->synthesize(reachSet, pSet, formula->get_left_subformula());

    if (P1.is_empty()) {
      return this->synthesize(reachSet, pSet, formula->get_right_subformula());
    }

    PolytopesUnion result
        = transition_and_synthesis(reachSet, P1, formula, time);
    result.add(
        this->synthesize(reachSet, pSet, formula->get_right_subformula()));

    return result;
  }

  // If none of the above condition holds, then it must holds that :
  // 			t_itvl.begin()<=time and t_itvl.end()==time
  return this->synthesize(reachSet, pSet, formula->get_right_subformula());
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
                                const std::shared_ptr<STL::Always> formula,
                                const int time) const
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
    //		std::cout << "Inside interval (a = 0)" << std::endl;

    // Refine wrt phi
    PolytopesUnion P
        = this->synthesize(reachSet, pSet, formula->get_subformula());

    if (P.is_empty()) {
      return P;
    }

    return transition_and_synthesis(reachSet, P, formula, time);
  }

  // If none of the above condition holds, then it must holds that :
  // 			t_itvl.begin()<=time and t_itvl.end()==time
  return this->synthesize(reachSet, pSet, formula->get_subformula());
}
