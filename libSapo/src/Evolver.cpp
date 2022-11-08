/**
 * @file Evolver.cpp
 * @author Alberto Casagrande <acasagrande@units.it>
 * @brief Transform geometric sets subject to a dynamical system
 * @version 0.1
 * @date 2022-07-02
 *
 * @copyright Copyright (c) 2022
 */

#include "Evolver.h"

#include <utility>

#ifdef WITH_THREADS
#include <mutex>
#include <shared_mutex>

#include "SapoThreads.h"
#endif // WITH_THREADS

#include "Bernstein.h"
#include "VarsGenerator.h"

/**
 * @brief Avoid \f$-0\f$
 *
 * This macro avoids \f$-0\f$ by replacing it by \f$0\f$.
 */
#define AVOID_NEG_ZERO(value) ((value) == 0 ? 0 : (value))

/**
 * @brief Replace some variables in expression numerators
 *
 * This function replaces all the occurrences of the variable `vars[i]`
 * in `expressions` numerators with `subs[i]`.
 *
 * @param expressions is the vector of expressions whose variable
 *                    occurrences must be replaced
 * @param vars is the vector variables whose occurrences must be
 *             replaced
 * @param subs is the vector of expressions that must replace
 *             variable occurrences
 * @return A vector of the containing `expressions` numerators
 *         in which each occurrence of the variables in `vars`
 *         have been replaced by the expressions in `subs`
 */
template<typename T>
std::pair<std::vector<SymbolicAlgebra::Expression<T>>,
          std::vector<SymbolicAlgebra::Expression<T>>>
replace_in_rational(
    const std::vector<SymbolicAlgebra::Expression<T>> &expressions,
    const std::vector<SymbolicAlgebra::Symbol<T>> &vars,
    const std::vector<SymbolicAlgebra::Expression<T>> &subs)
{
  using namespace SymbolicAlgebra;

  Expression<>::replacement_type repl;

  for (unsigned int k = 0; k < vars.size(); ++k) {
    repl[vars[k]] = subs[k];
  }

  std::pair<std::vector<Expression<T>>, std::vector<Expression<T>>> results;
  for (auto ex_it = std::begin(expressions); ex_it != std::end(expressions);
       ++ex_it) {
    results.first.push_back(ex_it->get_numerator().replace(repl));
    results.second.push_back(ex_it->get_denominator().replace(repl));
  }

  return results;
}

/**
 * @brief Replace some variables in expressions by other expressions
 *
 * This function replaces all the occurrences of the variable `vars[i]`
 * in any expression in `expressions` with `subs[i]`.
 *
 * @param expressions is the vector of expressions whose variable
 *                    occurrences must be replaced
 * @param vars is the vector variables whose occurrences must be
 *             replaced
 * @param subs is the vector of expressions that must replace
 *             variable occurrences
 * @return A vector of the containing the input expression
 *         in which each occurrence of the variables in `vars`
 *         have been replaced by the expressions in `subs`
 */
std::vector<SymbolicAlgebra::Expression<>>
replace_in(const std::vector<SymbolicAlgebra::Expression<>> &expressions,
           const std::vector<SymbolicAlgebra::Symbol<>> &vars,
           const std::vector<SymbolicAlgebra::Expression<>> &subs)
{
  using namespace SymbolicAlgebra;

  Expression<>::replacement_type repl;

  for (unsigned int k = 0; k < vars.size(); ++k) {
    repl[vars[k]] = subs[k];
  }

  std::vector<Expression<>> results;
  for (auto ex_it = std::begin(expressions); ex_it != std::end(expressions);
       ++ex_it) {
    results.push_back(Expression<>(*ex_it).replace(repl));
  }

  return results;
}

/**
 * @brief Compute the Bernstein coefficients for a direction
 *
 * @param alpha is the vector of the variables
 * @param f is the function whose Bernstein coefficients must be computed
 * @param direction is the direction along
 * @return the Bernstein coefficients for `f` on `direction`
 */
std::vector<SymbolicAlgebra::Expression<>>
compute_Bern_coefficients(const std::vector<SymbolicAlgebra::Symbol<>> &alpha,
                          const std::vector<SymbolicAlgebra::Expression<>> &f,
                          const LinearAlgebra::Vector<double> &direction)
{
  SymbolicAlgebra::Expression<> Lfog = 0;
  // upper facets
  for (unsigned int k = 0; k < direction.size(); k++) {
    if (direction[k] != 0) {
      Lfog += direction[k] * f[k];
    }
  }

  return get_Bernstein_coefficients(alpha, Lfog);
}

/**
 * @brief Compute the variable substitutions for a parallelotope
 *
 * @param P is a parallelotope.
 * @param q are the variables associated to the parallelotope's base vertex.
 * @param beta are the variables associated to the parallelotope's lengths.
 * @return the symbolic equations representing the the variable
 *         substitutions for `P`.
 */
SymbolicAlgebra::Expression<>::replacement_type
get_subs_from(const Parallelotope &P,
              const std::vector<SymbolicAlgebra::Symbol<>> &q,
              const std::vector<SymbolicAlgebra::Symbol<>> &beta)
{
  using namespace SymbolicAlgebra;

  const LinearAlgebra::Vector<double> &base_vertex = P.base_vertex();
  const LinearAlgebra::Vector<double> &lengths = P.lengths();

  Expression<>::replacement_type repl;

  for (unsigned int k = 0; k < q.size(); k++) {
    repl[q[k]] = base_vertex[k];
    repl[beta[k]] = lengths[k];
  }

  return repl;
}

/**
 * @brief Build the generator function of a parallelotope
 *
 * This function build the generator functions of a parallelotope.
 * In particular, it returns the symbolic vector:
 * \f[\textrm{base} + ((\textrm{alpha} \circ \lambda)^T \cdot D)^T\f]
 * where \f$\circ\f$ is the Hadamard product, and
 * \f$\textrm{base}\f$, \f$\lambda\f$, and $D$ are the base vertex,
 * the generator lengths vector, the generators matrix
 * of the parallelotope \f$P\f$, respectively.
 *
 * @param alpha is a vector of variables
 * @param P is the considered parallelotope
 * @return The generator function of `P`
 */
std::vector<SymbolicAlgebra::Expression<>> build_instantiated_generator_functs(
    const std::vector<SymbolicAlgebra::Symbol<>> &alpha,
    const Parallelotope &P)
{
  using namespace LinearAlgebra;

  std::vector<SymbolicAlgebra::Expression<>> gen_functs;
  for (auto it = std::begin(P.base_vertex()); it != std::end(P.base_vertex());
       ++it) {
    gen_functs.push_back(*it);
  }

  const std::vector<Vector<double>> &generators = P.generators();

  for (unsigned int i = 0; i < generators.size(); i++) {
    // some of the non-null rows of the generator matrix
    // correspond to 0-length dimensions in degenerate
    // parallelotopes and must be avoided
    if (P.lengths()[i] != 0) {
      Vector<double> vector = P.lengths()[i] * generators[i];
      for (unsigned int j = 0; j < vector.size(); j++) {
        if (vector[j] != 0) {
          gen_functs[j] += alpha[i] * vector[j];
        }
      }
    }
  }

  return gen_functs;
}

/**
 * @brief Build the generator function of a symbolic parallelotope
 *
 * This function build the generator functions of a symbolic
 * parallelotope. In particular, it returns the symbolic vector:
 * \f[\textrm{base} + ((\textrm{alpha} \circ \lambda)^T \cdot D)^T\f]
 * where \f$\circ\f$ is the Hadamard product, and
 * \f$\textrm{base}f$, \f$\lambda\f$, and \f$D\f$ are the base vertex,
 * the generator lengths vector, the generators matrix
 * of a generic parallelotope, respectively.
 *
 * @param base is the base vertex variable vector
 * @param alpha is a vector of variables
 * @param lambda is the generator length variable vector
 * @return The generator function of a parallelotope having `base`
 *       as base vertex,`lambda` as generator length vector, and
 *       `D` as generator matrix
 */
template<typename T>
std::vector<SymbolicAlgebra::Expression<T>>
build_generator_functs(const std::vector<SymbolicAlgebra::Symbol<T>> &base,
                       const std::vector<SymbolicAlgebra::Symbol<T>> &alpha,
                       const std::vector<SymbolicAlgebra::Symbol<T>> &lambda,
                       const LinearAlgebra::Dense::Matrix<T> &D)
{
  using namespace LinearAlgebra;

  std::vector<SymbolicAlgebra::Expression<>> gen_functs;
  gen_functs.reserve(base.size());
  for (auto it = std::begin(base); it != std::end(base); ++it) {
    gen_functs.push_back(*it);
  }

  for (unsigned int i = 0; i < D.size(); i++) {
    const std::vector<T> &generator = D[i];
    for (unsigned int j = 0; j < generator.size(); j++) {
      if (generator[j] != 0) {
        gen_functs[j] += alpha[i] * lambda[i] * generator[j];
      }
    }
  }

  return gen_functs;
}

/**
 * @brief A container whose accesses are synchronized
 *
 * The objects of this class store values that can be
 * accessed in a synchronized way according to a mutex.
 * The value themselves are updated exclusively if the
 * call `COND::operator()` on the stored value and the
 * possible new value, returns `true`.
 *
 * @tparam COND
 */
template<typename T, typename COND>
class CondSyncUpdater
{
#ifdef WITH_THREADS
  mutable std::shared_timed_mutex
      mutex; //!< The mutex that synchronizes the accesses
#endif
  T _value;  //!< the value stored in the object
  COND _cmp; //!< the condition that must be satisfied to update the value

public:
  /**
   * @brief Construct a new Cond Sync Update object
   */
  CondSyncUpdater()
  {
    // initialize _value to the minimum wrt the order
    _value = (_cmp(-std::numeric_limits<T>::infinity(),
                   std::numeric_limits<T>::infinity())
                  ? std::numeric_limits<T>::infinity()
                  : -std::numeric_limits<T>::infinity());
  }

  /**
   * @brief Get the value store in the updater
   *
   * @return the value stored in this object
   */
  inline operator T() const
  {
#ifdef WITH_THREADS
    std::shared_lock<std::shared_timed_mutex> readlock(mutex);
#endif
    return _value;
  }

  /**
   * @brief Update the value if the condition holds
   *
   * This method updates the value stored by the updater
   * if and only if the call
   * `COND::operator()(_value, value)` returns `true`.
   * @param value is the possible new value contained
   *        by the updater
   */
  void update(const T &value)
  {

#ifdef WITH_THREADS
    std::unique_lock<std::shared_timed_mutex> writelock(mutex);
#endif

    if (_cmp(value, _value)) {
      _value = value;
    }
  }
};

/**
 * @brief A class to find the minimum and the maximum Bernstein coefficients
 */
class MinMaxCoeffFinder
{
  /**
   * @brief Evaluate a Bernstein coefficient
   *
   * @param bernCoeff is the symbolical representation of Bernstein
   *    coefficient
   * @return The numerical evaluation of Bernstein coefficient
   */
  double eval_coeff(const SymbolicAlgebra::Expression<> &bernCoeff) const;

public:
  /**
   * @brief Create MaxCoeffType objects.
   */
  MinMaxCoeffFinder() {}

  /**
   * @brief Find the interval containing the Bernstein coefficients.
   *
   * @param b_coeffs is a list of symbolical Bernstein coefficients.
   * @return The pair minimum-maximum among all the Bernstein
   *          coefficients in `b_coeffs`.
   */
  virtual std::pair<double, double>
  operator()(const std::vector<SymbolicAlgebra::Expression<>> &b_coeffs) const;

  virtual ~MinMaxCoeffFinder() {}
};

/**
 * @brief  A class to find the minimum and the maximum parametric Bernstein
 * coefficients
 */
class ParamMinMaxCoeffFinder : public MinMaxCoeffFinder
{
  const std::vector<SymbolicAlgebra::Symbol<>> &params;
  const Polytope &paraSet;
  /**
   * @brief Evaluate the parametric Bernstein coefficient upper-bound
   *
   * @param bernCoeff is the symbolical representation of Bernstein
   *                  coefficient
   * @return The numerical evaluation of parametric Bernstein
   *         coefficient upper-bound
   */
  double maximize_coeff(const SymbolicAlgebra::Expression<> &bernCoeff) const;

  /**
   * @brief Evaluate the parametric Bernstein coefficient lower-bound
   *
   * @param bernCoeff is the symbolical representation of Bernstein
   *                  coefficient
   * @return The numerical evaluation of parametric Bernstein coefficient
   *         lower-bound
   */
  double minimize_coeff(const SymbolicAlgebra::Expression<> &bernCoeff) const;

public:
  /**
   * @brief Constructor
   *
   * @param params is the list of parameters
   * @param paraSet is the set of admissible values for parameters
   */
  ParamMinMaxCoeffFinder(const std::vector<SymbolicAlgebra::Symbol<>> &params,
                         const Polytope &paraSet):
      MinMaxCoeffFinder(),
      params(params), paraSet(paraSet)
  {
  }

  /**
   * @brief Find the interval containing the Bernstein coefficients
   *
   * @param b_coeffs is a list of symbolical Bernstein coefficients
   * @return The pair minimum-maximum among all the Bernstein
   *          coefficients in `b_coeffs`
   */
  std::pair<double, double>
  operator()(const std::vector<SymbolicAlgebra::Expression<>> &b_coeffs) const;

  ~ParamMinMaxCoeffFinder() {}
};

double MinMaxCoeffFinder::eval_coeff(
    const SymbolicAlgebra::Expression<> &bernCoeff) const
{
  double value = bernCoeff.evaluate();

  // TODO: The following conditional evaluation avoids -0
  //       values. Check the difference between -0 and 0.
  return AVOID_NEG_ZERO(value);
}

double ParamMinMaxCoeffFinder::maximize_coeff(
    const SymbolicAlgebra::Expression<> &bernCoeff) const
{
  return paraSet.maximize(params, bernCoeff).optimum();
}

double ParamMinMaxCoeffFinder::minimize_coeff(
    const SymbolicAlgebra::Expression<> &bernCoeff) const
{
  return paraSet.minimize(params, bernCoeff).optimum();
}

std::pair<double, double> MinMaxCoeffFinder::operator()(
    const std::vector<SymbolicAlgebra::Expression<>> &b_coeffs) const
{
  auto b_coeff_it = b_coeffs.begin();

  double max_value = eval_coeff(*b_coeff_it);
  double min_value = max_value;

  for (++b_coeff_it; b_coeff_it != b_coeffs.end(); ++b_coeff_it) {
    double coeff_value = eval_coeff(*b_coeff_it);

    if (coeff_value > max_value) {
      max_value = coeff_value;
    }

    if (coeff_value < min_value) {
      min_value = coeff_value;
    }
  }

  return std::pair<double, double>(min_value, max_value);
}

std::pair<double, double> ParamMinMaxCoeffFinder::operator()(
    const std::vector<SymbolicAlgebra::Expression<>> &b_coeffs) const
{
  auto b_coeff_it = b_coeffs.begin();

  double max_value = maximize_coeff(*b_coeff_it);
  double min_value = minimize_coeff(*b_coeff_it);

  for (++b_coeff_it; b_coeff_it != b_coeffs.end(); ++b_coeff_it) {
    double coeff_value = maximize_coeff(*b_coeff_it);

    if (coeff_value > max_value) {
      max_value = coeff_value;
    }

    coeff_value = minimize_coeff(*b_coeff_it);
    if (coeff_value < min_value) {
      min_value = coeff_value;
    }
  }

  return std::pair<double, double>(min_value, max_value);
}

template<typename T>
inline typename SymbolicAlgebra::Expression<T>::interpretation_type
make_interpretation(const std::vector<SymbolicAlgebra::Symbol<T>> &symbols,
                    const std::vector<T> &values)
{
  if (symbols.size() != values.size()) {
    throw std::runtime_error("make_interpretation: symbols and values "
                             "sizes are not the same");
  }

  typename SymbolicAlgebra::Expression<T>::interpretation_type inter;
  for (size_t i = 0; i < symbols.size(); ++i) {
    inter[symbols[i]] = values[i];
  }

  return inter;
}

LinearAlgebra::Vector<double> get_approx_center(const Bundle &bundle)
{
  using namespace LinearAlgebra;
  std::vector<double> approx_c;

  if (bundle.dim() == 0) {
    return approx_c;
  }

  for (const auto &b_template: bundle.templates()) {
    if (approx_c.size() == 0) {
      approx_c = (bundle.get_parallelotope(b_template)).center();
    } else {
      approx_c = approx_c + (bundle.get_parallelotope(b_template)).center();
    }
  }

  return approx_c / (double)bundle.num_of_templates();
}

std::vector<SymbolicAlgebra::Expression<double>>
average_dynamics(const DynamicalSystem<double> &ds,
                 const Polytope &parameter_set)
{
  using namespace SymbolicAlgebra;

  std::vector<double> approx_c = get_approx_center(Bundle(parameter_set));

  auto inter{make_interpretation(ds.parameters(), approx_c)};

  std::vector<SymbolicAlgebra::Expression<>> avg_ds;
  for (const auto &dyn: ds.dynamics()) {
    avg_ds.emplace_back(dyn.apply(inter));
  }

  return avg_ds;
}

/**
 * @brief Evaluate a vector of expressions over an interpretation
 *
 * @tparam T is the type of the expression constants
 * @param expressions is the vector of expressions to be evaluated
 * @param inter is the interpretation to be applied
 * @return the vector of evaluations of the expressions in the
 *         input vector
 */
template<typename T>
LinearAlgebra::Vector<T> evaluate(
    const LinearAlgebra::Vector<SymbolicAlgebra::Expression<T>> &expressions,
    const typename SymbolicAlgebra::Expression<T>::interpretation_type &inter)
{
  LinearAlgebra::Vector<T> result;

  for (const auto &e: expressions) {
    result.emplace_back(e.apply(inter).evaluate());
  }

  return result;
}

std::vector<LinearAlgebra::Vector<double>>
get_new_directions(const Bundle &bundle, const DynamicalSystem<double> &ds,
                   const Polytope &parameter_set)
{
  std::vector<LinearAlgebra::Vector<double>> A;

  auto avg_ds = average_dynamics(ds, parameter_set);

  auto approx_c = get_approx_center(bundle);

  const auto lengths = bundle.edge_lengths();
  for (size_t i = 0; i < bundle.size(); ++i) {
    using namespace LinearAlgebra;

    const auto &dir = bundle.get_direction(i);
    Vector<double> delta = (lengths[i] * dir) / norm_2(dir);

    auto inter = make_interpretation(ds.variables(), approx_c + delta);

    const auto up{evaluate(ds.dynamics(), inter)};

    inter = make_interpretation(ds.variables(), approx_c - delta);

    const auto down{evaluate(ds.dynamics(), inter)};

    const auto new_dir{up - down};

    A.push_back(new_dir);
  }

  return A;
}

template<>
Bundle Evolver<double>::operator()(const Bundle &bundle,
                                   const Polytope &parameter_set)
{
  using namespace std;
  using namespace SymbolicAlgebra;
  using namespace LinearAlgebra;

  if (bundle.dim() != _ds.variables().size()) {
    throw std::domain_error("The bundle and the dynamic laws must have the "
                            "same number of dimensions.");
  }

  if (parameter_set.dim() != _ds.parameters().size()) {
    throw std::domain_error("The vector of parameters and the parameter "
                            "set must have same number of dimensions.");
  }

  std::vector<LinearAlgebra::Vector<double>> new_dirs;
  if (this->dynamic_directions) { // if dynamic directions is enabled

    // compute new directions
    new_dirs = get_new_directions(bundle, _ds, parameter_set);
  } else { // if dynamic directions is *not* enabled

    // use old directions
    new_dirs = bundle.directions();
  }

  vector<CondSyncUpdater<double, std::less<double>>> max_coeffs(bundle.size());
  vector<CondSyncUpdater<double, std::greater<double>>> min_coeffs(
      bundle.size());
  MinMaxCoeffFinder *itvl_finder;

  if (_ds.parameters().size() == 0) {
    itvl_finder = new MinMaxCoeffFinder();
  } else {
    itvl_finder = new ParamMinMaxCoeffFinder(_ds.parameters(), parameter_set);
  }

  std::vector<Symbol<double>> alpha
      = get_symbol_vector<double>("alpha", _ds.dim());

  auto refine_coeff_itvl =
      [&max_coeffs, &min_coeffs, &bundle, &new_dirs, &alpha,
       &itvl_finder](const DynamicalSystem<double> &ds,
                     const std::vector<unsigned int> &bundle_template,
                     const evolver_mode t_mode) {
        Parallelotope P = bundle.get_parallelotope(bundle_template);

        const auto &genFun = build_instantiated_generator_functs(alpha, P);
        const auto genFun_f
            = replace_in(ds.dynamics(), ds.variables(), genFun);

        unsigned int dir_b;

        // for each direction
        const size_t num_of_dirs
            = (t_mode == ONE_FOR_ONE ? bundle_template.size() : bundle.size());

        for (unsigned int j = 0; j < num_of_dirs; j++) {
          dir_b = (t_mode == ONE_FOR_ONE ? bundle_template[j] : j);

          const auto coefficients
              = compute_Bern_coefficients(alpha, genFun_f, new_dirs[dir_b]);

          const auto coeff_itvl = (*itvl_finder)(coefficients);

          min_coeffs[dir_b].update(coeff_itvl.first);
          max_coeffs[dir_b].update(coeff_itvl.second);
        }
      };

  try {
#ifdef WITH_THREADS
    ThreadPool::BatchId batch_id = thread_pool.create_batch();

    for (auto t_it = std::begin(bundle.templates());
         t_it != std::end(bundle.templates()); ++t_it) {

      // submit the task to the thread pool
      thread_pool.submit_to_batch(batch_id, refine_coeff_itvl, std::ref(_ds),
                                  std::ref(*t_it), this->mode);
    }

    // join to the pool threads
    thread_pool.join_threads(batch_id);

    // close the batch
    thread_pool.close_batch(batch_id);
#else  // WITH_THREADS
    for (auto t_it = std::begin(bundle.templates());
         t_it != std::end(bundle.templates()); ++t_it) {
      refine_coeff_itvl(std::ref(_ds), std::ref(*t_it), _mode);
    }
#endif // WITH_THREADS
    delete itvl_finder;
  } catch (SymbolicAlgebra::symbol_evaluation_error &e) {

    delete itvl_finder;
    std::ostringstream oss;

    oss << "The symbol \"" << e.get_symbol_name() << "\" is unknown";
    throw std::domain_error(oss.str());
  }

  LinearAlgebra::Vector<double> lower_bounds, upper_bounds;
  for (auto it = std::begin(max_coeffs); it != std::end(max_coeffs); ++it) {
    upper_bounds.push_back(*it);
  }
  for (auto it = std::begin(min_coeffs); it != std::end(min_coeffs); ++it) {
    lower_bounds.push_back(*it);
  }

  Bundle res(std::move(new_dirs), std::move(lower_bounds),
             std::move(upper_bounds), bundle.templates());

  for (const auto &len: bundle.edge_lengths()) {
    if (len > _edge_threshold) {
      std::ostringstream ss;

      ss << "Evolver<double>::operator(): one of the "
            "computed bundle edge lengths is larger than "
            "the set threshold, i.e., "
         << _edge_threshold;

      throw std::runtime_error(ss.str());
    }
  }

  if (this->mode == ONE_FOR_ONE) {
    res.canonize();
  }

  return res;
}

template<typename T>
inline void init_genFun_if_necessary(
    std::pair<std::vector<SymbolicAlgebra::Expression<T>>,
              std::vector<SymbolicAlgebra::Expression<T>>> &genFun_f,
    const Parallelotope &P, const DynamicalSystem<T> &ds,
    const std::vector<SymbolicAlgebra::Symbol<T>> &base,
    const std::vector<SymbolicAlgebra::Symbol<T>> &alpha,
    const std::vector<SymbolicAlgebra::Symbol<T>> &lambda)
{
  if (genFun_f.first.size() == 0) { // not initialized yet
    auto genFun = build_generator_functs(base, alpha, lambda, P.generators());

    genFun_f = replace_in_rational(ds.dynamics(), ds.variables(), genFun);
  }
}

template<typename T>
const std::vector<SymbolicAlgebra::Expression<T>> &get_symbolic_coefficients(
    CachedEvolver<T> *evolver, const Parallelotope &P,
    std::pair<std::vector<SymbolicAlgebra::Expression<T>>,
              std::vector<SymbolicAlgebra::Expression<T>>> &genFun_f,
    const std::vector<SymbolicAlgebra::Symbol<T>> &base,
    const std::vector<SymbolicAlgebra::Symbol<T>> &alpha,
    const std::vector<SymbolicAlgebra::Symbol<T>> &lambda,
    const LinearAlgebra::Vector<double> &bundle_dir)
{
  if (evolver->coefficients_in_cache(P,
                                     bundle_dir)) { // Bernstein coefficients
                                                    // have been computed

    return evolver->get_cached_coefficients(P, bundle_dir);
  }

  init_genFun_if_necessary(genFun_f, P, evolver->dynamical_system(), base,
                           alpha, lambda);

  auto coeff = compute_Bern_coefficients(alpha, genFun_f, bundle_dir);

  return evolver->set_cache_coefficients(P, bundle_dir, std::move(coeff));
}

template<typename T>
inline void
init_genFun_if_necessary(std::vector<SymbolicAlgebra::Expression<T>> &genFun_f,
                         const Parallelotope &P, const DynamicalSystem<T> &ds,
                         const std::vector<SymbolicAlgebra::Symbol<T>> &base,
                         const std::vector<SymbolicAlgebra::Symbol<T>> &alpha,
                         const std::vector<SymbolicAlgebra::Symbol<T>> &lambda)
{
  if (genFun_f.size() == 0) { // not initialized yet
    auto genFun = build_generator_functs(base, alpha, lambda, P.generators());

    genFun_f = replace_in(ds.dynamics(), ds.variables(), genFun);
  }
}

template<typename T>
const std::vector<SymbolicAlgebra::Expression<T>> &get_symbolic_coefficients(
    CachedEvolver<T> *evolver, const Parallelotope &P,
    std::vector<SymbolicAlgebra::Expression<T>> &genFun_f,
    const std::vector<SymbolicAlgebra::Symbol<T>> &base,
    const std::vector<SymbolicAlgebra::Symbol<T>> &alpha,
    const std::vector<SymbolicAlgebra::Symbol<T>> &lambda,
    const LinearAlgebra::Vector<double> &bundle_dir)
{
  if (evolver->coefficients_in_cache(P,
                                     bundle_dir)) { // Bernstein coefficients
                                                    // have been computed

    return evolver->get_cached_coefficients(P, bundle_dir);
  }

  init_genFun_if_necessary(genFun_f, P, evolver->dynamical_system(), base,
                           alpha, lambda);

  return evolver->set_cache_coefficients(
      P, bundle_dir, compute_Bern_coefficients(alpha, genFun_f, bundle_dir));
}

template<>
Bundle CachedEvolver<double>::operator()(const Bundle &bundle,
                                         const Polytope &parameter_set)
{
  using namespace std;
  using namespace SymbolicAlgebra;
  using namespace LinearAlgebra;

  if (bundle.dim() != _ds.variables().size()) {
    throw std::domain_error("The bundle and the dynamic laws must have the "
                            "same number of dimensions.");
  }

  if (parameter_set.dim() != _ds.parameters().size()) {
    throw std::domain_error("The vector of parameters and the parameter "
                            "set must have same number of dimensions.");
  }

  std::vector<LinearAlgebra::Vector<double>> new_dirs;
  if (this->dynamic_directions) { // if dynamic directions is enabled

    // compute new directions
    new_dirs = get_new_directions(bundle, _ds, parameter_set);
  } else { // if dynamic directions is *not* enabled

    // use old directions
    new_dirs = bundle.directions();
  }

  vector<CondSyncUpdater<double, std::less<double>>> max_coeffs(bundle.size());
  vector<CondSyncUpdater<double, std::greater<double>>> min_coeffs(
      bundle.size());
  MinMaxCoeffFinder *itvl_finder;

  if (_ds.parameters().size() == 0) {
    itvl_finder = new MinMaxCoeffFinder();
  } else {
    itvl_finder = new ParamMinMaxCoeffFinder(_ds.parameters(), parameter_set);
  }

  auto alpha = get_symbol_vector<double>("alpha", _ds.dim());
  auto lambda = get_symbol_vector<double>("lambda", _ds.dim());
  auto base = get_symbol_vector<double>("base", _ds.dim());

  auto refine_coeff_itvl =
      [&max_coeffs, &min_coeffs, &bundle, &new_dirs, &alpha, &lambda, &base,
       &itvl_finder](CachedEvolver<double> *evolver,
                     const std::vector<unsigned int> &bundle_template) {
        Parallelotope P = bundle.get_parallelotope(bundle_template);

        Expression<double>::interpretation_type P_interpretation;

        for (size_t i = 0; i < P.dim(); ++i) {
          P_interpretation[base[i]] = P.base_vertex()[i];
          P_interpretation[lambda[i]] = P.lengths()[i];
        }

        std::vector<SymbolicAlgebra::Expression<double>> genFun_f;

        // for each direction
        const size_t num_of_dirs
            = (evolver->mode == ONE_FOR_ONE ? bundle_template.size()
                                            : bundle.size());

        for (unsigned int j = 0; j < num_of_dirs; j++) {
          auto dir_b = (evolver->mode == ONE_FOR_ONE ? bundle_template[j] : j);

          auto symb_coeff = get_symbolic_coefficients(
              evolver, P, genFun_f, base, alpha, lambda, new_dirs[dir_b]);

          std::vector<Expression<double>> coefficients;

          coefficients.reserve(symb_coeff.size());
          for (auto &coeff: symb_coeff) {
            coefficients.push_back(coeff.apply(P_interpretation));
          }

          auto coeff_itvl = (*itvl_finder)(coefficients);

          min_coeffs[dir_b].update(coeff_itvl.first);
          max_coeffs[dir_b].update(coeff_itvl.second);
        }
      };

  try {
#ifdef WITH_THREADS
    ThreadPool::BatchId batch_id = thread_pool.create_batch();

    for (auto t_it = std::begin(bundle.templates());
         t_it != std::end(bundle.templates()); ++t_it) {

      // submit the task to the thread pool
      thread_pool.submit_to_batch(batch_id, refine_coeff_itvl, this,
                                  std::ref(*t_it));
    }

    // join to the pool threads
    thread_pool.join_threads(batch_id);

    // close the batch
    thread_pool.close_batch(batch_id);
#else  // WITH_THREADS
    for (auto t_it = std::begin(bundle.templates());
         t_it != std::end(bundle.templates()); ++t_it) {
      refine_coeff_itvl(this, std::ref(*t_it));
    }
#endif // WITH_THREADS
    delete itvl_finder;
  } catch (SymbolicAlgebra::symbol_evaluation_error &e) {

    delete itvl_finder;
    std::ostringstream oss;

    oss << "The symbol \"" << e.get_symbol_name() << "\" is unknown";
    throw std::domain_error(oss.str());
  }

  LinearAlgebra::Vector<double> lower_bounds, upper_bounds;
  for (auto it = std::begin(max_coeffs); it != std::end(max_coeffs); ++it) {
    upper_bounds.push_back(*it);
  }
  for (auto it = std::begin(min_coeffs); it != std::end(min_coeffs); ++it) {
    lower_bounds.push_back(*it);
  }

  Bundle res(std::move(new_dirs), std::move(lower_bounds),
             std::move(upper_bounds), bundle.templates());

  for (const auto &len: bundle.edge_lengths()) {
    if (len > _edge_threshold) {
      std::ostringstream ss;

      ss << "CachedEvolver<double>::operator(): one of the "
            "computed bundle edge lengths is larger than "
            "the set threshold, i.e., "
         << _edge_threshold;

      throw std::runtime_error(ss.str());
    }
  }

  if (this->mode == ONE_FOR_ONE) {
    res.canonize();
  }

  return res;
}

SetsUnion<Polytope> synthesize(const Evolver<double> &evolver,
                               const Bundle &bundle,
                               const SetsUnion<Polytope> &parameter_set,
                               const std::shared_ptr<STL::Atom> atom)
{
  using namespace std;
  using namespace SymbolicAlgebra;

  auto &ds = evolver.dynamical_system();

  if (bundle.dim() != ds.variables().size()) {
    throw std::domain_error("The bundle and the dynamical system must "
                            "have the same number of dimensions.");
  }

  if (parameter_set.dim() != ds.parameters().size()) {
    throw std::domain_error("The parameter set and the parameter vector "
                            "must have the same number of dimensions.");
  }

  if (ds.parameters().size() == 0) {
    throw std::domain_error("Nothing to synthesize: the parameter set "
                            "is empty.");
  }

  SetsUnion<Polytope> result = parameter_set;

  std::vector<Symbol<>> alpha = get_symbol_vector<double>("f", bundle.dim());

  for (auto t_it = std::begin(bundle.templates());
       t_it != std::end(bundle.templates());
       ++t_it) { // for each parallelotope

    Parallelotope P = bundle.get_parallelotope(*t_it);
    std::vector<Expression<>> genFun
        = build_instantiated_generator_functs(alpha, P);

    const auto fog = replace_in(ds.dynamics(), ds.variables(), genFun);

    // compose sigma(f(gamma(x)))
    Expression<>::replacement_type repl;
    for (unsigned int j = 0; j < ds.dim(); j++) {
      repl[ds.variable(j)] = fog[j];
    }

    Expression<> sofog = atom->get_expression();
    sofog.replace(repl);

    // compute the Bernstein control points
    auto controlPts = get_Bernstein_coefficients(alpha, sofog);

    Polytope constraints(ds.parameters(), controlPts);
    result = ::intersect(result, constraints);
  }

  return result;
}