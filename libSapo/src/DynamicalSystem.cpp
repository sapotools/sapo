/**
 * @file DynamicalSystem.cpp
 * @author Alberto Casagrande <acasagrande@units.it>
 * @brief  Represent dynamic system laws
 * @version 0.1
 * @date 2022-05-21
 *
 * @copyright Copyright (c) 2016-2022
 */

#include "DynamicalSystem.h"

#include <utility>

#ifdef WITH_THREADS
#include <mutex>
#include <shared_mutex>

#include "SapoThreads.h"
#endif // WITH_THREADS

#include "BaseConverter.h"
#include "VarsGenerator.h"

/**
 * @brief Avoid \f$-0\f$
 *
 * This macro avoids \f$-0\f$ by replacing it by \f$0\f$.
 */
#define AVOID_NEG_ZERO(value) ((value) == 0 ? 0 : (value))

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
 * @param alpha
 * @param f
 * @param direction
 * @return std::vector<SymbolicAlgebra::Expression<>>
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

  return BaseConverter(alpha, Lfog).getBernCoeffsMatrix();
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
 * \f[q + ((\textrm{alpha} \circ \beta)^T \cdot D)^T\f]
 * where \f$\circ\f$ is the Hadamard product, and
 * \f$q\f$, \f$\beta\f$, and $D$ are the base vector,
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

template<>
Bundle
DynamicalSystem<double>::transform(const Bundle &bundle,
                                   const Polytope &parameter_set,
                                   const transformation_mode t_mode) const
{

  using namespace std;
  using namespace SymbolicAlgebra;
  using namespace LinearAlgebra;

  if (bundle.dim() != _variables.size()) {
    throw std::domain_error("The bundle and the dynamic laws must have the "
                            "same number of dimensions.");
  }

  if (_parameters.size() != parameter_set.dim()) {
    throw std::domain_error("The vector of parameters and the parameter "
                            "set must have same number of dimensions.");
  }

  vector<CondSyncUpdater<double, std::less<double>>> max_coeffs(bundle.size());
  vector<CondSyncUpdater<double, std::greater<double>>> min_coeffs(
      bundle.size());
  MinMaxCoeffFinder *itvl_finder;

  if (_parameters.size() == 0) {
    itvl_finder = new MinMaxCoeffFinder();
  } else {
    itvl_finder = new ParamMinMaxCoeffFinder(_parameters, parameter_set);
  }

  std::vector<Symbol<double>> alpha = get_symbol_vector("f", dim());

  auto refine_coeff_itvl = [&max_coeffs, &min_coeffs, &bundle, &alpha,
                            &itvl_finder,
                            &t_mode](const DynamicalSystem<double> *ds,
                                     const unsigned int template_num) 
  {
    Parallelotope P = bundle.get_parallelotope(template_num);

    const auto &genFun = build_instantiated_generator_functs(alpha, P);
    const auto genFun_f = replace_in(ds->dynamics(), ds->variables(), genFun);

    const auto &template_i = bundle.get_template(template_num);

    unsigned int dir_b;

    // for each direction
    const size_t num_of_dirs
        = (t_mode == ONE_FOR_ONE ? template_i.size() : bundle.size());

    for (unsigned int j = 0; j < num_of_dirs; j++) {
      if (t_mode == ONE_FOR_ONE) {
        dir_b = template_i[j];
      } else {
        dir_b = j;
      }
      auto coefficients
          = compute_Bern_coefficients(alpha, genFun_f, bundle.get_direction(dir_b));

      auto coeff_itvl = (*itvl_finder)(coefficients);

      min_coeffs[dir_b].update(coeff_itvl.first);
      max_coeffs[dir_b].update(coeff_itvl.second);
    }
  };

  try {
#ifdef WITH_THREADS
    ThreadPool::BatchId batch_id = thread_pool.create_batch();

    for (unsigned int i = 0; i < bundle.num_of_templates(); i++) {
      // submit the task to the thread pool
      thread_pool.submit_to_batch(batch_id, refine_coeff_itvl, this, i);
    }

    // join to the pool threads
    thread_pool.join_threads(batch_id);

    // close the batch
    thread_pool.close_batch(batch_id);
#else  // WITH_THREADS
    for (unsigned int i = 0; i < bundle.num_of_templates(); i++) {
      refine_coeff_itvl(this, i);
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

  Bundle res(bundle.directions(), lower_bounds, upper_bounds,
             bundle.templates());

  if (t_mode == ONE_FOR_ONE) {
    res.canonize();
  }

  return res;
}



PolytopesUnion synthesize(const DynamicalSystem<double> &ds,
                          const Bundle &bundle,
                          const PolytopesUnion &parameter_set,
                          const std::shared_ptr<STL::Atom> atom)
{
  using namespace std;
  using namespace SymbolicAlgebra;

  if (bundle.dim() != ds.variables().size()) {
    throw std::domain_error("The bundle and the dynamical system must "
                            "have the same number of dimensions.");
  }

  if (parameter_set.dim() != ds.parameters().size()) {
    throw std::domain_error("The parameter set and the parameter vector "
                            "must have the same number of dimensions.");
  }

  PolytopesUnion result = parameter_set;

  std::vector<Symbol<>> alpha = get_symbol_vector("f", bundle.dim());

  for (unsigned int i = 0; i < bundle.num_of_templates();
       i++) { // for each parallelotope

    Parallelotope P = bundle.get_parallelotope(i);
    std::vector<Expression<>> genFun
        = build_instantiated_generator_functs(alpha, P);

    const std::vector<Expression<>> fog
        = replace_in(ds.dynamics(), ds.variables(), genFun);

    // compose sigma(f(gamma(x)))
    Expression<>::replacement_type repl;
    for (unsigned int j = 0; j < ds.dim(); j++) {
      repl[ds.variable(j)] = fog[j];
    }

    Expression<> sofog = atom->get_expression();
    sofog.replace(repl);

    // compute the Bernstein control points
    std::vector<Expression<>> controlPts
        = BaseConverter(alpha, sofog).getBernCoeffsMatrix();

    Polytope constraints(ds.parameters(), controlPts);
    result = ::intersect(result, constraints);
  }

  return result;
}