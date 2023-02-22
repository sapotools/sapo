/**
 * @file Evolver.cpp
 * @author Alberto Casagrande <acasagrande@units.it>
 * @brief Transform geometric sets subject to a discrete system
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
#include "ErrorHandling.h"

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
 * @tparam T is the numeric type of the coefficients
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

template<typename T>
std::vector<SymbolicAlgebra::Expression<T>>
replace_in(const std::vector<SymbolicAlgebra::Expression<T>> &expressions,
           const std::vector<SymbolicAlgebra::Symbol<T>> &vars,
           const std::vector<SymbolicAlgebra::Expression<T>> &subs)
{
  using namespace SymbolicAlgebra;

  typename Expression<T>::replacement_type repl;

  for (unsigned int k = 0; k < vars.size(); ++k) {
    repl[vars[k]] = subs[k];
  }

  std::vector<Expression<T>> results;
  for (auto ex_it = std::begin(expressions); ex_it != std::end(expressions);
       ++ex_it) {
    results.push_back(Expression<T>(*ex_it).replace(repl));
  }

  return results;
}

/**
 * @brief Compute the Bernstein coefficients for a direction
 *
 * @tparam T is the numeric type of the coefficients
 * @param alpha is the vector of the variables
 * @param f is the function whose Bernstein coefficients must be computed
 * @param direction is the direction along
 * @return the Bernstein coefficients for `f` on `direction`
 */
template<typename T>
std::vector<SymbolicAlgebra::Expression<T>> compute_Bernstein_coefficients(
    const std::vector<SymbolicAlgebra::Symbol<T>> &alpha,
    const std::vector<SymbolicAlgebra::Expression<T>> &f,
    const LinearAlgebra::Vector<T> &direction)
{
  SymbolicAlgebra::Expression<T> Lfog = 0;
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
 * @tparam T is the numeric type of the coefficients
 * @param P is a parallelotope
 * @param q are the variables associated to the parallelotope's base vertex
 * @param beta are the variables associated to the parallelotope's lengths
 * @return the symbolic equations representing the the variable
 *         substitutions for `P`
 */
template<typename T>
typename SymbolicAlgebra::Expression<T>::replacement_type
get_subs_from(const Parallelotope &P,
              const std::vector<SymbolicAlgebra::Symbol<T>> &q,
              const std::vector<SymbolicAlgebra::Symbol<T>> &beta)
{
  using namespace SymbolicAlgebra;

  const LinearAlgebra::Vector<T> &base_vertex = P.base_vertex();
  const LinearAlgebra::Vector<T> &lengths = P.lengths();

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
 * @tparam T is a numeric type
 * @param alpha is a vector of variables
 * @param P is the considered parallelotope
 * @return The generator function of `P`
 */
template<typename T>
std::vector<SymbolicAlgebra::Expression<T>>
build_generator_functions(const std::vector<SymbolicAlgebra::Symbol<T>> &alpha,
                          const Parallelotope &P)
{
  using namespace LinearAlgebra;

  std::vector<SymbolicAlgebra::Expression<T>> gen_functs;
  std::copy(P.base_vertex().begin(), P.base_vertex().end(),
            std::back_inserter(gen_functs));

  const std::vector<Vector<T>> &generators = P.generators();

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
 * @tparam T is a numeric type
 * @param base is the base vertex variable vector
 * @param alpha is a vector of variables
 * @param lambda is the generator length variable vector
 * @return The generator function of a parallelotope having `base`
 *       as base vertex,`lambda` as generator length vector, and
 *       `D` as generator matrix
 */
template<typename T>
std::vector<SymbolicAlgebra::Expression<T>> build_symbolic_generator_functions(
    const std::vector<SymbolicAlgebra::Symbol<T>> &base,
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
 * @tparam T is the numeric type of the value
 * @tparam COND is the type of the condition to update
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
template<typename T>
class MinMaxCoeffFinder
{
  /**
   * @brief Evaluate a Bernstein coefficient
   *
   * @param coefficient is the symbolical representation of Bernstein
   *    coefficient
   * @return The numerical evaluation of Bernstein coefficient
   */
  T eval_coeff(const SymbolicAlgebra::Expression<T> &coefficient) const;

public:
  /**
   * @brief Create MaxCoeffType objects.
   */
  MinMaxCoeffFinder() {}

  /**
   * @brief Find the interval containing the Bernstein coefficients.
   *
   * @param coefficients is the vector of symbolical Bernstein coefficients
   * @return The pair minimum-maximum among all the Bernstein
   *          coefficients in `coefficients`
   */
  virtual std::pair<T, T> operator()(
      const std::vector<SymbolicAlgebra::Expression<T>> &coefficients) const;

  virtual ~MinMaxCoeffFinder() {}
};

/**
 * @brief  A class to find the minimum and the maximum parametric Bernstein
 * coefficients
 */
template<typename T>
class ParamMinMaxCoeffFinder : public MinMaxCoeffFinder<T>
{
  const std::vector<SymbolicAlgebra::Symbol<T>> &params;
  const Polytope &paraSet;

  /**
   * @brief Evaluate the parametric Bernstein coefficient upper-bound
   *
   * @param coefficient is the symbolical representation of Bernstein
   *                  coefficient
   * @return The numerical evaluation of parametric Bernstein
   *         coefficient upper-bound
   */
  T maximize_coeff(const SymbolicAlgebra::Expression<T> &coefficient) const;

  /**
   * @brief Evaluate the parametric Bernstein coefficient lower-bound
   *
   * @param coefficient is the symbolical representation of Bernstein
   *                  coefficient
   * @return The numerical evaluation of parametric Bernstein coefficient
   *         lower-bound
   */
  T minimize_coeff(const SymbolicAlgebra::Expression<T> &coefficient) const;

public:
  /**
   * @brief Constructor
   *
   * @param params is the list of parameters
   * @param paraSet is the set of admissible values for parameters
   */
  ParamMinMaxCoeffFinder(const std::vector<SymbolicAlgebra::Symbol<T>> &params,
                         const Polytope &paraSet):
      MinMaxCoeffFinder<T>(),
      params(params), paraSet(paraSet)
  {
  }

  /**
   * @brief Find the interval containing the Bernstein coefficients.
   *
   * @param coefficients is the vector of symbolical Bernstein coefficients
   * @return The pair minimum-maximum among all the Bernstein
   *          coefficients in `coefficients`
   */
  std::pair<T, T> operator()(
      const std::vector<SymbolicAlgebra::Expression<T>> &coefficients) const;

  ~ParamMinMaxCoeffFinder() {}
};

template<typename T>
T MinMaxCoeffFinder<T>::eval_coeff(
    const SymbolicAlgebra::Expression<T> &coefficient) const
{
  return AVOID_NEG_ZERO(coefficient.evaluate());
}

template<typename T>
T ParamMinMaxCoeffFinder<T>::maximize_coeff(
    const SymbolicAlgebra::Expression<T> &coefficient) const
{
  return paraSet.maximize(params, coefficient).objective_value();
}

template<typename T>
T ParamMinMaxCoeffFinder<T>::minimize_coeff(
    const SymbolicAlgebra::Expression<T> &coefficient) const
{
  return paraSet.minimize(params, coefficient).objective_value();
}

template<typename T>
std::pair<T, T> MinMaxCoeffFinder<T>::operator()(
    const std::vector<SymbolicAlgebra::Expression<T>> &coefficients) const
{
  auto b_coeff_it = coefficients.begin();

  T max_value = eval_coeff(*b_coeff_it);
  T min_value = max_value;

  for (++b_coeff_it; b_coeff_it != coefficients.end(); ++b_coeff_it) {
    T coeff_value = eval_coeff(*b_coeff_it);

    if (coeff_value > max_value) {
      max_value = coeff_value;
    }

    if (coeff_value < min_value) {
      min_value = coeff_value;
    }
  }

  return std::pair<T, T>(std::move(min_value), std::move(max_value));
}

template<typename T>
std::pair<T, T> ParamMinMaxCoeffFinder<T>::operator()(
    const std::vector<SymbolicAlgebra::Expression<T>> &coefficients) const
{
  auto b_coeff_it = coefficients.begin();

  T max_value = maximize_coeff(*b_coeff_it);
  T min_value = minimize_coeff(*b_coeff_it);

  for (++b_coeff_it; b_coeff_it != coefficients.end(); ++b_coeff_it) {
    T coeff_value = maximize_coeff(*b_coeff_it);

    if (coeff_value > max_value) {
      max_value = coeff_value;
    }

    coeff_value = minimize_coeff(*b_coeff_it);
    if (coeff_value < min_value) {
      min_value = coeff_value;
    }
  }

  return std::pair<T, T>(std::move(min_value), std::move(max_value));
}

template<typename T>
inline typename SymbolicAlgebra::Expression<T>::interpretation_type
build_interpretation(const std::vector<SymbolicAlgebra::Symbol<T>> &symbols,
                     const std::vector<T> &values)
{
  if (symbols.size() != values.size()) {
    SAPO_ERROR("symbols and values sizes are not the same", std::domain_error);
  }

  typename SymbolicAlgebra::Expression<T>::interpretation_type inter;
  for (size_t i = 0; i < symbols.size(); ++i) {
    inter[symbols[i]] = values[i];
  }

  return inter;
}

template<typename T>
LinearAlgebra::Vector<T> get_approx_center(const Bundle &bundle)
{
  using namespace LinearAlgebra;
  std::vector<T> approx_c;

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

  return approx_c / T(bundle.num_of_templates());
}

template<typename T>
std::vector<SymbolicAlgebra::Expression<T>>
average_dynamics(const DynamicalSystem<T> &ds, const Polytope &parameter_set)
{
  using namespace SymbolicAlgebra;

  std::vector<T> approx_c = get_approx_center<T>(Bundle(parameter_set));

  auto inter{build_interpretation(ds.parameters(), approx_c)};

  std::vector<Expression<T>> avg_ds;
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

template<typename T>
inline LinearAlgebra::Vector<T> compute_image(
    const std::vector<SymbolicAlgebra::Expression<T>> &dynamical_system,
    const std::vector<SymbolicAlgebra::Symbol<T>> &variables,
    const LinearAlgebra::Vector<T> &point)
{
  auto point_interpretation = build_interpretation(variables, point);
  return evaluate(dynamical_system, point_interpretation);
}

template<typename T>
std::vector<LinearAlgebra::Vector<T>> compute_images(
    const std::vector<SymbolicAlgebra::Expression<T>> &dynamical_system,
    const std::vector<SymbolicAlgebra::Symbol<T>> &variables,
    const std::vector<LinearAlgebra::Vector<T>> &points)
{
  std::vector<LinearAlgebra::Vector<T>> images;

  for (const auto &point: points) {
    images.push_back(compute_image(dynamical_system, variables, point));
  }

  return images;
}

/**
 * @brief Compute a unit vector orthogonal to n-1 vectors
 *
 * This function computes a unit vector orthogonal to n-1 linearly independent
 * vectors in $\mathbb{R}^n$, $v_1$, $v_{2}$, ..., $v_{n-1}$ by evaluating the
 * determinant of the matrix $M = [e; v_1; v_2; \ldots; v_{n-1}]$ where
 * $e=[e_1, \ldots, e_{n-1}]$ and $e_1$, ..., $e_{n-1}$ is a canonical base for
 * $\mathbb{R}^n$. This determinant has the form $d_1*e_1 + \ldots +
 * d_{n}*e_{n}$ and can be written as the vector $d = [d_1, \ldots, d_n]$. By
 * construction, $d$ is orthogonal to any $v_i$ and the function returns
 * $d/||d||_2$. The determinant of $M$ can be computed as
 * $$det(M) = \sum_{j=1}^n e_j*(-1)^{1+j}*det(M_{1,j})$$
 * where $M_{1,j}$ is the minor of the entry in the 1st row and $j$th column.
 * Thus, $d_j = (-1)^{1+j}*det(M_{1,j})$. We know that
 * $det(M_{1,j})=det(M_{1,j}^T)$ and $M_{1,j}^T$ is the matrix
 * $M' = [v_1; v_2; \ldots; v_{n-1}]^T$ without the $j$th row; let us denote
 * it by ${M'}_{j}$.
 * This function iteratively builds ${M'}_{n-1}$, ${M'}_{n-2}$, ...,
 * ${M'}_{1}$ and computes their determinants to return $d$.
 *
 * @param vectors is a vector of n-1 linearly independent vectors in
 * $\mathbb{R}^n$
 * @return a unit vector orthogonal to all the vectors in vectors
 */
LinearAlgebra::Vector<double> compute_orthogonal_direction(
    const std::vector<LinearAlgebra::Vector<double>> &vectors)
{
  using namespace LinearAlgebra;
  using namespace LinearAlgebra::Dense;

  // orth_vector will eventually store the aimed orthogonal vector
  LinearAlgebra::Vector<double> orth_vector(vectors[0].size());

  // build ${M'}_{n-1}$ and save the missing row in missing
  auto Mp_j = transpose(vectors);
  auto missing = std::move(Mp_j[orth_vector.size() - 1]);
  Mp_j.pop_back();

  // for all the rows in $M'$
  for (size_t j = orth_vector.size() - 1; j > 0; --j) {
    orth_vector[j] = determinant(Mp_j) * (j % 2 == 0 ? 1 : -1);

    // swap the missing row and the (j-1)-th
    swap(Mp_j[j - 1], missing);
  }
  orth_vector[0] = determinant(Mp_j);

  return orth_vector / norm_2(orth_vector);
}

/**
 * @brief Extract an element from a vector
 *
 * This function extracts the i-th element from a
 * C++ vector in time \f$O(1)\f$. This is done in
 * two phases: 1) the i-th element and the last
 * one are swapped 2) the vector size is decreased
 * by one.
 *
 * After the execution of the function:
 * - the element that was the last element of the
 *   vector becomes the i-th element
 * - the element that was the i-th element is
 *   removed from the vector
 * - all the remaining elements maintain their
 *   respective positions
 *
 * The two functions `extract_value` and
 * `insert_value` must be such that
 * ```
 * insert_value(v,i,extract_value(v,i))==v
 * ```
 *
 * @tparam T is the element type
 * @param v is a vector
 * @param i is the position of the element to be
 *      removed
 * @return the value that was in `v[i]` before
 *      the execution of this function
 */
template<typename T>
T extract_value(std::vector<T> &v, const size_t i)
{
  T value;
  std::swap(v[v.size() - 1], value);
  std::swap(v[i], value);
  v.pop_back();

  return value;
}

/**
 * @brief Insert an element in a vector
 *
 * This function inserts a value in a C++ vector
 * in position i-th in time \f$O(1)\f$. This is
 * done in two phases: 1) the new value is inserted
 * at the end of the vector 2) the i-th value and
 * the last one in the vector are swapped.
 *
 * After the execution of the function:
 * - the element that was the i-th position becomes
 *   the last one in the vector
 * - the inserted value is the i-th element in
 *   in the vector
 * - all the remaining values maintain their
 *   respective positions
 *
 * The two functions `extract_value` and
 * `insert_value` must be such that
 * ```
 * insert_value(v,i,extract_value(v,i))==v
 * ```
 *
 * @tparam T is the element type
 * @param v is a vector
 * @param i is the position in which the new value
 *     must be placed
 * @param value is the value to be inserted in the
 *     i-th position
 */
template<typename T>
void insert_value(std::vector<T> &v, const size_t i, T &&value)
{
  v.push_back(std::move(value));
  std::swap(v[v.size() - 1], v[i]);
}

template<typename T>
inline void fix_new_dir_verse(const LinearAlgebra::Vector<T> &dir_image,
                              LinearAlgebra::Vector<T> &new_dir)
{
  using namespace LinearAlgebra;

  if (new_dir * dir_image < 0) {
    new_dir *= T(-1);
  }
}

/**
 * @brief Get the basis vertices of a parallelotope
 *
 * Every parallelotope is described by a vector of generator
 * vectors \f$G\f$, their lengths \f$\lambda\f$, and a
 * base vertices \f$b\f$. The parallelotope basis is
 * the set of vectors \f$\lambda_i*G_i\f$.
 * The parallelotope basis vertices are vertices of the
 * parallelotope in the form \f$\lambda_i*G_i+b\f$.
 *
 * @param p is a parallelotope
 * @return the vector of parallelotope basis vertices
 */
std::vector<LinearAlgebra::Vector<double>>
get_parallelotope_basis_vertices(const Parallelotope &p)
{
  using namespace LinearAlgebra;

  std::vector<Vector<double>> basis_vertices;
  for (size_t j = 0; j < p.dim(); ++j) {
    auto basis_vertex = p.generators()[j] * p.lengths()[j] + p.base_vertex();
    basis_vertices.push_back(std::move(basis_vertex));
  }

  return basis_vertices;
}

/**
 * @brief Get the parallelotope basis image
 *
 * Every parallelotope is described by a vector of generator
 * vectors \f$G\f$, the associated lengths \f$\lambda\f$, and a
 * base vertices \f$b\f$. The parallelotope basis is
 * the set of vectors \f$\lambda_i*G_i\f$.
 * The parallelotope basis vertices are vertices of the
 * parallelotope in the form \f$\lambda_i*G_i+b\f$.
 * The image of the parallelotope basis through a function
 * \f$f\f$ is the set of vectors \f$f(lambda_i*G_i+b)-f(b)\f$.
 *
 * @param p is a parallelotope
 * @return the vector of the parallelotope basis image
 */
std::vector<LinearAlgebra::Vector<double>> get_parallelotope_basis_images(
    const std::vector<SymbolicAlgebra::Expression<double>> &dynamical_system,
    const std::vector<SymbolicAlgebra::Symbol<double>> &variables,
    const Parallelotope &p)
{
  using namespace LinearAlgebra;

  auto edge_vertices = get_parallelotope_basis_vertices(p);

  const auto base_vertex_image
      = compute_image(dynamical_system, variables, p.base_vertex());

  std::vector<Vector<double>> basis_images;
  for (const auto &edge_vertex: edge_vertices) {
    auto vertex_image
        = compute_image(dynamical_system, variables, edge_vertex);
    vertex_image -= base_vertex_image;
    basis_images.push_back(std::move(vertex_image));
  }

  return basis_images;
}

/**
 * @brief Compute new directions for a bundle
 *
 * This method computes new directions for a bundle subject to a
 * parametric dynamical system. This is achieved by computing
 * the images of each parallelotope edge on the base-vertex side
 * so to map the original parallelotope in a new parallelotope
 * whose directions are the new dictions.
 *
 * @param bundle is a bundle
 * @param ds is a dynamical system
 * @param parameter_set is a parameter set
 * @return a vector of new directions for the bundle
 */
std::vector<LinearAlgebra::Vector<double>>
compute_new_directions(const Bundle &bundle, const DynamicalSystem<double> &ds,
                       const Polytope &parameter_set)
{
  using namespace LinearAlgebra;

  // if no direction is dynamic return the old direction vector
  if (bundle.adaptive_directions().size() == 0) {
    return bundle.directions();
  }

  std::vector<Vector<double>> new_dirs(bundle.directions());

  auto avg_ds = average_dynamics(ds, parameter_set);

  for (const auto &bundle_template: bundle.templates()) {
    if (bundle_template.is_adaptive()) {
      Parallelotope p = bundle.get_parallelotope(bundle_template);

      auto basis_images
          = get_parallelotope_basis_images(avg_ds, ds.variables(), p);

      // if an basis has length 0, then the corresponding generator
      // does not change
      for (size_t i = 0; i < p.dim(); ++i) {
        if (p.lengths()[i] == 0) {
          basis_images[i] = p.generators()[i];
        }
      }

      for (const auto &adaptive_index: bundle_template.adaptive_indices()) {
        const auto &direction_index = bundle_template[adaptive_index];
        if (bundle.is_direction_adaptive(direction_index)) {
          // remove the basis vector of the considered direction
          Vector<double> basis_direction
              = extract_value(basis_images, adaptive_index);

          // compute one of vectors orthogonal to sampling images
          new_dirs[direction_index]
              = compute_orthogonal_direction(basis_images);

          // insert again the basis vector of the considered direction
          insert_value(basis_images, adaptive_index,
                       std::move(basis_direction));

          fix_new_dir_verse(bundle.get_direction(direction_index),
                            new_dirs[direction_index]);
        }
      }
    }
  }

  return new_dirs;
}

/**
 * @brief Bundle bound refiner
 *
 * When a bundle is transformed according to a parametric
 * dynamical system, the bounds of its image on a set of
 * new directions must be computed. For each of the new
 * directions, we must evaluate the maximum and the minimum
 * among the Bernstein coefficients of a set of generator
 * functions which depend on the considered direction, on
 * the bundle parallelotopes, and on the dynamical system.
 * This class manages the computation of bundle image bounds
 * and provides an practical interface to process bundle
 * templates one by one.
 *
 * @tparam T is the constant numeric type
 */
template<typename T>
class BoundRefiner
{
  const Bundle _bundle;                        //!< the considered bundle
  const DynamicalSystem<T> &_dynamical_system; //!< the dynamical system

  std::vector<CondSyncUpdater<T, std::greater<T>>>
      _lower_bound; //!< the vector of identified lower bounds
  std::vector<CondSyncUpdater<T, std::less<T>>>
      _upper_bound; //!< the vector of identified upper bounds

  const std::vector<LinearAlgebra::Vector<T>>
      &_new_directions; //!< the vector of bundle new directions

  const std::vector<SymbolicAlgebra::Symbol<T>>
      _alpha; //!< the dynamical law variables

  const std::vector<SymbolicAlgebra::Symbol<T>>
      _lambda; //!< the variables representing parallelotope edge lengths
  const std::vector<SymbolicAlgebra::Symbol<T>>
      _base; //!< the variables representing parallelotope base

  MinMaxCoeffFinder<T>
      *_minmax_finder; //!< an object minimize and maximize Bernstein
                       //!< coefficient in the parameter set

  /**
   * @brief Parallelotope processor
   *
   * This class refines the bundle image bounds considering one
   * of the bundle templates.
   */
  class ParallelotopeProcessor
  {
    const Parallelotope
        _parallelotope; //!< the parallelotope of the considered template

    std::vector<SymbolicAlgebra::Expression<T>>
        _generator_functions; //!< the generator functions that maps [0,1]^n in
                              //!< the parallelotope image

    typename SymbolicAlgebra::Expression<T>::interpretation_type
        _parallelotope_interpretation; //!< the interpretation of the lambda
                                       //!< and beta variables for the
                                       //!< considered parallelotope

    BernsteinCache<T>
        *_cache; //!< a pointer to the symbolic Bernstein coefficient cache

    /**
     * @brief Get the symbolic Bernstein coefficients of a direction
     *
     * This method search for the symbolic Bernstein coefficients of a
     * direction in the cache. If it does not contain them, they are computed
     * and stored in the cache. Finally, a reference to the stored symbolic
     * Bernstein coefficients is returned.
     *
     * @param alpha is the vector of alpha variables to be used
     * @param direction is the direction whose symbolic Bernstein coefficients
     * are aimed
     * @return a reference to the computed symbolic Bernstein coefficients of a
     * direction
     */
    const std::vector<SymbolicAlgebra::Expression<T>> &
    get_symbolic_coefficients(
        const std::vector<SymbolicAlgebra::Symbol<T>> alpha,
        const LinearAlgebra::Vector<T> &direction)
    {
      if (_cache->coefficients_in_cache(_parallelotope,
                                        direction)) { // Bernstein coefficients
                                                      // have been computed

        return _cache->get_coefficients(_parallelotope, direction);
      }

      auto coefficients = compute_Bernstein_coefficients(
          alpha, _generator_functions, direction);
      return _cache->save_coefficients(_parallelotope, direction,
                                       std::move(coefficients));
    }

  public:
    /**
     * @brief A constructor
     *
     * @param refiner is the bound refiner that is about to process a bundle
     * parallelotope
     * @param bundle_template is one of the bundle templates
     * @param cache is a pointer to the symbolic Bernstein coefficient cache
     */
    ParallelotopeProcessor(BoundRefiner &refiner,
                           const BundleTemplate &bundle_template,
                           BernsteinCache<T> *cache):
        _parallelotope(refiner._bundle.get_parallelotope(bundle_template)),
        _cache(bundle_template.is_adaptive() ? nullptr : cache)
    {
      std::vector<SymbolicAlgebra::Expression<T>> genFun;
      if (_cache == nullptr) {
        genFun = build_generator_functions(refiner._alpha, _parallelotope);
      } else {
        genFun = build_symbolic_generator_functions(
            refiner._base, refiner._alpha, refiner._lambda,
            _parallelotope.generators());
      }

      _generator_functions
          = replace_in(refiner._dynamical_system.dynamics(),
                       refiner._dynamical_system.variables(), genFun);

      for (size_t i = 0; i < _parallelotope.dim(); ++i) {
        _parallelotope_interpretation[refiner._base[i]]
            = _parallelotope.base_vertex()[i];
        _parallelotope_interpretation[refiner._lambda[i]]
            = _parallelotope.lengths()[i];
      }
    }

    /**
     * @brief Get the direction bounds in the bundle image
     *
     * @param alpha is the alpha variable vector
     * @param direction is the direction whose bounds are aimed
     * @param minmax_finder is the object that search among Bernstein
     * coefficients for their maximum and minimum
     * @return a pair minimum-maximum bounds for the specified direction in the
     * image of the bundle through the dynamical system
     */
    std::pair<T, T>
    get_direction_bounds(const std::vector<SymbolicAlgebra::Symbol<T>> alpha,
                         const LinearAlgebra::Vector<T> &direction,
                         MinMaxCoeffFinder<T> *minmax_finder)
    {
      std::vector<SymbolicAlgebra::Expression<double>> coefficients;

      if (_cache == nullptr) {
        coefficients = compute_Bernstein_coefficients(
            alpha, _generator_functions, direction);
      } else {
        auto &symbolic_coefficients
            = get_symbolic_coefficients(alpha, direction);

        coefficients.reserve(symbolic_coefficients.size());
        for (auto &symbolic_coefficient: symbolic_coefficients) {
          coefficients.push_back(
              symbolic_coefficient.apply(_parallelotope_interpretation));
        }
      }

      return (*minmax_finder)(coefficients);
    }
  };

public:
  BoundRefiner(const Bundle &bundle,
               const DynamicalSystem<T> &dynamical_system,
               const Polytope &parameter_set,
               const std::vector<LinearAlgebra::Vector<T>> &new_directions):
      _bundle(bundle),
      _dynamical_system(dynamical_system), _lower_bound(bundle.size()),
      _upper_bound(bundle.size()), _new_directions(new_directions),
      _alpha(get_symbol_vector<T>("alpha", dynamical_system.dim())),
      _lambda(get_symbol_vector<T>("lambda", dynamical_system.dim())),
      _base(get_symbol_vector<T>("base", dynamical_system.dim()))
  {
    if (dynamical_system.parameters().size() == 0) {
      _minmax_finder = new MinMaxCoeffFinder<T>();
    } else {
      _minmax_finder = new ParamMinMaxCoeffFinder<T>(
          dynamical_system.parameters(), parameter_set);
    }
  }

  /**
   * @brief Process a bundle template
   *
   * @param bundle_template is a template of the considered bundle
   * @param mode is the bound computation mode, i.e., one-for-one or
   * all-for-one
   * @param cache is the symbolic Bernstein coefficient cache
   * @return a reference to the current object
   */
  BoundRefiner<T> &
  process_template(const BundleTemplate &bundle_template,
                   const typename Evolver<T>::evolver_mode &mode,
                   BernsteinCache<T> *cache = nullptr)
  {
    ParallelotopeProcessor processor(*this, bundle_template, cache);

    // for each direction
    const size_t num_of_directions
        = (mode == Evolver<T>::ONE_FOR_ONE ? bundle_template.dim()
                                           : _bundle.size());

    for (size_t j = 0; j < num_of_directions; j++) {
      auto direction_index
          = (mode == Evolver<T>::ONE_FOR_ONE ? bundle_template[j] : j);

      const auto &direction = _new_directions[direction_index];

      auto coefficients
          = processor.get_direction_bounds(_alpha, direction, _minmax_finder);

      _lower_bound[direction_index].update(coefficients.first);
      _upper_bound[direction_index].update(coefficients.second);
    }

    return *this;
  }

  /**
   * @brief Get the bundle image upper bounds
   *
   * @return the bundle image upper bounds as identified in the processed
   * templates
   */
  LinearAlgebra::Vector<T> get_upper_bounds() const
  {
    LinearAlgebra::Vector<T> upper_bounds;

    std::copy(_upper_bound.begin(), _upper_bound.end(),
              std::back_inserter(upper_bounds));

    return upper_bounds;
  }

  /**
   * @brief Get the bundle image lower bounds
   *
   * @return the bundle image lower bounds as identified in the processed
   * templates
   */
  LinearAlgebra::Vector<T> get_lower_bounds() const
  {
    LinearAlgebra::Vector<T> lower_bounds;

    std::copy(_lower_bound.begin(), _lower_bound.end(),
              std::back_inserter(lower_bounds));

    return lower_bounds;
  }

  /**
   * @brief The destroyer
   */
  ~BoundRefiner()
  {
    delete _minmax_finder;
  }

  friend class ParallelotopeProcessor;
};

template<>
Bundle Evolver<double>::operator()(const Bundle &bundle,
                                   const Polytope &parameter_set)
{
  using namespace std;
  using namespace SymbolicAlgebra;
  using namespace LinearAlgebra;

  if (bundle.dim() != _ds.variables().size()) {
    SAPO_ERROR("the bundle and the dynamic laws must have the "
               "same number of dimensions",
               std::domain_error);
  }

  if (parameter_set.dim() != _ds.parameters().size()) {
    SAPO_ERROR("the vector of parameters and the parameter "
               "set must have same number of dimensions",
               std::domain_error);
  }

  auto new_directions = compute_new_directions(bundle, _ds, parameter_set);

  BoundRefiner bound_refiner(bundle, _ds, parameter_set, new_directions);

  auto refine_bounds
      = [&bound_refiner](Evolver<double> *evolver,
                         const BundleTemplate &bundle_template) {
          bound_refiner.process_template(bundle_template, evolver->mode,
                                         evolver->_cache);
        };

  try {
#ifdef WITH_THREADS
    ThreadPool::BatchId batch_id = thread_pool.create_batch();

    for (auto t_it = std::begin(bundle.templates());
         t_it != std::end(bundle.templates()); ++t_it) {

      // submit the task to the thread pool
      thread_pool.submit_to_batch(batch_id, refine_bounds, this,
                                  std::ref(*t_it));
    }

    // join to the pool threads
    thread_pool.join_threads(batch_id);

    // close the batch
    thread_pool.close_batch(batch_id);
#else  // WITH_THREADS
    for (auto t_it = std::begin(bundle.templates());
         t_it != std::end(bundle.templates()); ++t_it) {
      refine_bounds(this, std::ref(*t_it));
    }
#endif // WITH_THREADS
  } catch (SymbolicAlgebra::symbol_evaluation_error &e) {
    std::ostringstream oss;

    oss << "the symbol \"" << e.get_symbol_name() << "\" is unknown";
    SAPO_ERROR(oss.str(), std::domain_error);
  }

  Bundle new_bundle(bundle.adaptive_directions(), std::move(new_directions),
                    bound_refiner.get_lower_bounds(),
                    bound_refiner.get_upper_bounds(), bundle.templates());

  for (const auto &len: bundle.edge_lengths()) {
    if (len > EDGE_MAX_LENGTH) {
      SAPO_ERROR("one of the computed bundle edge lengths is "
                 "larger than the set threshold (i.e., "
                     << EDGE_MAX_LENGTH << ")",
                 std::runtime_error);
    }
  }

  if (this->mode == ONE_FOR_ONE) {
    new_bundle.canonize();
  }

  return new_bundle;
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
    SAPO_ERROR("the bundle and the dynamical system must "
               "have the same number of dimensions",
               std::domain_error);
  }

  if (parameter_set.dim() != ds.parameters().size()) {
    SAPO_ERROR("the parameter set and the parameter vector "
               "must have the same number of dimensions",
               std::domain_error);
  }

  if (ds.parameters().size() == 0) {
    SAPO_ERROR("the parameter set is empty and nothing "
               "can be synthesized",
               std::domain_error);
  }

  SetsUnion<Polytope> result = parameter_set;

  std::vector<Symbol<>> alpha = get_symbol_vector<double>("f", bundle.dim());

  for (auto t_it = std::begin(bundle.templates());
       t_it != std::end(bundle.templates());
       ++t_it) { // for each parallelotope

    Parallelotope P = bundle.get_parallelotope(*t_it);
    std::vector<Expression<>> genFun = build_generator_functions(alpha, P);

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