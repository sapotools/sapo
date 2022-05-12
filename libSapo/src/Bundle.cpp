/**
 * @file Bundle.cpp
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @author Alberto Casagrande <acasagrande@units.it>
 * @brief Represent and manipulate bundles of parallelotopes 
 *        whose intersection represents a polytope
 * @version 0.2
 * @date 2022-05-04
 * 
 * @copyright Copyright (c) 2016-2022
 */

#include "Bundle.h"

#ifdef WITH_THREADS
#include <mutex>
#include <shared_mutex>

#include "SapoThreads.h"
#endif // WITH_THREADS

#include <limits>
#include <string>
#include <algorithm>
#include <functional>

#include "LinearAlgebra.h"

#define _USE_MATH_DEFINES //!< This macro enables the use of cmath constants

#include <cmath>

/**
 * @brief Avoid \f$-0\f$
 * 
 * This macro avoids \f$-0\f$ by replacing it by \f$0\f$.
 */
#define AVOID_NEG_ZERO(value) ((value) == 0 ? 0 : (value))

Bundle::Bundle(const Bundle &orig):
    directions(orig.directions), lower_bounds(orig.lower_bounds),
    upper_bounds(orig.upper_bounds), templates(orig.templates)
{
}

Bundle::Bundle(Bundle &&orig)
{
  swap(*this, orig);
}

void swap(Bundle &A, Bundle &B)
{
  std::swap(A.directions, B.directions);
  std::swap(A.upper_bounds, B.upper_bounds);
  std::swap(A.lower_bounds, B.lower_bounds);
  std::swap(A.templates, B.templates);
}

/**
 * Orthogonal proximity of v1 and v2, i.e.,
 * how close is the angle between v1 and v2 is to pi/2
 *
 * @param[in] v1 vector
 * @param[in] v2 vector
 * @returns orthogonal proximity
 */
double orthProx(LinearAlgebra::Vector<double> v1, LinearAlgebra::Vector<double> v2)
{
  return std::abs(LinearAlgebra::angle(v1, v2) - M_PI_2);
}

Bundle::Bundle(const std::vector<LinearAlgebra::Vector<double>> &directions,
               const LinearAlgebra::Vector<double> &lower_bounds,
               const LinearAlgebra::Vector<double> &upper_bounds,
               const std::vector<LinearAlgebra::Vector<int>> &templates):
    directions(directions),
    lower_bounds(lower_bounds), upper_bounds(upper_bounds),
    templates(templates)
{
}

Bundle::Bundle(std::vector<LinearAlgebra::Vector<double>> &&directions,
               LinearAlgebra::Vector<double> &&lower_bounds, LinearAlgebra::Vector<double> &&upper_bounds,
               std::vector<LinearAlgebra::Vector<int>> &&templates):
    directions(std::move(directions)),
    lower_bounds(std::move(lower_bounds)),
    upper_bounds(std::move(upper_bounds)), templates(std::move(templates))
{
  using namespace std;

  if (directions.size() == 0) {
    throw std::domain_error("Bundle::Bundle: directions must be non empty");
  }
  if (directions.size() != upper_bounds.size()) {
    throw std::domain_error("Bundle::Bundle: directions and upper_bounds "
                            "must have the same size");
  }
  if (directions.size() != lower_bounds.size()) {
    throw std::domain_error("Bundle::Bundle: directions and lower_bounds "
                            "must have the same size");
  }
  if (templates.size() > 0) {
    for (unsigned int i = 0; i < templates.size(); i++) {
      if (templates[i].size() != this->dim()) {
        throw std::domain_error("Bundle::Bundle: templates must have "
                                "as many columns as directions");
      }
    }
  } else {
    throw std::domain_error("Bundle::Bundle: templates must be non empty");
  }
}

Bundle &Bundle::operator=(Bundle &&orig)
{
  swap(*this, orig);

  return *this;
}

/**
 * Generate the polytope represented by the bundle
 *
 * @returns polytope represented by the bundle
 */
Bundle::operator Polytope() const
{
  using namespace std;
  using namespace LinearAlgebra;

  std::vector<Vector<double>> A;
  Vector<double> b;
  for (unsigned int i = 0; i < this->size(); i++) {
    A.push_back(this->directions[i]);
    b.push_back(this->upper_bounds[i]);
  }
  for (unsigned int i = 0; i < this->size(); i++) {
    A.push_back(-this->directions[i]);
    b.push_back(AVOID_NEG_ZERO(-this->lower_bounds[i]));
  }

  return Polytope(std::move(A), std::move(b));
}

Parallelotope Bundle::get_parallelotope(unsigned int i) const
{
  using namespace std;

  if (i > this->templates.size()) {
    throw std::domain_error("Bundle::get_parallelotope: i must be a valid row "
                            "for the template matrix");
  }

  vector<double> lbound, ubound;
  vector<LinearAlgebra::Vector<double>> Lambda;

  vector<int>::const_iterator it = std::begin(this->templates[i]);
  // upper facets
  for (unsigned int j = 0; j < this->dim(); j++) {
    const int idx = *it;
    Lambda.push_back(this->directions[idx]);
    ubound.push_back(this->upper_bounds[idx]);
    lbound.push_back(this->lower_bounds[idx]);

    ++it;
  }

  // TODO: Since Lambdas are always the same, check whether
  //       storing the LUP factorizations of Lambdas may give
  //       some speed-up.
  return Parallelotope(Lambda, lbound, ubound);
}

Bundle Bundle::get_canonical() const
{
  Bundle res(*this);

  res.canonize();

  return res;
}

Bundle& Bundle::canonize()
{
  // get current polytope
  Polytope bund = *this;
  for (unsigned int i = 0; i < this->size(); ++i) {
    lower_bounds[i] = bund.maximize(this->directions[i]).optimum();
    upper_bounds[i] = bund.minimize(this->directions[i]).optimum();
  }
  return *this;
}

/**
 * Check whether a vector is permutation of a sorted vector
 *
 * @param[in] v1 first vector
 * @param[in] sorted a sorted vector
 * @returns true is v1 is a permutation of v2
 */
bool isPermutationOfSorted(std::vector<int> v1,
                           const std::vector<int> &sorted)
{
  if (v1.size() != sorted.size()) {
    return false;
  }

  std::sort(std::begin(v1), std::end(v1));

  auto v1_it = std::begin(v1);
  auto v2_it = std::begin(sorted);
  while (v1_it != std::end(v1)) {
    if (*v1_it != *v2_it) {
      return false;
    }

    ++v1_it;
    ++v2_it;
  }

  return true;
}

/**
 * Check if v1 is a permutation of v2
 *
 * @param[in] v1 first vector
 * @param[in] v2 second vector
 * @returns true is v1 is a permutation of v2
 */
bool isPermutation(const std::vector<int> &v1, std::vector<int> v2)
{
  if (v1.size() != v2.size()) {
    return false;
  }

  std::sort(std::begin(v2), std::end(v2));

  return isPermutationOfSorted(v1, v2);
}

/**
 * @brief Test whether a row is permutation of another row
 * 
 * This method test whether the `i`-th row in the matrix `M`
 * is the permutation of another row in `M`.
 * 
 * @tparam T is the scalar type of `M`
 * @param M is the matrix whose rows must be tested
 * @param i is the index of the row to be tested
 * @return `true` if and only if the `i`-th row of `M`
 *         is the permutation of any other row in `M`
 */
template<typename T>
bool is_permutation_of_other_rows(const std::vector<std::vector<T>> &M,
                                  const unsigned int &i)
{
  using namespace std;

  vector<T> M_i = M[i];
  sort(std::begin(M_i), std::end(M_i));
  for (unsigned int j = 0; j < M.size(); j++) {
    if (j != i) {
      if (isPermutationOfSorted(M[j], M_i)) {
        return true;
      }
    }
  }

  return false;
}

/**
 * Maximum distance accumulation of a vector w.r.t. a set of vectors
 *
 * @param[in] vIdx index of the reference vector
 * @param[in] dirsIdx indexes of vectors to be considered
 * @param[in] dists pre-computed distances
 * @returns distance accumulation
 */
double maxOffsetDist(const int vIdx, const std::vector<int> &dirsIdx,
                     const LinearAlgebra::Vector<double> &dists)
{

  if (dirsIdx.empty()) {
    return 0;
  }

  double dist = dists[vIdx];
  for (unsigned int i = 0; i < dirsIdx.size(); i++) {
    dist = dist * dists[dirsIdx[i]];
  }
  return dist;
}

/**
 * Maximum distance accumulation of a set of vectors
 *
 * @param[in] dirsIdx indexes of vectors to be considered
 * @param[in] dists pre-computed distances
 * @returns distance accumulation
 */
double maxOffsetDist(const std::vector<int> &dirsIdx,
                     const LinearAlgebra::Vector<double> &dists)
{

  double dist = 1;
  for (unsigned int i = 0; i < dirsIdx.size(); i++) {
    dist = dist * dists[dirsIdx[i]];
  }
  return dist;
}

/**
 * Maximum distance accumulation of matrix
 *
 * @param[in] T matrix from which fetch the vectors
 * @param[in] dists pre-computed distances
 * @returns distance accumulation
 */
double maxOffsetDist(const std::vector<LinearAlgebra::Vector<int>> &T,
                     const LinearAlgebra::Vector<double> &dists)
{
  double maxdist = std::numeric_limits<double>::lowest();
  for (unsigned int i = 0; i < T.size(); i++) {
    maxdist = std::max(maxdist, maxOffsetDist(T[i], dists));
  }
  return maxdist;
}

/**
 * Maximum orthogonal proximity of a vector w.r.t. a set of vectors
 *
 * @param[in] directions is the direction matrix
 * @param[in] vIdx index of the reference vector
 * @param[in] dirsIdx indexes of vectors to be considered
 * @returns maximum orthogonal proximity
 */
double maxOrthProx(const std::vector<LinearAlgebra::Vector<double>> &directions,
                   const int vIdx, const std::vector<int> &dirsIdx)
{

  if (dirsIdx.empty()) {
    return 0;
  }

  double maxProx = 0;
  for (auto d_it = std::begin(dirsIdx); d_it != std::end(dirsIdx); ++d_it) {
    maxProx = std::max(maxProx, orthProx(directions[vIdx], directions[*d_it]));
  }
  return maxProx;
}

/**
 * Maximum orthogonal proximity within a set of vectors
 *
 * @param[in] directions is the direction matrix
 * @param[in] dirsIdx indexes of vectors to be considered
 * @returns maximum orthogonal proximity
 */
double maxOrthProx(const std::vector<LinearAlgebra::Vector<double>> &directions,
                   const std::vector<int> &dirsIdx)
{
  double maxProx = 0;
  for (unsigned int i = 0; i < dirsIdx.size(); i++) {
    for (unsigned int j = i + 1; j < dirsIdx.size(); j++) {
      maxProx = std::max(
          maxProx, orthProx(directions[dirsIdx[i]], directions[dirsIdx[j]]));
    }
  }
  return maxProx;
}

/**
 * Maximum orthogonal proximity of all the vectors of a matrix
 *
 * @param[in] directions is the direction matrix
 * @param[in] T collection of vectors
 * @returns maximum orthogonal proximity
 */
double maxOrthProx(const std::vector<LinearAlgebra::Vector<double>> &directions,
                   const std::vector<LinearAlgebra::Vector<int>> &T)
{
  double maxorth = std::numeric_limits<double>::lowest();
  for (auto T_it = std::begin(T); T_it != std::end(T); ++T_it) {
    maxorth = std::max(maxorth, maxOrthProx(directions, *T_it));
  }
  return maxorth;
}

/**
 * @brief A draft horse function to split a bundle
 * 
 * This recursive function is a private draft horse to 
 * split a bundle whose maximal magnitude, the maximal 
 * lenght of its generators, is greater than 
 * `max_magnitude`. The bundle is split into a
 * list of sub-bundles whose maximal magnitude is 
 * \f$\textrm{max_magnitude}*\textrm{split_ratio}\f$.
 * Each sub-bundle in the output list is built recursivelly 
 * by considering one bundle direction per time and:
 * 1. splitting the direction in the opportune number of 
 *    chunks, so to satisfy the sub-bundle maximal
 *    magnitude request;
 * 2. for each of the chuncks, setting the sub-bundle
 *    boundaries of the selected direction according 
 *    with the considered chunck and recursivelly 
 *    consider the remaining directions.
 * 
 * If \f$m\f$ in the maximal magnitude of the input bundle 
 * and the bundle itself has \f$d\f$ directions,
 * this method may produces upto
 * \f$(\frac{m}{\textrm{max_magnitude}})^d\f$ 
 * sub-bundles.
 * 
 * @param[out] res is the list of the resulting sub-bundles
 * @param[in, out] lower_bounds is the lower boundary vector  
 *                              of the in-building sub-bundles
 * @param[in, out] upper_bounds is the upper boundary vector 
 *                              of the in-building sub-bundles
 * @param[in] idx is the index of the direction to be split
 * @param[in] splitting is the to-be-split bundle
 * @param[in] max_magnitude is the maximal magnitude that
 *                          triggers the split
 * @param[in] split_ratio is the ratio between maximal
 *                        magnitude of the output bundles and
 *                        that that triggers the split
 * @return a reference to the list of sub-bundles produced
 *         splitting the input bundle
 */
std::list<Bundle> &
split_bundle(std::list<Bundle> &res, LinearAlgebra::Vector<double> &lower_bounds,
             LinearAlgebra::Vector<double> &upper_bounds, const unsigned int idx,
             const Bundle &splitting, const double &max_magnitude,
             const double &split_ratio)
{
  if (idx == splitting.get_directions().size()) {
    res.emplace_back(splitting.get_directions(), lower_bounds, upper_bounds,
                     splitting.get_templates());

    return res;
  }

  if (std::abs(splitting.get_upper_bound(idx) - splitting.get_lower_bound(idx))
      > max_magnitude) {
    double lower_bound = splitting.get_lower_bound(idx);

    do {
      const double upper_bound = std::min(
          lower_bound + split_ratio * max_magnitude,
          splitting.get_upper_bound(idx));

      lower_bounds[idx] = lower_bound;
      upper_bounds[idx] = upper_bound;
      split_bundle(res, lower_bounds, upper_bounds, idx + 1, splitting,
                   max_magnitude, split_ratio);

      lower_bound = upper_bound;
    } while (splitting.get_upper_bound(idx) != lower_bound);
  } else {
    lower_bounds[idx] = splitting.get_lower_bound(idx);
    upper_bounds[idx] = splitting.get_upper_bound(idx);
    split_bundle(res, lower_bounds, upper_bounds, idx + 1, splitting,
                 max_magnitude, split_ratio);
  }
  return res;
}

std::list<Bundle> Bundle::split(const double max_magnitude,
                                const double split_ratio) const
{
  std::list<Bundle> split_list;

  LinearAlgebra::Vector<double> upper_bounds(this->size());
  LinearAlgebra::Vector<double> lower_bounds(this->size());

  split_bundle(split_list, lower_bounds, upper_bounds, 0, *this,
               max_magnitude, split_ratio);

  return split_list;
}

Bundle Bundle::decompose(double dec_weight, int max_iters)
{
  using namespace std;
  using namespace LinearAlgebra;

  vector<double> offDists = this->edge_lengths();

  // get current template and try to improve it
  vector<Vector<int>> curT = this->templates;

  // get current template and try to improve it
  vector<Vector<int>> bestT = this->templates;
  int temp_card = this->templates.size();

  int i = 0;
  while (i < max_iters) {

    vector<Vector<int>> tmpT = curT;

    // generate random coordinates to swap
    unsigned int i1 = rand() % temp_card;
    int j1 = rand() % this->dim();

    // swap them
    tmpT[i1][j1] = rand() % this->size();

    if (!is_permutation_of_other_rows(tmpT, i1)) {
      std::vector<Vector<double>> A;
      for (unsigned int j = 0; j < this->dim(); j++) {
        A.push_back(this->directions[tmpT[i1][j]]);
      }

      Dense::LUP_Factorization<double> fact(A);
      try {
        fact.solve(Vector<double>(this->dim(), 0));

        double w1 = dec_weight * maxOffsetDist(tmpT, offDists)
                    + (1 - dec_weight) * maxOrthProx(this->directions, tmpT);
        double w2 = dec_weight * maxOffsetDist(bestT, offDists)
                    + (1 - dec_weight) * maxOrthProx(this->directions, bestT);

        if (w1 < w2) {
          bestT = tmpT;
        }
        curT = tmpT;
      } catch (...) {
        // The system Ax=b cannot be solved
      }
    }
    i++;
  }

  return Bundle(directions, lower_bounds, upper_bounds, bestT);
}

/**
 * @brief Replace some variables in expressions by other expressions
 * 
 * This function replaces all the occurences of the variable `vars[i]` 
 * in any expression in `expressions` with `subs[i]`.
 * 
 * @param expressions is the vector of expressions whose variable 
 *                    occurences must be replaced 
 * @param vars is the vector variables whose occurences must be 
 *             replaced 
 * @param subs is the vector of expressions that must replace 
 *             variable occurences
 * @return A vector of the containing the input expression 
 *         in which each occurence of the variables in `vars`
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
  for (auto ex_it = std::begin(expressions); ex_it != std::end(expressions); ++ex_it) {
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
compute_Bern_coeffs(const std::vector<SymbolicAlgebra::Symbol<>> &alpha,
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

double Bundle::MinMaxCoeffFinder::eval_coeff(
    const SymbolicAlgebra::Expression<> &bernCoeff) const
{
  double value = bernCoeff.evaluate();

  // TODO: The following conditional evaluation avoids -0
  //       values. Check the difference between -0 and 0.
  return AVOID_NEG_ZERO(value);
}

double Bundle::ParamMinMaxCoeffFinder::maximize_coeff(
    const SymbolicAlgebra::Expression<> &bernCoeff) const
{
  return paraSet.maximize(params, bernCoeff).optimum();
}

double Bundle::ParamMinMaxCoeffFinder::minimize_coeff(
    const SymbolicAlgebra::Expression<> &bernCoeff) const
{
  return paraSet.minimize(params, bernCoeff).optimum();
}

std::pair<double, double> Bundle::MinMaxCoeffFinder::find_coeffs_itvl(
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

std::pair<double, double> Bundle::ParamMinMaxCoeffFinder::find_coeffs_itvl(
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
std::vector<SymbolicAlgebra::Expression<>> build_instanciated_generator_functs(
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
    // correspond to 0-length dimensions in degenerous
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
 * accessed in a syncronized way according to a mutex.
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
  mutable std::shared_timed_mutex mutex; //!< The mutex that synchronizes the accesses
#endif
  T _value;  //!< the value stored in the object
  COND _cmp;      //!< the condition that must be satisfied to update the value

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
 * @brief Perform the parametric synthesis for an atom
 *
 * This method computes a set of parameters such that the
 * transformation from the current bundle satisfies the
 * provided atom.
 *
 * @param[in] variables is the vector of variables
 * @param[in] parameters is the vector of parameter
 * @param[in] dynamics is the vector of dynamic law expressions
 * @param[in] parameter_set is the initial parameter set
 * @param[in] atom is the specification to be satisfied
 * @return a subset of `pSet` such that the transformation
 *         from the current bundle through the dynamic laws
 *         when the parameters are in the returned set satisfies
 *         the atomic formula
 */
PolytopesUnion
Bundle::synthesize(const std::vector<SymbolicAlgebra::Symbol<>> &variables,
                   const std::vector<SymbolicAlgebra::Symbol<>> &parameters,
                   const std::vector<SymbolicAlgebra::Expression<>> &dynamics,
                   const PolytopesUnion &parameter_set,
                   const std::shared_ptr<STL::Atom> atom) const
{
  using namespace std;
  using namespace SymbolicAlgebra;

  PolytopesUnion result = parameter_set;

  std::vector<Symbol<>> alpha = get_symbol_vector("f", this->dim());

  for (unsigned int i = 0; i < this->num_of_templates();
       i++) { // for each parallelotope

    Parallelotope P = this->get_parallelotope(i);
    std::vector<Expression<>> genFun
        = build_instanciated_generator_functs(alpha, P);

    const std::vector<Expression<>> fog
        = replace_in(dynamics, variables, genFun);

    // compose sigma(f(gamma(x)))
    Expression<>::replacement_type repl;
    for (unsigned int j = 0; j < variables.size(); j++) {
      repl[variables[j]] = fog[j];
    }

    Expression<> sofog = atom->get_expression();
    sofog.replace(repl);

    // compute the Bernstein control points
    std::vector<Expression<>> controlPts
        = BaseConverter(alpha, sofog).getBernCoeffsMatrix();

    Polytope constraints(parameters, controlPts);
    result = ::intersect(result, constraints);
  }

  return result;
}

/**
 * @brief Transform the bundle according to a dynamic law
 *
 * @param[in] variables is the vector of the variables
 * @param[in] dynamics is the vector of the dynamic law expressions
 * @param[in] max_finder is a pointer to an MinMaxCoeffFinder object
 * @param[in] mode transformation mode, i.e., OFO or AFO
 * @returns the transformed bundle
 */
Bundle
Bundle::transform(const std::vector<SymbolicAlgebra::Symbol<>> &variables,
                  const std::vector<SymbolicAlgebra::Expression<>> &dynamics,
                  const Bundle::MinMaxCoeffFinder *max_finder,
                  Bundle::transfomation_mode mode) const
{

  using namespace std;
  using namespace SymbolicAlgebra;
  using namespace LinearAlgebra;

  vector<CondSyncUpdater<double, std::less<double>>> max_coeffs(this->size());
  vector<CondSyncUpdater<double, std::greater<double>>> min_coeffs(this->size());

  std::vector<Symbol<>> alpha = get_symbol_vector("f", dim());

  auto refine_coeff_itvl = [&max_coeffs, &min_coeffs, &variables, &alpha,
                            &dynamics, &max_finder,
                            &mode](const Bundle *bundle,
                                   const unsigned int template_num) {
    Parallelotope P = bundle->get_parallelotope(template_num);

    const std::vector<SymbolicAlgebra::Expression<>> &genFun
        = build_instanciated_generator_functs(alpha, P);
    const std::vector<SymbolicAlgebra::Expression<>> genFun_f
        = replace_in(dynamics, variables, genFun);

    const std::vector<int> &template_i = bundle->templates[template_num];

    unsigned int dir_b;

    // for each direction
    const size_t num_of_dirs
        = (mode == Bundle::OFO ? template_i.size()
                               : bundle->directions.size());

    for (unsigned int j = 0; j < num_of_dirs; j++) {
      if (mode == Bundle::OFO) {
        dir_b = template_i[j];
      } else {
        dir_b = j;
      }
      std::vector<SymbolicAlgebra::Expression<>> bernCoeffs
          = compute_Bern_coeffs(alpha, genFun_f, bundle->directions[dir_b]);

      auto coeff_itvl = max_finder->find_coeffs_itvl(bernCoeffs);

      min_coeffs[dir_b].update(coeff_itvl.first);
      max_coeffs[dir_b].update(coeff_itvl.second);
    }
  };

#ifdef WITH_THREADS
  ThreadPool::BatchId batch_id = thread_pool.create_batch();

  for (unsigned int i = 0; i < this->num_of_templates(); i++) {
    // submit the task to the thread pool
    thread_pool.submit_to_batch(batch_id, refine_coeff_itvl, this, i);
  }

  // join to the pool threads
  thread_pool.join_threads(batch_id);

  // close the batch
  thread_pool.close_batch(batch_id);
#else  // WITH_THREADS
  for (unsigned int i = 0; i < this->num_of_templates(); i++) {
    refine_coeff_itvl(this, i);
  }
#endif // WITH_THREADS

  LinearAlgebra::Vector<double> lower_bounds, upper_bounds;
  for (auto it = std::begin(max_coeffs); it != std::end(max_coeffs); ++it) {
    upper_bounds.push_back(*it);
  }
  for (auto it = std::begin(min_coeffs); it != std::end(min_coeffs); ++it) {
    lower_bounds.push_back(*it);
  }

  Bundle res(this->directions, lower_bounds, upper_bounds, this->templates);

  if (mode == Bundle::OFO) {
    res.canonize();
  }

  return res;
}

/**
 * Compute the distances between the half-spaced of the parallelotopes
 *
 * @returns vector of distances
 */
LinearAlgebra::Vector<double> Bundle::edge_lengths()
{
  using namespace LinearAlgebra;

  Vector<double> dist(this->size());
  for (unsigned int i = 0; i < this->size(); i++) {
    dist[i] = std::abs(this->upper_bounds[i] - this->lower_bounds[i])
              / norm_2(this->directions[i]);
  }
  return dist;
}

/**
 * Determine belonging of an element in a vector
 *
 * @param[in] n element to be searched
 * @param[in] v vector in which to look for
 * @returns true is n belongs to v
 */
bool isIn(int n, std::vector<int> v)
{

  for (unsigned int i = 0; i < v.size(); i++) {
    if (n == v[i]) {
      return true;
    }
  }
  return false;
}

/**
 * Determine belonging of a vector in a set of vectors
 *
 * @param[in] v vector to be searched
 * @param[in] vlist set of vectors in which to look for
 * @returns true is v belongs to vlist
 */
bool isIn(std::vector<int> v, std::vector<LinearAlgebra::Vector<int>> vlist)
{
  for (unsigned int i = 0; i < vlist.size(); i++) {
    if (isPermutation(v, vlist[i])) {
      return true;
    }
  }
  return false;
}

/**
 * @brief Search for linearly dependent vector in an array 
 * 
 * This function search for the index `i` such that the vector 
 * `A[i]` is linearly dependent to `v`. If `i` does exist, 
 * then this function returns it. Otherwise, it returns the 
 * `A`-invalid index `A.size()`.
 * 
 * @param v is the query vector
 * @param A is the array of vectors among which the linear 
 *        dependent vector is searched
 * @return if there exists a non-empty set of `i` such that
 *         `A[i]` and `v` are linear dependent, the minimum
 *         value in such a set. Otherwise, `A.size()`
 */
unsigned int
get_a_linearly_dependent_in(const LinearAlgebra::Vector<double> &v,
                            const std::vector<LinearAlgebra::Vector<double>> &A)
{
  for (unsigned int i = 0; i < A.size(); ++i) {
    if (LinearAlgebra::are_linearly_dependent(v, A[i])) {
      return i;
    }
  }
  return A.size();
}

/**
 * @brief Test whether a template direction have been placed before
 *        a index threshold  
 * 
 * During the intersection of two bundles we have to decide whether 
 * which templates are contained by both of them to avoid to add 
 * twice the same template to the resulting bundle. This is done 
 * by verifing that all the directions of the second bundle have 
 * been mapped to directions of the intersection bundle that lay 
 * before the number of directions of the first bundle.
 * 
 * @param bundle_template is the template of the second bundle 
 *                        whose presence in the first bundle must 
 *                        be tested
 * @param new_indices are the indices in the intersection 
 *                    bundle of the second bundle directions
 * @param first_bundle_size is the number of directions in the 
 *                          first bundle
 * @return `false` if and only if some of the indices of the 
 *         `bundle_template` directions have been mapped by 
 *         `new_indices` in an index greater or equal to 
 *         `first_bundle_size`
 */
bool copy_required(const LinearAlgebra::Vector<int> &bundle_template,
                   const std::vector<unsigned int> &new_indices,
                   const unsigned int first_bundle_size)
{
  for (auto &id: bundle_template) {
    if (new_indices[id] >= first_bundle_size) {
      return true;
    }
  }

  return false;
}

Bundle &Bundle::intersect(const Bundle &A)
{
  unsigned int old_size = this->size();
  std::vector<unsigned int> new_ids(A.size());

  // for each direction in A
  for (unsigned int i = 0; i < A.size(); ++i) {
    const LinearAlgebra::Vector<double> &A_dir = A.directions[i];
    new_ids[i] = get_a_linearly_dependent_in(A_dir, this->directions);

    // if the direction is not present in this object
    if (new_ids[i] == this->size()) {

      // add the direction and the corresponding boundaries
      this->directions.push_back(A_dir);
      this->lower_bounds.push_back(A.lower_bounds[i]);
      this->upper_bounds.push_back(A.upper_bounds[i]);
    } else { // if the direction is already included in this object

      // compute the dependency coefficent
      const unsigned int &new_i = new_ids[i];
      const double dep_coeff = this->directions[new_i][0] / A_dir[0];

      // if necessary, increase the lower bound
      this->lower_bounds[new_i]
          = std::max(this->lower_bounds[new_i], A.lower_bounds[i] * dep_coeff);

      // if necessary, decrease the upper bound
      this->upper_bounds[new_i]
          = std::min(this->upper_bounds[new_i], A.upper_bounds[i] * dep_coeff);
    }
  }

  // check whether some of the templates must be copied
  for (unsigned int i = 0; i < A.templates.size(); ++i) {
    const std::vector<int> &A_template = A.templates[i];
    if (copy_required(A_template, new_ids, old_size)) {

      // if this is the case, copy and update the template indices
      std::vector<int> t_copy(A_template);
      for (unsigned int j = 0; j < t_copy.size(); ++j) {
        t_copy[j] = new_ids[t_copy[j]];
      }

      // add the new template to the intesected bundle
      this->templates.push_back(t_copy);
    }
  }

  return *this;
}

void Bundle::add_templates_for(
    std::set<unsigned int> &to_be_copied_directions)
{
  using namespace LinearAlgebra;
  using namespace LinearAlgebra::Dense;

  while (!to_be_copied_directions.empty()) {
    // fill T with all the bundle directions starting
    // from those not included in a template
    Matrix<double> T = this->directions;
    Vector<int> row_pos(T.size());
    std::iota(std::begin(row_pos), std::end(row_pos), 0);

    unsigned int j = 0;
    for (unsigned int i: to_be_copied_directions) {
      std::swap(T[i], T[j]);
      std::swap(row_pos[i], row_pos[j]);
      ++j;
    }

    // the LUP factorization find the first n non-linearly
    // dependent rows in T and the permutation can be used
    // to discover them
    auto fP = LUP_Factorization<double>(T).P();
    row_pos = fP(row_pos);
    Vector<int> new_template(this->dim());

    std::copy(std::begin(row_pos), std::begin(row_pos) + new_template.size(),
              std::begin(new_template));

    this->templates.push_back(new_template);
    for (int i: new_template) {
      to_be_copied_directions.erase((unsigned int)i);
    }
  }
}

/**
 * @brief Get the intersection between this bundle and a linear set
 *
 * This method intersects the current instance of the `Bundle` class
 * and a possible unbounded linear set provided as a `LinearSystem`.
 * The result is stored in the current object and a reference to it
 * is returned.
 *
 * @param ls is the intersecting linear system
 * @return a reference to the updated object
 */
Bundle &Bundle::intersect(const LinearSystem &ls)
{
  using namespace LinearAlgebra;
  using namespace LinearAlgebra::Dense;

  std::set<unsigned int> outside_templates;
  const Matrix<double> &A = ls.A();

  std::vector<unsigned int> new_ids(A.size());
  // for each row in the linear system
  for (unsigned int i = 0; i < A.size(); ++i) {
    const LinearAlgebra::Vector<double> &A_row = A[i];

    new_ids[i] = get_a_linearly_dependent_in(A_row, this->directions);

    // if the direction is not present in this object
    if (new_ids[i] == this->size()) {
      // compute a non-influent lower bound for the row
      double lower_bound = ((Polytope)*this).minimize(A_row).optimum();

      // add the direction and the corresponding boundaries
      outside_templates.insert(new_ids[i]);
      this->directions.push_back(A_row);
      this->lower_bounds.push_back(lower_bound);
      this->upper_bounds.push_back(ls.getb(i));
    } else {
      // compute the dependency coefficent
      const unsigned int &new_i = new_ids[i];
      const double dep_coeff = this->directions[new_i] / A_row;

      auto b_rescaled = ls.getb(i) * dep_coeff;
      if (dep_coeff>0) {
        // if necessary, decrease the upper bound
        if (this->upper_bounds[new_i] > b_rescaled) {
          this->upper_bounds[new_i] = b_rescaled;
        }
      } else {
        // or increase the lower bound
        if (this->lower_bounds[new_i] < b_rescaled) {
          this->lower_bounds[new_i] = b_rescaled;
        }
      }
    }
  }

  // add missing templates
  add_templates_for(outside_templates);

  return *this;
}

Bundle::~Bundle()
{
  // TODO Auto-generated destructor stub
}
