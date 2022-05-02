/**
 * @file Bundle.cpp
 * Represent and manipulate bundles of parallelotopes whose intersection
 * represents a polytope
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
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

#define _USE_MATH_DEFINES

#include <cmath>

#define AVOID_NEG_ZERO(value) ((value) == 0 ? 0 : (value))

/**
 * Copy constructor that instantiates the bundle
 *
 * @param[in] orig is the model for the new bundle
 */
Bundle::Bundle(const Bundle &orig):
    directions(orig.directions), lower_bounds(orig.lower_bounds),
    upper_bounds(orig.upper_bounds), templates(orig.templates),
    constraintDirections(orig.constraintDirections),
    constraintOffsets(orig.constraintOffsets)
{
}

/**
 * Move constructor that instantiates the bundle
 *
 * @param[in, out] orig is the model for the new bundle
 */
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
  std::swap(A.constraintDirections, B.constraintDirections);
  std::swap(A.constraintOffsets, B.constraintOffsets);
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

/**
 * Constructor
 *
 * @param[in] directions matrix of directions
 * @param[in] lower_bounds lower offsets
 * @param[in] upper_bounds upper offsets
 * @param[in] templates templates matrix
 */
Bundle::Bundle(const std::vector<LinearAlgebra::Vector<double>> &directions,
               const LinearAlgebra::Vector<double> &lower_bounds,
               const LinearAlgebra::Vector<double> &upper_bounds,
               const std::vector<LinearAlgebra::Vector<int>> &templates):
    Bundle(directions, lower_bounds, upper_bounds, templates, {}, {})
{
}

/**
 * Move constructor
 *
 * @param[in] directions matrix of directions
 * @param[in] lower_bounds lower offsets
 * @param[in] upper_bounds upper offsets
 * @param[in] templates templates matrix
 */
Bundle::Bundle(std::vector<LinearAlgebra::Vector<double>> &&directions,
               LinearAlgebra::Vector<double> &&lower_bounds, LinearAlgebra::Vector<double> &&upper_bounds,
               std::vector<LinearAlgebra::Vector<int>> &&templates):
    directions(std::move(directions)),
    lower_bounds(std::move(lower_bounds)),
    upper_bounds(std::move(upper_bounds)), templates(std::move(templates))
{
}

/**
 * Constructor
 *
 * @param[in] directions matrix of directions
 * @param[in] lower_bounds lower offsets
 * @param[in] upper_bounds upper offsets
 * @param[in] templates templatess matrix
 * @param[in] constrDirs directions that are constrained by assumptions
 * @param[in] constrOffsets offsets of assumptions
 */
Bundle::Bundle(const std::vector<LinearAlgebra::Vector<double>> &directions,
               const LinearAlgebra::Vector<double> &lower_bounds,
               const LinearAlgebra::Vector<double> &upper_bounds,
               const std::vector<LinearAlgebra::Vector<int>> &templates,
               const std::vector<LinearAlgebra::Vector<double>> &constrDirs,
               const LinearAlgebra::Vector<double> &constrOffsets):
    directions(directions),
    lower_bounds(lower_bounds), upper_bounds(upper_bounds),
    templates(templates), constraintDirections(constrDirs),
    constraintOffsets(constrOffsets)
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

/**
 * Get the i-th parallelotope of the bundle
 *
 * @param[in] i parallelotope index to fetch
 * @returns i-th parallelotope
 */
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

/**
 * Canonize the current bundle pushing the constraints toward the symbolic
 * polytope
 *
 * @returns canonized bundle
 */
Bundle Bundle::get_canonical() const
{
  // get current polytope
  Polytope bund = *this;
  LinearAlgebra::Vector<double> c_up_offset(this->size()), c_lo_offset(this->size());
  for (unsigned int i = 0; i < this->size(); ++i) {
    c_up_offset[i] = bund.maximize(this->directions[i]).optimum();
    c_lo_offset[i] = bund.minimize(this->directions[i]).optimum();
  }
  return Bundle(this->directions, c_lo_offset, c_up_offset, templates);
}

/**
 * Check whether a vector is permutation of a sorted vector
 *
 * @param[in] v1 first vector
 * @param[in] v2 a sorted vector
 * @returns true is v1 is a permutation of v2
 */
bool isPermutationOfSorted(std::vector<int> v1,
                           const std::vector<int> &v2_sorted)
{
  if (v1.size() != v2_sorted.size()) {
    return false;
  }

  std::sort(std::begin(v1), std::end(v1));

  auto v1_it = std::begin(v1);
  auto v2_it = std::begin(v2_sorted);
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

std::list<Bundle> &
split_bundle(std::list<Bundle> &res, LinearAlgebra::Vector<double> &lower_bounds,
             LinearAlgebra::Vector<double> &upper_bounds, const unsigned int idx,
             const Bundle &splitting, const double &max_bundle_magnitude,
             const double &split_magnitude_ratio, int)
{
  if (idx == splitting.get_directions().size()) {
    res.emplace_back(splitting.get_directions(), lower_bounds, upper_bounds,
                     splitting.get_templates());

    return res;
  }

  if (std::abs(splitting.get_upper_bound(idx) - splitting.get_lower_bound(idx))
      > max_bundle_magnitude) {
    double lower_bound = splitting.get_lower_bound(idx);

    do {
      const double upper_bound = std::min(
          lower_bound + split_magnitude_ratio * max_bundle_magnitude,
          splitting.get_upper_bound(idx));

      lower_bounds[idx] = lower_bound;
      upper_bounds[idx] = upper_bound;
      split_bundle(res, lower_bounds, upper_bounds, idx + 1, splitting,
                   max_bundle_magnitude, split_magnitude_ratio, 2);

      lower_bound = upper_bound;
    } while (splitting.get_upper_bound(idx) != lower_bound);
  } else {
    lower_bounds[idx] = splitting.get_lower_bound(idx);
    upper_bounds[idx] = splitting.get_upper_bound(idx);
    split_bundle(res, lower_bounds, upper_bounds, idx + 1, splitting,
                 max_bundle_magnitude, split_magnitude_ratio, 3);
  }
  return res;
}

std::list<Bundle> Bundle::split(const double max_bundle_magnitude,
                                const double split_magnitude_ratio) const
{
  std::list<Bundle> split_list;

  LinearAlgebra::Vector<double> upper_bounds(this->size());
  LinearAlgebra::Vector<double> lower_bounds(this->size());

  split_bundle(split_list, lower_bounds, upper_bounds, 0, *this,
               max_bundle_magnitude, split_magnitude_ratio, 3);

  return split_list;
}

// TODO: the following method probably does not work; it
//       should be fixed
/**
 * Decompose the current symbolic polytope
 *
 * @param[in] dec_weight weight parameter in [0,1] for decomposition (0 for
 * distance, 1 for orthogonality)
 * @param[in] max_iter maximum number of randomly generated templates
 * @returns new bundle decomposing current symbolic polytope
 */
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

std::vector<SymbolicAlgebra::Expression<>>
sub_vars(const std::vector<SymbolicAlgebra::Expression<>> &ex_list,
         const std::vector<SymbolicAlgebra::Symbol<>> &vars,
         const std::vector<SymbolicAlgebra::Expression<>> &expressions)
{
  using namespace SymbolicAlgebra;

  Expression<>::replacement_type repl;

  for (unsigned int k = 0; k < vars.size(); ++k) {
    repl[vars[k]] = expressions[k];
  }

  std::vector<Expression<>> results;
  for (auto ex_it = std::begin(ex_list); ex_it != std::end(ex_list); ++ex_it) {
    results.push_back(Expression<>(*ex_it).replace(repl));
  }

  return results;
}

std::vector<SymbolicAlgebra::Expression<>>
compute_Bern_coeffs(const std::vector<SymbolicAlgebra::Symbol<>> &alpha,
                    const std::vector<SymbolicAlgebra::Expression<>> &f,
                    const LinearAlgebra::Vector<double> &dir_vector)
{
  SymbolicAlgebra::Expression<> Lfog = 0;
  // upper facets
  for (unsigned int k = 0; k < dir_vector.size(); k++) {
    if (dir_vector[k] != 0) {
      Lfog += dir_vector[k] * f[k];
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
 * $$
 * q + ((alpha \circ beta)^T \cdot G)^T
 * $$
 * where $\cdot$ is the row-column product, $\circ$ is the
 * Hadamard product, and $q$, $beta$, and $G$ are the base vector,
 * the vector lengths, the versors matrix of the
 * considered parallelotope $P$, respectively.
 *
 * @param alpha is a vector of variables.
 * @param P is the considered parallelotope.
 * @return The generator function of `P`.
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

  const std::vector<Vector<double>> &versors = P.versors();

  for (unsigned int i = 0; i < versors.size(); i++) {
    // some of the non-null rows of the generator matrix
    // correspond to 0-length dimensions in degenerous
    // parallelotopes and must be avoided
    if (P.lengths()[i] != 0) {
      Vector<double> vector = P.lengths()[i] * versors[i];
      for (unsigned int j = 0; j < vector.size(); j++) {
        if (vector[j] != 0) {
          gen_functs[j] += alpha[i] * vector[j];
        }
      }
    }
  }

  return gen_functs;
}

template<typename COND>
class CondSyncUpdate
{
#ifdef WITH_THREADS
  mutable std::shared_timed_mutex mutex;
#endif
  double _value;
  COND _cmp;

public:
  CondSyncUpdate()
  {
    // initialize _value to the minimum wrt the order
    _value = (_cmp(std::numeric_limits<double>::lowest(),
                   std::numeric_limits<double>::max())
                  ? std::numeric_limits<double>::max()
                  : std::numeric_limits<double>::lowest());
  }

  inline operator double() const
  {
#ifdef WITH_THREADS
    std::shared_lock<std::shared_timed_mutex> readlock(mutex);
#endif
    return _value;
  }

  void update(const double &value)
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
 * @param[in] pSet is the initial parameter set
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
        = sub_vars(dynamics, variables, genFun);

    // compose sigma(f(gamma(x)))
    Expression<>::replacement_type repl;
    for (unsigned int j = 0; j < variables.size(); j++) {
      repl[variables[j]] = fog[j];
    }

    Expression<> sofog = atom->getPredicate();
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

  vector<CondSyncUpdate<std::less<double>>> max_coeffs(this->size());
  vector<CondSyncUpdate<std::greater<double>>> min_coeffs(this->size());

  std::vector<Symbol<>> alpha = get_symbol_vector("f", dim());

  auto refine_coeff_itvl = [&max_coeffs, &min_coeffs, &variables, &alpha,
                            &dynamics, &max_finder,
                            &mode](const Bundle *bundle,
                                   const unsigned int template_num) {
    Parallelotope P = bundle->get_parallelotope(template_num);

    const std::vector<SymbolicAlgebra::Expression<>> &genFun
        = build_instanciated_generator_functs(alpha, P);
    const std::vector<SymbolicAlgebra::Expression<>> genFun_f
        = sub_vars(dynamics, variables, genFun);

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

      // for each asserted direction, check that the new offset
      // does not violate the constraint
      for (unsigned assertIndex = 0;
           assertIndex < bundle->constraintDirections.size(); assertIndex++) {
        if (bundle->directions[dir_b]
            == bundle->constraintDirections[assertIndex]) {
          max_coeffs[dir_b].update(bundle->constraintOffsets[assertIndex]);
        } else if (bundle->directions[dir_b]
                   == -bundle->constraintDirections[assertIndex]) {
          min_coeffs[dir_b].update(-bundle->constraintOffsets[assertIndex]);
        }
      }
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

  Bundle res
      = Bundle(this->directions, lower_bounds, upper_bounds, this->templates,
               this->constraintDirections, this->constraintOffsets);

  if (mode == Bundle::OFO) {
    return res.get_canonical();
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

LinearSystem get_linear_system(const Bundle &A)
{
  using namespace LinearAlgebra;

  std::vector<Vector<double>> A_dirs = A.get_directions();
  A_dirs.reserve(2 * A_dirs.size());
  for (unsigned int i = 0; i < A_dirs.size(); ++i) {
    A_dirs.push_back(-A_dirs[i]);
  }

  Vector<double> A_bounds = A.get_upper_bounds();
  A_bounds.reserve(2 * A_bounds.size());
  for (unsigned int i = 0; i < A_bounds.size(); ++i) {
    A_bounds.push_back(AVOID_NEG_ZERO(-A.get_lower_bound(i)));
  }

  return LinearSystem(std::move(A_dirs), std::move(A_bounds));
}

unsigned int
get_a_linearly_dependent_row_in(const LinearAlgebra::Vector<double> &v,
                                const std::vector<LinearAlgebra::Vector<double>> &A)
{
  for (unsigned int i = 0; i < A.size(); ++i) {
    if (LinearAlgebra::are_linearly_dependent(v, A[i])) {
      return i;
    }
  }
  return A.size();
}

bool require_copy(const LinearAlgebra::Vector<int> &bundle_template,
                  const std::vector<unsigned int> &new_ids,
                  const unsigned int old_size)
{
  for (auto &id: bundle_template) {
    if (new_ids[id] >= old_size) {
      return true;
    }
  }

  return false;
}

/**
 * @brief Get the intersection between two bundles
 *
 * This method intersects the current object and another bundle.
 * The result is stored in the current object and a reference
 * to it is returned.
 *
 * @param A is the intersecting bundle.
 * @return a reference to the updated object.
 */
Bundle &Bundle::intersect(const Bundle &A)
{
  unsigned int old_size = this->size();
  std::vector<unsigned int> new_ids(A.size());

  // for each direction in A
  for (unsigned int i = 0; i < A.size(); ++i) {
    const LinearAlgebra::Vector<double> &A_dir = A.directions[i];
    new_ids[i] = get_a_linearly_dependent_row_in(A_dir, this->directions);

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

  // copy missing templates
  for (unsigned int i = 0; i < A.templates.size(); ++i) {
    const std::vector<int> &A_template = A.templates[i];
    if (require_copy(A_template, new_ids, old_size)) {
      std::vector<int> t_copy(A_template);
      for (unsigned int j = 0; j < t_copy.size(); ++j) {
        t_copy[j] = new_ids[t_copy[j]];
      }
      this->templates.push_back(t_copy);
    }
  }

  return *this;
}

/**
 * @brief Add to the bundle templates for some directions
 *
 * @param missing_template_dirs the indices of the directions whose
 *            template we are missing
 */
void Bundle::add_templates_for(
    std::set<unsigned int> &missing_template_directions)
{
  using namespace LinearAlgebra;
  using namespace LinearAlgebra::Dense;

  while (!missing_template_directions.empty()) {
    // fill T with all the bundle directions starting
    // from those not included in a template
    Matrix<double> T = this->directions;
    Vector<int> row_pos;
    std::iota(std::begin(row_pos), std::end(row_pos), 1);

    unsigned int j = 0;
    for (unsigned int i: missing_template_directions) {
      std::swap(T[i], T[j]);
      std::swap(row_pos[i], row_pos[j]);
      ++j;
    }

    // the LUP factorization find the first n non-linearly
    // dependent rows in T and the permutation can be used
    // to discover them
    auto fP = LUP_Factorization<double>(T).P();
    row_pos = fP(row_pos);
    Vector<int> new_temp(this->dim());

    std::copy(std::begin(row_pos), std::begin(row_pos) + new_temp.size(),
              std::begin(new_temp));

    this->templates.push_back(new_temp);
    for (int i: new_temp) {
      missing_template_directions.erase((unsigned int)i);
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
 * @param ls is the intersecting bundle.
 * @return a reference to the updated object.
 */
Bundle &Bundle::intersect(const LinearSystem &ls)
{
  using namespace LinearAlgebra;
  using namespace LinearAlgebra::Dense;

  std::set<unsigned int> outside_templates;
  const Matrix<double> &A = ls.getA();

  std::vector<unsigned int> new_ids(A.size());
  // for each row in the linear system
  for (unsigned int i = 0; i < A.size(); ++i) {
    const LinearAlgebra::Vector<double> &A_row = A[i];

    new_ids[i] = get_a_linearly_dependent_row_in(A_row, this->directions);

    // if the direction is not present in this object
    if (new_ids[i] == this->size()) {
      // compute a non-influent lower bound for the row
      LinearSystem this_ls = get_linear_system(*this);

      double lower_bound = this_ls.minimize(A_row).optimum();

      // add the direction and the corresponding boundaries
      outside_templates.insert(this->directions.size());
      this->directions.push_back(A_row);
      this->lower_bounds.push_back(lower_bound);
      this->upper_bounds.push_back(ls.getb(i));
    } else {
      // compute the dependency coefficent
      const unsigned int &new_i = new_ids[i];
      const double dep_coeff = this->directions[new_i] / A_row;

      // if necessary, decrease the upper bound
      this->upper_bounds[new_i]
          = std::min(this->upper_bounds[new_i], ls.getb(i) * dep_coeff);
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
