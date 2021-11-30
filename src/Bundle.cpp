/**
 * @file Bundle.cpp
 * Represent and manipulate bundles of parallelotopes whose intersection
 * represents a polytope
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#include "Bundle.h"


#include <limits>
#include <string>
#include <algorithm>

#include "LinearAlgebra.h"

#define _USE_MATH_DEFINES

#include <cmath>

/**
 * Copy constructor that instantiates the bundle
 *
 * @param[in] orig is the model for the new bundle
 */
Bundle::Bundle(const Bundle &orig):
    dir_matrix(orig.dir_matrix), offp(orig.offp), offm(orig.offm),
    t_matrix(orig.t_matrix), Theta(orig.Theta)
{
}

/**
 * Swap constructor that instantiates the bundle
 *
 * @param[in] orig is the model for the new bundle
 */
Bundle::Bundle(Bundle &&orig)
{
  swap(*this, orig);
}

void swap(Bundle &A, Bundle &B)
{
  std::swap(A.dir_matrix, B.dir_matrix);
  std::swap(A.offp, B.offp);
  std::swap(A.offm, B.offm);
  std::swap(A.t_matrix, B.t_matrix);
  std::swap(A.Theta, B.Theta);
}

/**
 * Orthogonal proximity of v1 and v2, i.e.,
 * how close is the angle between v1 and v2 is to pi/2
 *
 * @param[in] v1 vector
 * @param[in] v2 vector
 * @returns orthogonal proximity
 */
double orthProx(std::vector<double> v1, std::vector<double> v2)
{
  return std::abs(angle(v1, v2) - M_PI_2);
}

/**
 * Constructor that instantiates the bundle with auto-generated variables
 *
 * @param[in] dir_matrix matrix of directions
 * @param[in] offp upper offsets
 * @param[in] offm lower offsets
 * @param[in] t_matrix templates matrix
 */
Bundle::Bundle(const Matrix &dir_matrix, const Vector &offp,
               const Vector &offm,
               const std::vector<std::vector<int>> &t_matrix):
    dir_matrix(dir_matrix),
    offp(offp), offm(offm), t_matrix(t_matrix)
{
  using namespace std;
  using namespace GiNaC;

  if (dir_matrix.size() == 0) {
    std::cerr << "Bundle::Bundle : dir_matrix must be non empty" << std::endl;

    exit(EXIT_FAILURE);
  }
  if (dir_matrix.size() != offp.size()) {
    std::cerr << "Bundle::Bundle : dir_matrix and offp "
              << "must have the same size" << std::endl;
    exit(EXIT_FAILURE);
  }
  if (dir_matrix.size() != offm.size()) {
    std::cerr << "Bundle::Bundle : dir_matrix and offm must have "
              << "the same size" << std::endl;
    exit(EXIT_FAILURE);
  }
  if (t_matrix.size() > 0) {
    for (unsigned int i = 0; i < t_matrix.size(); i++) {
      if (t_matrix[i].size() != this->dim()) {
        std::cerr << "Bundle::Bundle : t_matrix must have " << this->dim()
                  << " columns" << std::endl;
        exit(EXIT_FAILURE);
      }
    }
  } else {
    std::cerr << "Bundle::Bundle : t_matrix must be non empty" << std::endl;
    exit(EXIT_FAILURE);
  }

  // initialize orthogonal proximity
  this->Theta = vector<vector<double>>(this->size(),
                                       vector<double>(this->size(), 0));

  for (unsigned int i = 0; i < this->size(); i++) {
    this->Theta[i][i] = 0;
    for (unsigned int j = i + 1; j < this->size(); j++) {
      double prox = orthProx(this->dir_matrix[i], this->dir_matrix[j]);
      this->Theta[i][j] = prox;
      this->Theta[j][i] = prox;
    }
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

  vector<vector<double>> A;
  vector<double> b;
  for (unsigned int i = 0; i < this->size(); i++) {
    A.push_back(this->dir_matrix[i]);
    b.push_back(this->offp[i]);
  }
  for (unsigned int i = 0; i < this->size(); i++) {
    A.push_back(-this->dir_matrix[i]);
    b.push_back(this->offm[i]);
  }

  return Polytope(A, b);
}

/**
 * Get the i-th parallelotope of the bundle
 *
 * @param[in] i parallelotope index to fetch
 * @returns i-th parallelotope
 */
Parallelotope Bundle::getParallelotope(unsigned int i) const
{
  using namespace std;

  if (i > this->t_matrix.size()) {
    cerr << "Bundle::getParallelotope : i must be between 0 and "
         << t_matrix.size() << endl;
    exit(EXIT_FAILURE);
  }

  vector<double> lbound, ubound;
  vector<vector<double>> Lambda;

  vector<int>::const_iterator it = std::begin(this->t_matrix[i]);
  // upper facets
  for (unsigned int j = 0; j < this->dim(); j++) {
    const int idx = *it;
    Lambda.push_back(this->dir_matrix[idx]);
    ubound.push_back(this->offp[idx]);
    lbound.push_back(this->offm[idx]);

    ++it;
  }

  // TODO: Since Lambdas are always the same, check whether
  //       storing the PLU factorizations of Lambdas may give
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
  std::vector<double> canoffp(this->size()), canoffm(this->size());
  for (unsigned int i = 0; i < this->size(); i++) {
    canoffp[i] = bund.maximize(this->dir_matrix[i]);
    canoffm[i] = bund.maximize(-this->dir_matrix[i]);
  }
  return Bundle(dir_matrix, canoffp, canoffm, t_matrix);
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
                     const std::vector<double> &dists)
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
                     const std::vector<double> &dists)
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
double maxOffsetDist(const std::vector<std::vector<int>> &T,
                     const std::vector<double> &dists)
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
 * @param[in] dir_matrix is the direction matrix
 * @param[in] vIdx index of the reference vector
 * @param[in] dirsIdx indexes of vectors to be considered
 * @returns maximum orthogonal proximity
 */
double maxOrthProx(const std::vector<std::vector<double>> &dir_matrix,
                   const int vIdx, const std::vector<int> &dirsIdx)
{

  if (dirsIdx.empty()) {
    return 0;
  }

  double maxProx = 0;
  for (auto d_it = std::begin(dirsIdx); d_it != std::end(dirsIdx); ++d_it) {
    maxProx = std::max(maxProx, orthProx(dir_matrix[vIdx], dir_matrix[*d_it]));
  }
  return maxProx;
}

/**
 * Maximum orthogonal proximity within a set of vectors
 *
 * @param[in] dir_matrix is the direction matrix
 * @param[in] dirsIdx indexes of vectors to be considered
 * @returns maximum orthogonal proximity
 */
double maxOrthProx(const std::vector<std::vector<double>> &dir_matrix,
                   const std::vector<int> &dirsIdx)
{
  double maxProx = 0;
  for (unsigned int i = 0; i < dirsIdx.size(); i++) {
    for (unsigned int j = i + 1; j < dirsIdx.size(); j++) {
      maxProx = std::max(
          maxProx, orthProx(dir_matrix[dirsIdx[i]], dir_matrix[dirsIdx[j]]));
    }
  }
  return maxProx;
}

/**
 * Maximum orthogonal proximity of all the vectors of a matrix
 *
 * @param[in] dir_matrix is the direction matrix
 * @param[in] T collection of vectors
 * @returns maximum orthogonal proximity
 */
double maxOrthProx(const std::vector<std::vector<double>> &dir_matrix,
                   const std::vector<std::vector<int>> &T)
{
  double maxorth = std::numeric_limits<double>::lowest();
  for (auto T_it = std::begin(T); T_it != std::end(T); ++T_it) {
    maxorth = std::max(maxorth, maxOrthProx(dir_matrix, *T_it));
  }
  return maxorth;
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

  vector<double> offDists = this->offsetDistances();

  // get current template and try to improve it
  vector<vector<int>> curT = this->t_matrix;

  // get current template and try to improve it
  vector<vector<int>> bestT = this->t_matrix;
  int temp_card = this->t_matrix.size();

  int i = 0;
  while (i < max_iters) {

    vector<vector<int>> tmpT = curT;

    // generate random coordinates to swap
    unsigned int i1 = rand() % temp_card;
    int j1 = rand() % this->dim();

    // swap them
    tmpT[i1][j1] = rand() % this->size();

    if (!is_permutation_of_other_rows(tmpT, i1)) {
      std::vector<std::vector<double>> A;
      for (unsigned int j = 0; j < this->dim(); j++) {
        A.push_back(this->dir_matrix[tmpT[i1][j]]);
      }

      DenseLinearAlgebra::PLU_Factorization<double> fact(A);
      try {
        fact.solve(std::vector<double>(this->dim(), 0));

        double w1 = dec_weight * maxOffsetDist(tmpT, offDists)
                    + (1 - dec_weight) * maxOrthProx(this->dir_matrix, tmpT);
        double w2 = dec_weight * maxOffsetDist(bestT, offDists)
                    + (1 - dec_weight) * maxOrthProx(this->dir_matrix, bestT);

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

  return Bundle(dir_matrix, offp, offm, bestT);
}

GiNaC::lst sub_vars(const GiNaC::lst &ex_list, const GiNaC::lst &vars,
                    const GiNaC::lst &expressions)
{
  GiNaC::lst sub;

  for (unsigned int k = 0; k < vars.nops(); ++k) {
    sub.append(vars[k] == expressions[k]);
  }

  GiNaC::lst results;
  for (auto ex_it = std::begin(ex_list); ex_it != std::end(ex_list); ++ex_it) {
    results.append(ex_it->subs(sub));
  }

  return results;
}

GiNaC::lst compute_Bern_coeffs(const GiNaC::lst &alpha, const GiNaC::lst &f,
                               const std::vector<double> &dir_vector)
{
  GiNaC::ex Lfog = 0;
  // upper facets
  for (unsigned int k = 0; k < dir_vector.size(); k++) {
    if (dir_vector[k] != 0) {
      Lfog = Lfog + dir_vector[k] * f[k];
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
GiNaC::lst get_subs_from(const Parallelotope &P, const GiNaC::lst &q,
                         const GiNaC::lst &beta)
{
  const std::vector<double> &base_vertex = P.base_vertex();
  const std::vector<double> &lengths = P.lengths();

  GiNaC::lst subs;

  for (unsigned int k = 0; k < q.nops(); k++) {
    subs.append(q[k] == base_vertex[k]);
    subs.append(beta[k] == lengths[k]);
  }

  return subs;
}

double Bundle::MaxCoeffFinder::coeff_eval_p(const GiNaC::ex &c) const
{
  using namespace GiNaC;

  return ex_to<numeric>(c).to_double();
}

double Bundle::MaxCoeffFinder::coeff_eval_m(const GiNaC::ex &bernCoeff) const
{
  using namespace GiNaC;

  double value = ex_to<numeric>(bernCoeff).to_double();

  // TODO: The following conditional evaluation avoids -0
  //       values. Check the difference between -0 and 0.
  return (value == 0 ? 0 : -value);
}

double
Bundle::ParamMaxCoeffFinder::coeff_eval_p(const GiNaC::ex &bernCoeff) const
{
  return paraSet.maximize(params, bernCoeff);
}

double
Bundle::ParamMaxCoeffFinder::coeff_eval_m(const GiNaC::ex &bernCoeff) const
{
  return paraSet.maximize(params, -bernCoeff);
};

Bundle::MaxCoeffFinder::MaxCoeffType
Bundle::MaxCoeffFinder::find_max_coeffs(const GiNaC::lst &b_coeffs) const
{
  // find the maximum coefficient
  GiNaC::lst::const_iterator b_coeff_it = b_coeffs.begin();

  double maxCoeffp = coeff_eval_p(*b_coeff_it);
  double maxCoeffm = coeff_eval_m(*b_coeff_it);

  for (++b_coeff_it; b_coeff_it != b_coeffs.end(); ++b_coeff_it) {
    double actCoeff = coeff_eval_p(*b_coeff_it);

    if (actCoeff > maxCoeffp) {
      maxCoeffp = actCoeff;
    }

    actCoeff = coeff_eval_m(*b_coeff_it);
    if (actCoeff > maxCoeffm) {
      maxCoeffm = actCoeff;
    }
  }

  return MaxCoeffType{maxCoeffp, maxCoeffm};
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
GiNaC::lst build_instanciated_generator_functs(const GiNaC::lst &alpha,
                                               const Parallelotope &P)
{
  GiNaC::lst gen_functs;
  for (auto it = std::begin(P.base_vertex()); it != std::end(P.base_vertex());
       ++it) {
    gen_functs.append(*it);
  }

  const std::vector<std::vector<double>> &versors = P.versors();

  for (unsigned int i = 0; i < versors.size(); i++) {
    // some of the non-null rows of the generator matrix
    // correspond to 0-length dimensions in degenerous
    // parallelotopes and must be avoided
    if (P.lengths()[i] != 0) {
      std::vector<double> vector = P.lengths()[i] * versors[i];
      for (unsigned int j = 0; j < vector.size(); j++) {
        gen_functs[j] = gen_functs[j] + alpha[i] * vector[j];
      }
    }
  }

  return gen_functs;
}

/**
 * Transform the bundle
 *
 * @param[in] vars variables appearing in the transforming function
 * @param[in] f transforming function
 * @param[in] max_finder is a pointer to an MaxCoeffFinder object.
 * @param[in] mode transformation mode (0=OFO,1=AFO)
 * @returns transformed bundle
 */
Bundle Bundle::transform(const GiNaC::lst &vars, const GiNaC::lst &f,
                         const Bundle::MaxCoeffFinder *max_finder,
                         int mode) const
{
  using namespace std;
  using namespace GiNaC;

  vector<double> newDp(this->size(), std::numeric_limits<double>::max());
  vector<double> newDm = newDp;

  GiNaC::lst alpha = get_symbol_lst("f", dim());

  vector<int> dirs_to_bound;
  if (mode) { // dynamic transformation
    for (unsigned int i = 0; i < this->dir_matrix.size(); i++) {
      dirs_to_bound.push_back(i);
    }
  }

  // for each parallelotope
  for (unsigned int i = 0; i < this->num_of_templates(); i++) {

    Parallelotope P = this->getParallelotope(i);
  
    const lst &genFun = build_instanciated_generator_functs(alpha, P);
    const lst genFun_f = sub_vars(f, vars, genFun);

    if (mode == 0) { // static mode
      dirs_to_bound = this->t_matrix[i];
    }

    // for each direction
    for (unsigned int j = 0; j < dirs_to_bound.size(); j++) {
      lst bernCoeffs
          = compute_Bern_coeffs(alpha, genFun_f, dir_matrix[dirs_to_bound[j]]);

      auto maxCoeff = max_finder->find_max_coeffs(bernCoeffs);

      const unsigned int &dir_b = dirs_to_bound[j];
      if (newDp[dir_b] > maxCoeff.p) {
        newDp[dir_b] = maxCoeff.p;
      }

      if (newDm[dir_b] > maxCoeff.m) {
        newDm[dir_b] = maxCoeff.m;
      }
    }
  }

  Bundle res = Bundle(this->dir_matrix, newDp, newDm, this->t_matrix);
  if (mode == 0) {
    return res.get_canonical();
  }

  return res;
}

/**
 * Compute the distances between the half-spaced of the parallelotopes
 *
 * @returns vector of distances
 */
std::vector<double> Bundle::offsetDistances()
{

  std::vector<double> dist(this->size());
  for (unsigned int i = 0; i < this->size(); i++) {
    dist[i] = std::abs(this->offp[i] - this->offm[i])
              / norm_2(this->dir_matrix[i]);
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
bool isIn(std::vector<int> v, std::vector<std::vector<int>> vlist)
{
  for (unsigned int i = 0; i < vlist.size(); i++) {
    if (isPermutation(v, vlist[i])) {
      return true;
    }
  }
  return false;
}

Bundle::~Bundle()
{
  // TODO Auto-generated destructor stub
}
