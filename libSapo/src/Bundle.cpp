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

#include <limits>
#include <string>
#include <algorithm>
#include <functional>
#include <sstream>

#include <glpk.h>

#define _USE_MATH_DEFINES //!< This macro enables the use of cmath constants

#include <cmath>

#ifdef WITH_THREADS
#include <mutex>
#include <shared_mutex>

#include "SapoThreads.h"
#endif // WITH_THREADS

#include "SetsUnion.h"
#include "LinearAlgebra.h"
#include "LinearAlgebraIO.h"

/**
 * @brief Avoid \f$-0\f$
 *
 * This macro avoids \f$-0\f$ by replacing it by \f$0\f$.
 */
#define AVOID_NEG_ZERO(value) ((value) == 0 ? 0 : (value))

Bundle::Bundle() {}

Bundle::Bundle(const Bundle &orig):
    _directions(orig._directions), _lower_bounds(orig._lower_bounds),
    _upper_bounds(orig._upper_bounds), _templates(orig._templates)
{
}

Bundle::Bundle(Bundle &&orig)
{
  swap(*this, orig);
}

void swap(Bundle &A, Bundle &B)
{
  std::swap(A._directions, B._directions);
  std::swap(A._upper_bounds, B._upper_bounds);
  std::swap(A._lower_bounds, B._lower_bounds);
  std::swap(A._templates, B._templates);
}

/**
 * Orthogonal proximity of v1 and v2, i.e.,
 * how close is the angle between v1 and v2 is to pi/2
 *
 * @param[in] v1 vector
 * @param[in] v2 vector
 * @returns orthogonal proximity
 */
double orthProx(LinearAlgebra::Vector<double> v1,
                LinearAlgebra::Vector<double> v2)
{
  return std::abs(LinearAlgebra::angle(v1, v2) - M_PI_2);
}

/**
 * @brief Test whether a vector is sorted
 *
 * @tparam T is the type of the values
 * @param V is the vector to test
 * @return `true` if and only if `V` is sorted
 */
template<typename T>
bool is_sorted(const std::vector<T> &V)
{
  typename std::vector<T>::const_iterator v_it = std::begin(V);
  typename std::vector<T>::const_iterator previous_it = std::begin(V);

  for (++v_it; v_it != std::end(V); ++previous_it, ++v_it) {
    if (*previous_it > *v_it) {
      return false;
    }
  }

  return true;
}

Bundle::Bundle(const std::vector<LinearAlgebra::Vector<double>> &directions,
               const LinearAlgebra::Vector<double> &lower_bounds,
               const LinearAlgebra::Vector<double> &upper_bounds,
               const std::set<std::vector<unsigned int>> &templates):
    _directions(directions),
    _lower_bounds(lower_bounds), _upper_bounds(upper_bounds), _templates()
{
  using namespace LinearAlgebra;

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
  if (templates.size() == 0) {
    throw std::domain_error("Bundle::Bundle: templates must be non empty");
  }

  for (auto t_it = std::begin(templates); t_it != std::end(templates);
       ++t_it) {
    if (t_it->size() != this->dim()) {
      throw std::domain_error("Bundle::Bundle: templates must have "
                              "as many columns as directions");
    }

    std::vector<unsigned int> new_template(*t_it);

    if (!is_sorted(new_template)) {
      // sort unsorted templates
      std::sort(std::begin(new_template), std::end(new_template));
    }

    _templates.insert(new_template);

    Dense::Matrix<double> A;
    for (auto d_it = std::begin(*t_it); d_it != std::end(*t_it); ++d_it) {
      if (*d_it >= directions.size()) {
        throw std::domain_error("Bundle::Bundle: templates must contains "
                                "as values indices of the directions vector");
      }

      A.emplace_back(directions[*d_it]);
    }

    if (Dense::rank(A) != A.size()) {
      std::ostringstream oss;

      oss << "Bundle::Bundle: template directions must be linearly "
             "independent, "
          << "but template directions " << *t_it << " are not.";
      throw std::domain_error(oss.str());
    }
  }
}

Bundle::Bundle(std::vector<LinearAlgebra::Vector<double>> &&directions,
               LinearAlgebra::Vector<double> &&lower_bounds,
               LinearAlgebra::Vector<double> &&upper_bounds,
               const std::set<std::vector<unsigned int>> &templates):
    _directions(std::move(directions)),
    _lower_bounds(std::move(lower_bounds)),
    _upper_bounds(std::move(upper_bounds)), _templates()
{
  using namespace std;
  using namespace LinearAlgebra;

  if (_directions.size() == 0) {
    throw std::domain_error("Bundle::Bundle: directions must be non empty");
  }
  if (_directions.size() != _upper_bounds.size()) {
    throw std::domain_error("Bundle::Bundle: directions and upper_bounds "
                            "must have the same size");
  }
  if (_directions.size() != _lower_bounds.size()) {
    throw std::domain_error("Bundle::Bundle: directions and lower_bounds "
                            "must have the same size");
  }
  if (templates.size() == 0) {
    throw std::domain_error("Bundle::Bundle: templates must be non empty");
  }

  for (auto t_it = std::begin(templates); t_it != std::end(templates);
       ++t_it) {
    if (t_it->size() != this->dim()) {
      throw std::domain_error("Bundle::Bundle: templates must have "
                              "as many columns as directions");
    }

    std::vector<unsigned int> new_template(*t_it);

    if (!is_sorted(new_template)) {
      // sort unsorted templates
      std::sort(std::begin(new_template), std::end(new_template));
    }

    _templates.insert(new_template);

    Dense::Matrix<double> A;
    for (auto d_it = std::begin(*t_it); d_it != std::end(*t_it); ++d_it) {
      if (*d_it >= _directions.size()) {
        throw std::domain_error("Bundle::Bundle: templates must contains "
                                "as values indices of the directions vector");
      }

      A.emplace_back(_directions[*d_it]);
    }

    if (Dense::rank(A) != A.size()) {
      std::ostringstream oss;

      oss << "Bundle::Bundle: template directions must be linearly "
             "independent, "
          << "but template directions " << *t_it << " are not.";
      throw std::domain_error(oss.str());
    }
  }
}

Bundle::Bundle(const std::vector<LinearAlgebra::Vector<double>> &directions,
               const LinearAlgebra::Vector<double> &lower_bounds,
               const LinearAlgebra::Vector<double> &upper_bounds):
    _directions(directions),
    _lower_bounds(lower_bounds), _upper_bounds(upper_bounds), _templates()
{
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

  // add missing templates
  std::set<unsigned int> missing;

  for (unsigned int i = 0; i < _directions.size(); ++i) {
    missing.insert(i);
  }

  add_templates_for(missing);
}

Bundle::Bundle(std::vector<LinearAlgebra::Vector<double>> &&directions,
               LinearAlgebra::Vector<double> &&lower_bounds,
               LinearAlgebra::Vector<double> &&upper_bounds):
    _directions(std::move(directions)),
    _lower_bounds(std::move(lower_bounds)),
    _upper_bounds(std::move(upper_bounds)), _templates()
{
  using namespace std;

  if (_directions.size() == 0) {
    throw std::domain_error("Bundle::Bundle: directions must be non empty");
  }
  if (_directions.size() != _upper_bounds.size()) {
    throw std::domain_error("Bundle::Bundle: directions and upper_bounds "
                            "must have the same size");
  }
  if (_directions.size() != _lower_bounds.size()) {
    throw std::domain_error("Bundle::Bundle: directions and lower_bounds "
                            "must have the same size");
  }

  // add missing templates
  std::set<unsigned int> missing;

  for (unsigned int i = 0; i < _directions.size(); ++i) {
    missing.insert(i);
  }

  add_templates_for(missing);
}

Bundle::Bundle(const Polytope &P)
{
  using namespace LinearAlgebra;

  std::vector<size_t> double_row(P.size(), false);

  for (size_t i = 0; i < P.size(); ++i) {
    if (!double_row[i]) {
      const auto &Ai = P.A(i);
      const auto nAi = -P.A(i);
      for (size_t j = i + 1; j < P.size(); ++j) {
        if (!double_row[j] && (P.A(j) == Ai || P.A(j) == nAi)) {
          double_row[j] = true;
        }
      }
      _directions.push_back(P.A(i));
      _upper_bounds.push_back(P.maximize(P.A(i)).optimum());
      _lower_bounds.push_back(P.minimize(P.A(i)).optimum());
    }
  }

  // add missing template
  std::set<unsigned int> missing;

  for (unsigned int i = 0; i < _directions.size(); ++i) {
    missing.insert(i);
  }

  add_templates_for(missing);
}

/**
 * @brief Copy operator
 *
 * @param orig is the original object to be copied
 * @return a reference to the updated object
 */
Bundle &Bundle::operator=(const Bundle &orig)
{
  using namespace LinearAlgebra;

  this->_directions = std::vector<Vector<double>>();
  this->_templates = std::set<Vector<unsigned int>>();

  this->_directions.reserve(orig._directions.size());

  for (auto d_it = std::begin(orig._directions);
       d_it != std::end(orig._directions); ++d_it) {
    this->_directions.emplace_back(*d_it);
  }

  for (auto t_it = std::begin(orig._templates);
       t_it != std::end(orig._templates); ++t_it) {
    this->_templates.emplace(*t_it);
  }

  this->_upper_bounds = orig._upper_bounds;
  this->_lower_bounds = orig._lower_bounds;

  return *this;
}

/**
 * @brief Copy operator
 *
 * @param orig is the original object to be copied
 * @return a reference to the updated object
 */
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
    A.push_back(this->_directions[i]);
    b.push_back(this->_upper_bounds[i]);
  }
  for (unsigned int i = 0; i < this->size(); i++) {
    A.push_back(-this->_directions[i]);
    b.push_back(AVOID_NEG_ZERO(-this->_lower_bounds[i]));
  }

  return Polytope(std::move(A), std::move(b));
}

Parallelotope Bundle::get_parallelotope(
    const std::vector<unsigned int> &bundle_template) const
{
  using namespace std;

  vector<double> lbound, ubound;
  vector<LinearAlgebra::Vector<double>> Lambda;

  vector<unsigned int>::const_iterator it = std::begin(bundle_template);
  // upper facets
  for (unsigned int j = 0; j < this->dim(); j++) {
    const unsigned int idx = *it;
    if (idx >= this->_directions.size()) {
      throw std::domain_error(
          "Bundle::get_parallelotope: the parameter is not a template for "
          "for the bundle");
    }
    Lambda.push_back(this->_directions[idx]);
    ubound.push_back(this->_upper_bounds[idx]);
    lbound.push_back(this->_lower_bounds[idx]);

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

Bundle &Bundle::canonize()
{
  if (this->size() == 0) {
    return *this;
  }

  // if the bundle is empty
  if (is_empty()) {
    return *this;
  }

  // get current polytope
  Polytope bund = *this;
  for (unsigned int i = 0; i < this->size(); ++i) {
    _lower_bounds[i] = bund.minimize(this->_directions[i]).optimum();
    _upper_bounds[i] = bund.maximize(this->_directions[i]).optimum();
  }
  return *this;
}

/// @private
void _bundle_free_lp_problem(glp_prob *lp, int *ia, int *ja, double *ar)
{
  glp_delete_prob(lp);
  glp_free_env();
  free(ia);
  free(ja);
  free(ar);
}

bool are_disjoint(const Bundle &A, const Bundle &B)
{
  if (A.dim() != B.dim()) {
    throw std::domain_error("The two bundles must have the "
                            "same dimension");
  }

  if (A.dim() == 0) {
    return false;
  }

  unsigned int num_rows = A.size() + B.size();
  unsigned int num_cols = A.dim();
  unsigned int size_lp = num_rows * num_cols;

  std::vector<double> obj_fun(num_cols, 0);
  obj_fun[0] = 1;

  int *ia, *ja;
  double *ar;

  ia = (int *)calloc(size_lp + 1, sizeof(int));
  ja = (int *)calloc(size_lp + 1, sizeof(int));
  ar = (double *)calloc(size_lp + 1, sizeof(double));

  glp_prob *lp;
  lp = glp_create_prob();
  glp_set_obj_dir(lp, GLP_MAX);

  // Turn off verbose mode
  glp_smcp lp_param;
  glp_init_smcp(&lp_param);
  lp_param.msg_lev = GLP_MSG_ERR;

  glp_add_rows(lp, num_rows);
  for (unsigned int i = 0; i < A.size(); i++) {
    if (A.get_lower_bound(i) > A.get_upper_bound(i)) {
      _bundle_free_lp_problem(lp, ia, ja, ar);

      return true;
    }
    glp_set_row_bnds(lp, i + 1, GLP_DB, A.get_lower_bound(i),
                     A.get_upper_bound(i));
  }
  for (unsigned int i = 0; i < B.size(); i++) {
    if (B.get_lower_bound(i) > B.get_upper_bound(i)) {
      _bundle_free_lp_problem(lp, ia, ja, ar);

      return true;
    }
    glp_set_row_bnds(lp, A.size() + i + 1, GLP_DB, B.get_lower_bound(i),
                     B.get_upper_bound(i));
  }

  glp_add_cols(lp, num_cols);
  for (unsigned int i = 0; i < num_cols; i++) {
    glp_set_col_bnds(lp, i + 1, GLP_FR, 0.0, 0.0);
    glp_set_obj_coef(lp, i + 1, obj_fun[i]);
  }

  unsigned int k = 1;
  for (unsigned int i = 0; i < A.size(); i++) {
    for (unsigned int j = 0; j < num_cols; j++) {
      ia[k] = i + 1;
      ja[k] = j + 1;
      ar[k] = A.get_direction(i)[j]; /* a[i+1,j+1] = A[i][j] */
      k++;
    }
  }
  for (unsigned int i = 0; i < B.size(); i++) {
    for (unsigned int j = 0; j < num_cols; j++) {
      ia[k] = A.size() + i + 1;
      ja[k] = j + 1;
      ar[k] = B.get_direction(i)[j]; /* a[i+1,j+1] = A[i][j] */
      k++;
    }
  }

  glp_load_matrix(lp, size_lp, ia, ja, ar);
  glp_exact(lp, &lp_param);

  auto status = glp_get_status(lp);
  bool res = (status == GLP_NOFEAS || status == GLP_INFEAS);

  _bundle_free_lp_problem(lp, ia, ja, ar);

  return res;
}

/**
 * @brief Test whether a bundle is subset of another bundle
 *
 * This method tests whether the current object is subset
 * of a bundle.
 *
 * @param[in] bundle is the tested bundle
 * @return `true` if and only if the current bundle is a
 *         subset of `bundle`
 */
bool Bundle::is_subset_of(const Bundle &bundle) const
{
  // if this object is empty
  if (this->is_empty()) {
    return true;
  }

  // if the parameter is empty and this object is not,
  // return false
  if (bundle.is_empty()) {
    return false;
  }

  Polytope P_this = *this;
  // for each direction in the bundle
  for (unsigned int dir_idx = 0; dir_idx < bundle.size(); ++dir_idx) {

    // if the minimum of this object on that direction is lesser than
    // the bundle minimum, this object is not a subset of the bundle
    if (P_this.minimize(bundle.get_direction(dir_idx)).optimum()
        < bundle.get_lower_bound(dir_idx)) {
      return false;
    }

    // if the maximum of this object on that direction is greater than
    // the bundle maximum, this object is not a subset of the bundle
    if (P_this.maximize(bundle.get_direction(dir_idx)).optimum()
        > bundle.get_upper_bound(dir_idx)) {
      return false;
    }
  }

  // if none of the previous conditions hold, this object is a
  // subset for the bundle
  return true;
}

/**
 * @brief Check whether a bundle satisfies a linear system
 *
 * This method checks whether all the points in the
 * current object are solutions for a linear system.
 *
 * @param ls is the considered linear system
 * @return `true` if and only if all the points of
 *          the current bundle are solutions for `ls`
 */
bool Bundle::satisfies(const LinearSystem &ls) const
{
  // if this object is empty
  if (this->is_empty()) {
    return true;
  }

  // if the parameter has no solutions and this object is not,
  // return false
  if (!ls.has_solutions()) {
    return false;
  }

  Polytope P_this = *this;
  // for each direction in the bundle
  for (unsigned int dir_idx = 0; dir_idx < ls.size(); ++dir_idx) {

    // if the maximum of this object on that direction is smaller than
    // the bundle maximum, this object does not include the bundle
    if (P_this.maximize(ls.A(dir_idx)).optimum() > ls.b(dir_idx)) {
      return false;
    }
  }

  // if none of the previous conditions hold, the bundle is a
  // subset for this object
  return true;
}

/**
 * Check whether a vector is permutation of a sorted vector
 *
 * @param[in] v1 first vector
 * @param[in] sorted a sorted vector
 * @returns true is v1 is a permutation of v2
 */
bool isPermutationOfSorted(std::vector<unsigned int> v1,
                           const std::vector<unsigned int> &sorted)
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
bool isPermutation(const std::vector<unsigned int> &v1,
                   std::vector<unsigned int> v2)
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
 * @param[in] distances pre-computed distances
 * @returns distance accumulation
 */
double maxOffsetDist(const int vIdx, const std::vector<int> &dirsIdx,
                     const LinearAlgebra::Vector<double> &distances)
{

  if (dirsIdx.empty()) {
    return 0;
  }

  double dist = distances[vIdx];
  for (unsigned int i = 0; i < dirsIdx.size(); i++) {
    dist = dist * distances[dirsIdx[i]];
  }
  return dist;
}

/**
 * Maximum distance accumulation of a set of vectors
 *
 * @param[in] dirsIdx indexes of vectors to be considered
 * @param[in] distances pre-computed distances
 * @returns distance accumulation
 */
double maxOffsetDist(const std::vector<unsigned int> &dirsIdx,
                     const LinearAlgebra::Vector<double> &distances)
{

  double dist = 1;
  for (unsigned int i = 0; i < dirsIdx.size(); i++) {
    dist = dist * distances[dirsIdx[i]];
  }
  return dist;
}

/**
 * Maximum distance accumulation of matrix
 *
 * @param[in] T matrix from which fetch the vectors
 * @param[in] distances pre-computed distances
 * @returns distance accumulation
 */
double maxOffsetDist(const std::vector<std::vector<unsigned int>> &T,
                     const LinearAlgebra::Vector<double> &distances)
{
  double max_dist = std::numeric_limits<double>::lowest();
  for (unsigned int i = 0; i < T.size(); i++) {
    max_dist = std::max(max_dist, maxOffsetDist(T[i], distances));
  }
  return max_dist;
}

/**
 * Maximum orthogonal proximity of a vector w.r.t. a set of vectors
 *
 * @param[in] directions is the direction matrix
 * @param[in] vIdx index of the reference vector
 * @param[in] dirsIdx indexes of vectors to be considered
 * @returns maximum orthogonal proximity
 */
double
maxOrthProx(const std::vector<LinearAlgebra::Vector<double>> &directions,
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
double
maxOrthProx(const std::vector<LinearAlgebra::Vector<double>> &directions,
            const std::vector<unsigned int> &dirsIdx)
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
double
maxOrthProx(const std::vector<LinearAlgebra::Vector<double>> &directions,
            const std::vector<std::vector<unsigned int>> &T)
{
  double max_orth = std::numeric_limits<double>::lowest();
  for (auto T_it = std::begin(T); T_it != std::end(T); ++T_it) {
    max_orth = std::max(max_orth, maxOrthProx(directions, *T_it));
  }
  return max_orth;
}

/**
 * @brief A draft horse function to split a bundle
 *
 * This recursive function is a private draft horse to
 * split a bundle whose maximal magnitude, the maximal
 * length of its generators, is greater than
 * `max_magnitude`. The bundle is split into a
 * list of sub-bundles whose maximal magnitude is
 * \f$\textrm{max_magnitude}*\textrm{split_ratio}\f$.
 * Each sub-bundle in the output list is built recursively
 * by considering one bundle direction per time and:
 * 1. splitting the direction in the opportune number of
 *    chunks, so to satisfy the sub-bundle maximal
 *    magnitude request;
 * 2. for each of the chunks, setting the sub-bundle
 *    boundaries of the selected direction according
 *    with the considered chunk and recursively
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
std::list<Bundle> &split_bundle(std::list<Bundle> &res,
                                LinearAlgebra::Vector<double> &lower_bounds,
                                LinearAlgebra::Vector<double> &upper_bounds,
                                const unsigned int idx,
                                const Bundle &splitting,
                                const double &max_magnitude,
                                const double &split_ratio)
{
  if (idx == splitting.directions().size()) {
    res.emplace_back(splitting.directions(), lower_bounds, upper_bounds,
                     splitting.templates());

    return res;
  }

  if (std::abs(splitting.get_upper_bound(idx) - splitting.get_lower_bound(idx))
      > max_magnitude) {
    double lower_bound = splitting.get_lower_bound(idx);

    do {
      const double upper_bound
          = std::min(lower_bound + split_ratio * max_magnitude,
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

  split_bundle(split_list, lower_bounds, upper_bounds, 0, *this, max_magnitude,
               split_ratio);

  return split_list;
}

Bundle Bundle::decompose(double dec_weight, int max_iters)
{
  using namespace std;
  using namespace LinearAlgebra;

  vector<double> off_distances = this->edge_lengths();

  // get current template and try to improve it
  vector<Vector<unsigned int>> curT(std::begin(this->_templates),
                                    std::end(this->_templates));

  // get current template and try to improve it
  vector<Vector<unsigned int>> bestT(std::begin(this->_templates),
                                     std::end(this->_templates));
  int temp_card = this->_templates.size();

  int i = 0;
  while (i < max_iters) {

    vector<Vector<unsigned int>> tmpT = curT;

    // generate random coordinates to swap
    unsigned int i1 = rand() % temp_card;
    int j1 = rand() % this->dim();

    // swap them
    tmpT[i1][j1] = rand() % this->size();

    if (!is_permutation_of_other_rows(tmpT, i1)) {
      std::vector<Vector<double>> A;
      for (unsigned int j = 0; j < this->dim(); j++) {
        A.push_back(this->_directions[tmpT[i1][j]]);
      }

      Dense::LUP_Factorization<double> fact(A);
      try {
        fact.solve(Vector<double>(this->dim(), 0));

        double w1 = dec_weight * maxOffsetDist(tmpT, off_distances)
                    + (1 - dec_weight) * maxOrthProx(this->_directions, tmpT);
        double w2 = dec_weight * maxOffsetDist(bestT, off_distances)
                    + (1 - dec_weight) * maxOrthProx(this->_directions, bestT);

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

  return Bundle(_directions, _lower_bounds, _upper_bounds,
                set<Vector<unsigned int>>(std::begin(bestT), std::end(bestT)));
}

/**
 * Compute the distances between the half-spaced of the parallelotopes
 *
 * @returns vector of distances
 */
LinearAlgebra::Vector<double> Bundle::edge_lengths() const
{
  using namespace LinearAlgebra;

  Vector<double> dist(this->size());
  for (unsigned int i = 0; i < this->size(); i++) {
    dist[i] = std::abs(this->_upper_bounds[i] - this->_lower_bounds[i])
              / norm_2(this->_directions[i]);
  }
  return dist;
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
unsigned int get_a_linearly_dependent_in(
    const LinearAlgebra::Vector<double> &v,
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
 * by verifying that all the directions of the second bundle have
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
bool copy_required(const std::vector<unsigned int> &bundle_template,
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

/**
 * @brief Add missing templates to a template vector
 *
 * This function is called as the last phase of intersection or
 * approximated union between two bundles. After identifying
 * shared directions and adding non-shared ones to the direction
 * vector, the templates involving non-shared
 *
 * @param dest_templates is the destination template set
 * @param source_templates is the source template set
 * @param new_direction_indices is the map from source indices
 *                              to destination indices
 */
void add_missing_templates(
    std::set<std::vector<unsigned int>> &dest_templates,
    const std::set<std::vector<unsigned int>> &source_templates,
    const std::vector<unsigned int> &new_direction_indices)
{
  // check whether some of the templates must be copied
  for (const std::vector<unsigned int> &temp: source_templates) {
    // if this is the case, copy and update the template indices
    std::vector<unsigned int> t_copy(temp);
    for (unsigned int j = 0; j < t_copy.size(); ++j) {
      t_copy[j] = new_direction_indices[t_copy[j]];
    }

    std::sort(std::begin(t_copy), std::end(t_copy));

    // add the new template to the intersected bundle
    dest_templates.insert(t_copy);
  }
}

Bundle &Bundle::intersect_with(const Bundle &A)
{
  using namespace LinearAlgebra;

  if (dim() != A.dim()) {
    throw std::domain_error("Bundle::intersect_with: the two "
                            "bundles differ in dimensions");
  }

  std::vector<unsigned int> new_ids(A.size());

  // for each direction in A
  for (unsigned int i = 0; i < A.size(); ++i) {
    const Vector<double> &A_dir = A._directions[i];
    new_ids[i] = get_a_linearly_dependent_in(A_dir, this->_directions);

    // if the direction is not present in this object
    if (new_ids[i] == this->size()) {

      // add the direction and the corresponding boundaries
      this->_directions.push_back(A_dir);
      this->_lower_bounds.push_back(A._lower_bounds[i]);
      this->_upper_bounds.push_back(A._upper_bounds[i]);
    } else { // if the direction is already included in this object

      // compute the dependency coefficient
      const unsigned int &new_i = new_ids[i];
      const double dep_coeff = this->_directions[new_i] / A_dir;

      // if necessary, increase the lower bound
      this->_lower_bounds[new_i] = std::max(this->_lower_bounds[new_i],
                                            A._lower_bounds[i] * dep_coeff);

      // if necessary, decrease the upper bound
      this->_upper_bounds[new_i] = std::min(this->_upper_bounds[new_i],
                                            A._upper_bounds[i] * dep_coeff);
    }
  }

  add_missing_templates(_templates, A.templates(), new_ids);

  return *this;
}

void Bundle::add_templates_for(std::set<unsigned int> &to_be_copied_directions)
{
  using namespace LinearAlgebra;
  using namespace LinearAlgebra::Dense;

  while (!to_be_copied_directions.empty()) {
    // fill T with all the bundle directions starting
    // from those not included in a template
    Matrix<double> T = this->_directions;
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
    Vector<unsigned int> new_template(this->dim());

    std::copy(std::begin(row_pos), std::begin(row_pos) + new_template.size(),
              std::begin(new_template));

    if (!is_sorted(new_template)) {
      // sort unsorted templates
      std::sort(std::begin(new_template), std::end(new_template));
    }

    this->_templates.insert(new_template);
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
Bundle &Bundle::intersect_with(const LinearSystem &ls)
{
  using namespace LinearAlgebra;
  using namespace LinearAlgebra::Dense;

  // If the linear system has no constraint or involves no variables
  if (ls.dim() == 0) {
    return *this;
  }

  if (dim() != ls.dim()) {
    throw std::domain_error("Bundle::intersect_with: the bundle "
                            "and the linear system differ in "
                            "dimensions");
  }

  std::set<unsigned int> outside_templates;
  const Matrix<double> &A = ls.A();

  std::vector<unsigned int> new_ids(A.size());
  // for each row in the linear system
  for (unsigned int i = 0; i < A.size(); ++i) {
    const LinearAlgebra::Vector<double> &A_row = A[i];

    new_ids[i] = get_a_linearly_dependent_in(A_row, this->_directions);

    // if the direction is not present in this object
    if (new_ids[i] == this->size()) {
      // compute a non-influent lower bound for the row
      double lower_bound = ((Polytope) * this).minimize(A_row).optimum();

      // add the direction and the corresponding boundaries
      outside_templates.insert(new_ids[i]);
      this->_directions.push_back(A_row);
      this->_lower_bounds.push_back(lower_bound);
      this->_upper_bounds.push_back(ls.b(i));
    } else {
      // compute the dependency coefficient
      const unsigned int &new_i = new_ids[i];
      const double dep_coeff = this->_directions[new_i] / A_row;

      auto b_rescaled = ls.b(i) * dep_coeff;
      if (dep_coeff > 0) {
        // if necessary, decrease the upper bound
        if (this->_upper_bounds[new_i] > b_rescaled) {
          this->_upper_bounds[new_i] = b_rescaled;
        }
      } else {
        // or increase the lower bound
        if (this->_lower_bounds[new_i] < b_rescaled) {
          this->_lower_bounds[new_i] = b_rescaled;
        }
      }
    }
  }

  // add missing templates
  add_templates_for(outside_templates);

  return *this;
}

Bundle &Bundle::expand_by(const double delta)
{
  if (is_empty()) {
    return *this;
  }

  for (auto b_it = std::begin(_lower_bounds); b_it != std::end(_lower_bounds);
       ++b_it) {
    *b_it -= delta;
  }

  for (auto b_it = std::begin(_upper_bounds); b_it != std::end(_upper_bounds);
       ++b_it) {
    *b_it += delta;
  }

  return *this;
}

Bundle::~Bundle() {}

Bundle over_approximate_union(const Bundle &b1, const Bundle &b2)
{
  using namespace LinearAlgebra;
  using namespace LinearAlgebra::Dense;

  if (b1.dim() != b2.dim()) {
    throw std::domain_error("over_approximate_union: the two "
                            "bundles differ in dimensions");
  }

  if (b1.is_empty()) {
    return b2;
  }

  if (b2.is_empty()) {
    return b1;
  }

  Polytope p2(b2);
  Bundle res(b1);
  const Matrix<double> &res_dirs = res.directions();

  // Updates res boundaries to include p2
  for (unsigned int i = 0; i < res_dirs.size(); ++i) {
    const LinearAlgebra::Vector<double> &res_dir = res_dirs[i];
    res._lower_bounds[i]
        = std::min(p2.minimize(res_dir).optimum(), res.get_lower_bound(i));
    res._upper_bounds[i]
        = std::max(p2.maximize(res_dir).optimum(), res.get_upper_bound(i));
  }

  const Matrix<double> &b2_dirs = b2.directions();

  Polytope p1(b1);
  std::vector<unsigned int> new_ids(b2.size());
  // for each row in the linear system
  for (unsigned int i = 0; i < b2_dirs.size(); ++i) {
    const LinearAlgebra::Vector<double> &b2_dir = b2_dirs[i];

    new_ids[i] = get_a_linearly_dependent_in(b2_dir, res.directions());

    // if the direction is not present in this object
    if (new_ids[i] == res.size()) {
      double lower_bound
          = std::min(p1.minimize(b2_dir).optimum(), b2.get_lower_bound(i));
      double upper_bound
          = std::max(p1.maximize(b2_dir).optimum(), b2.get_upper_bound(i));

      // add the direction and the corresponding boundaries
      res._directions.push_back(b2_dir);
      res._lower_bounds.push_back(lower_bound);
      res._upper_bounds.push_back(upper_bound);
    } else {
      // compute the dependency coefficient
      const unsigned int &new_i = new_ids[i];
      const double dep_coeff = res.get_direction(new_i) / b2_dir;

      double lb_rescaled = b2.get_lower_bound(i) * dep_coeff;
      double ub_rescaled = b2.get_upper_bound(i) * dep_coeff;
      if (dep_coeff < 0) {
        std::swap(lb_rescaled, ub_rescaled);
      }

      // if necessary, update the boundaries
      if (res.get_upper_bound(new_i) < ub_rescaled) {
        res._upper_bounds[new_i] = ub_rescaled;
      }
      if (res.get_lower_bound(new_i) > lb_rescaled) {
        res._lower_bounds[new_i] = lb_rescaled;
      }
    }
  }

  add_missing_templates(res._templates, b2.templates(), new_ids);

  return res;
}

SetsUnion<Bundle> subtract_and_close(const Bundle &b1, const Bundle &b2)
{
  SetsUnion<Bundle> su;
  if (b2.includes(b1)) {
    return su;
  }

  if (are_disjoint(b1, b2)) {
    su.add(b1);

    return su;
  }

  Polytope p1 = b1;

  for (unsigned int i = 0; i < b2.size(); ++i) {
    auto new_bound = p1.maximize(b2.get_direction(i)).optimum();
    if (new_bound > b2.get_upper_bound(i)) {
      Bundle new_b1 = b1;
      new_b1._directions.push_back(b2.get_direction(i));
      new_b1._lower_bounds.push_back(b2.get_upper_bound(i));
      new_b1._upper_bounds.push_back(new_bound);

      su.add(std::move(new_b1.canonize()));
    }

    new_bound = p1.minimize(b2.get_direction(i)).optimum();
    if (new_bound < b2.get_lower_bound(i)) {
      Bundle new_b1 = b1;
      new_b1._directions.push_back(b2.get_direction(i));
      new_b1._lower_bounds.push_back(new_bound);
      new_b1._upper_bounds.push_back(b2.get_lower_bound(i));

      su.add(std::move(new_b1.canonize()));
    }
  }

  return su;
}
