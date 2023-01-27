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
#include <algorithm> //!< std::sort
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

#include "ErrorHandling.h"

/**
 * @brief Avoid \f$-0\f$
 *
 * This macro avoids \f$-0\f$ by replacing it by \f$0\f$.
 */
#define AVOID_NEG_ZERO(value) ((value) == 0 ? 0 : (value))

BundleTemplate::BundleTemplate(const BundleTemplate &orig):
    _dir_indices(orig._dir_indices),
    _dynamic_directions(orig._dynamic_directions)
{
}

BundleTemplate::BundleTemplate(BundleTemplate &&orig):
    BundleTemplate(std::move(orig._dir_indices), orig._dynamic_directions)
{
}

BundleTemplate::BundleTemplate(const std::vector<unsigned int> &dir_indices,
                               const bool dynamic_directions):
    _dir_indices(dir_indices),
    _dynamic_directions(dynamic_directions)
{
  std::sort(std::begin(_dir_indices), std::end(_dir_indices));
}

BundleTemplate::BundleTemplate(std::vector<unsigned int> &&dir_indices,
                               const bool dynamic_directions):
    _dir_indices(std::move(dir_indices)),
    _dynamic_directions(dynamic_directions)
{
  std::sort(std::begin(_dir_indices), std::end(_dir_indices));
}

bool std::less<BundleTemplate>::operator()(const BundleTemplate &a,
                                           const BundleTemplate &b) const
{
  std::less<std::vector<unsigned int>> cmp;

  return cmp(a.get_direction_indices(), b.get_direction_indices());
}

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

/**
 * @brief Canonize a template
 *
 * A template is in canonical form if it is sorted. This function
 * sorts the direction indices of a template and brings it in
 * canonical form.
 *
 * @param bundle_template is the template to be brought in canonical form
 * @return a reference to the updated template
 */
std::vector<unsigned int> &
canonize_template(std::vector<unsigned int> &bundle_template)
{
  if (!is_sorted(bundle_template)) {
    std::sort(std::begin(bundle_template), std::end(bundle_template));
  }

  return bundle_template;
}

/**
 * @brief Canonize a set of templates
 *
 * A template is in canonical form if it is sorted. This function
 * builds a set of templates in canonical form that is equivalent to
 * the set in input.
 *
 * @param templates is the set of templates to be brought in canonical form
 * @return the set of templates brought in canonical form
 */
std::set<std::vector<unsigned int>>
canonize_templates(const std::set<std::vector<unsigned int>> &templates)
{
  std::set<std::vector<unsigned int>> new_templates;

  for (auto bundle_template: templates) {
    canonize_template(bundle_template);

    new_templates.insert(std::move(bundle_template));
  }

  return new_templates;
}

/**
 * @brief Collect directions in templates
 *
 * @param[out] template_directions is the set in which the template directions
 * are placed
 * @param[in] templates is the set of templates whose directions is aimed for
 */
void collect_directions_in_templates(
    std::set<size_t> &template_directions,
    const std::set<std::vector<unsigned int>> &templates)
{
  for (const auto &bundle_template: templates) {
    template_directions.insert(std::begin(bundle_template),
                               std::end(bundle_template));
  }
}

/**
 * @brief Collect directions in templates
 *
 * @param templates is the set of templates whose directions is aimed for
 * @return the set of directions mentioned in the templates
 */
std::set<size_t> collect_directions_in_templates(
    const std::set<std::vector<unsigned int>> &templates)
{
  std::set<size_t> dirs;

  collect_directions_in_templates(dirs, templates);

  return dirs;
}

/**
 * @brief Find new positions for all the positions in a set
 *
 * @param positions is the set of position which must be replaced
 * @return a map assigning a new position to each of the old
 *        positions
 */
std::map<size_t, size_t>
find_new_position_for(const std::set<size_t> &positions)
{
  std::map<size_t, size_t> new_positions;
  size_t new_id = 0;

  // for each of the positions in the set of the positions to
  // be reassigned
  for (const auto &dir_id: positions) {

    // assign a new position
    new_positions[dir_id] = new_id++;
  }

  return new_positions;
}

template<typename T>
void mark_linearly_dep_dirs(
    std::vector<size_t> &new_pos,
    const std::vector<LinearAlgebra::Vector<T>> &directions,
    const size_t &dir_idx)
{
  using namespace LinearAlgebra;

  const size_t new_idx = new_pos[dir_idx];

  const auto &dir = directions[dir_idx];
  for (size_t j = dir_idx + 1; j < directions.size(); ++j) {
    if (new_pos[j] == j && are_linearly_dependent(dir, directions[j])) {
      new_pos[j] = new_idx;
    }
  }
}

/**
 * @brief Filter duplicate directions
 *
 * @tparam T is the scalar type of the directions
 * @param directions is the direction vector
 * @param lower_bounds is the vector of the direction lower bounds
 * @param upper_bounds is the vector of the direction upper bounds
 * @return a vector mapping each index of the input direction vector in the
 *       corresponding index in the final direction vector
 */
template<typename T>
std::vector<size_t>
filter_duplicated_directions(std::vector<LinearAlgebra::Vector<T>> &directions,
                             LinearAlgebra::Vector<T> &lower_bounds,
                             LinearAlgebra::Vector<T> &upper_bounds)
{
  std::vector<LinearAlgebra::Vector<T>> new_directions;
  LinearAlgebra::Vector<T> new_lower_bounds, new_upper_bounds;

  std::vector<size_t> new_pos(directions.size());
  std::iota(std::begin(new_pos), std::end(new_pos), 0);

  for (size_t i = 0; i < directions.size(); ++i) {

    // if directions[i] is the first instance of this direction
    // in directions
    if (new_pos[i] == i) {

      // get the new index for positions[i] in new_directions
      // and set in new_pos[i]
      const size_t new_idx = new_directions.size();
      new_pos[i] = new_idx;

      // add directions[i], lower_bounds[i], and upper_bounds[i]
      // in new_directions, new_lower_bounds, and new_upper_bounds,
      // respectively
      new_directions.push_back(directions[i]);
      new_lower_bounds.push_back(lower_bounds[i]);
      new_upper_bounds.push_back(upper_bounds[i]);

      // mark in new_pos directions in direction vector that are
      // linearly dependent to directions[i]
      mark_linearly_dep_dirs(new_pos, directions, i);
    } else {
      // directions[i] has been already inserted in new_directions
      // in position new_pos[i]
      using namespace LinearAlgebra;

      // get the index of directions[i] in new_directions
      const size_t new_idx = new_pos[i];

      const T coeff = new_directions[new_idx] / directions[i];

      // if necessary, increase the lower bound
      new_lower_bounds[new_idx]
          = std::max(new_lower_bounds[new_idx], lower_bounds[i] * coeff);

      // if necessary, decrease the upper bound
      new_upper_bounds[new_idx]
          = std::min(new_upper_bounds[new_idx], upper_bounds[i] * coeff);
    }
  }

  // replace old directions, lower_bounds, and upper_bounds
  std::swap(directions, new_directions);
  std::swap(lower_bounds, new_lower_bounds);
  std::swap(upper_bounds, new_upper_bounds);

  return new_pos;
}

/**
 * @brief Reorganize directions according to a set of new positions
 *
 * This function reorganizes the directions according to a new set
 * of positions and removes all the directions that do not have a
 * new position.
 *
 * @tparam T is the scalar type of the directions
 * @param directions is the vector of the directions to be reorganized
 * @param new_positions is the map of the new positions
 */
template<typename T>
void resort_directions(std::vector<LinearAlgebra::Vector<T>> &directions,
                       LinearAlgebra::Vector<T> &lower_bounds,
                       LinearAlgebra::Vector<T> &upper_bounds,
                       const std::map<size_t, size_t> &new_positions)
{
  std::vector<LinearAlgebra::Vector<T>> new_directions(new_positions.size());
  LinearAlgebra::Vector<T> new_lower_bounds(new_positions.size()),
      new_upper_bounds(new_positions.size());

  for (const auto &new_pos: new_positions) {
    std::swap(directions[new_pos.first], new_directions[new_pos.second]);
    std::swap(lower_bounds[new_pos.first], new_lower_bounds[new_pos.second]);
    std::swap(upper_bounds[new_pos.first], new_upper_bounds[new_pos.second]);
  }

  std::swap(directions, new_directions);
  std::swap(lower_bounds, new_lower_bounds);
  std::swap(upper_bounds, new_upper_bounds);
}

/**
 * @brief Update the direction positions in a set of templates
 *
 * This function returns a set of templates which are built by updating the
 * indices of those in the input set indices according to the vector of the new
 * direction positions. Moreover, the new templates are brought to the
 * canonical form.
 *
 * @param templates is the set of the templates whose direction position
 *               must be updated
 * @param new_positions is a vector mapping each old direction index in the new
 *               direction index
 * @return the updated set of templates
 */
std::set<std::vector<unsigned int>> update_templates_directions(
    const std::set<std::vector<unsigned int>> &templates,
    const std::vector<size_t> &new_positions)
{
  std::set<std::vector<unsigned int>> new_templates;

  for (const auto &bundle_template: templates) {
    std::vector<unsigned int> new_template(bundle_template);
    for (auto &dir_id: new_template) {
      dir_id = new_positions[dir_id];
    }
    canonize_template(new_template);

    new_templates.insert(std::move(new_template));
  }

  return new_templates;
}

/**
 * @brief Update the direction positions in a set of templates
 *
 * This function returns a set of templates which are built by updating the
 * indices of those in the input set indices according to the map of the new
 * direction positions. Moreover, the new templates are brought to the
 * canonical form.
 *
 * @param templates is the set of the templates whose direction position
 *               must be updated
 * @param new_positions is the map of new direction positions
 * @return the updated set of templates
 */
std::set<std::vector<unsigned int>> update_templates_directions(
    const std::set<std::vector<unsigned int>> &templates,
    const std::map<size_t, size_t> &new_positions)
{
  std::set<std::vector<unsigned int>> new_templates;

  for (const auto &bundle_template: templates) {
    std::vector<unsigned int> new_template(bundle_template);
    for (auto &dir_id: new_template) {
      dir_id = new_positions.at(dir_id);
    }
    canonize_template(new_template);

    new_templates.insert(std::move(new_template));
  }

  return new_templates;
}

template<typename T>
void add_missing_templates(
    std::set<std::vector<unsigned int>> &templates,
    const std::vector<LinearAlgebra::Vector<T>> &directions,
    std::set<size_t> missing_dirs)
{
  using namespace LinearAlgebra;
  using namespace LinearAlgebra::Dense;

  while (!missing_dirs.empty()) {
    // fill T with all the bundle directions starting
    // from those not included in a template
    Matrix<double> D = directions;
    Vector<size_t> row_pos(D.size());
    std::iota(std::begin(row_pos), std::end(row_pos), 0);

    size_t j = 0;
    for (const auto &i: missing_dirs) {
      std::swap(D[i], D[j]);
      std::swap(row_pos[i], row_pos[j]);
      ++j;
    }

    auto indep_rows = find_first_independent_rows(D);
    Vector<unsigned int> new_template;

    for (const auto &value: indep_rows) {
      new_template.push_back(row_pos[value]);
    }

    canonize_template(new_template);

    templates.insert(new_template);
    for (const auto &idx: new_template) {
      missing_dirs.erase(idx);
    }
  }
}

template<typename T>
void add_missing_templates(
    std::set<BundleTemplate> &templates,
    const std::vector<LinearAlgebra::Vector<T>> &directions,
    std::set<size_t> missing_dirs,
    const Bundle::template_less_fate_type template_less_fate
    = Bundle::STATIC_TEMPLATES)
{
  std::set<std::vector<unsigned int>> raw_templates;

  add_missing_templates(raw_templates, directions, missing_dirs);

  for (auto &raw_template: raw_templates) {
    templates.emplace(std::move(raw_template),
                      template_less_fate == Bundle::DYNAMIC_TEMPLATES);
  }
}

/**
 * @brief Test whether a template is a complete basis for the space
 *
 * @tparam T is the scalar type for the directions
 * @param directions is the vector of directions
 * @param bundle_template is the template to be tested
 * @return true if and only if the set of the directions associated to
 *       bundle_template is a complete basis for the space
 */
template<typename T>
bool template_is_a_complete_basis(
    const std::vector<LinearAlgebra::Vector<T>> &directions,
    const std::vector<unsigned int> &bundle_template)
{
  using namespace LinearAlgebra;

  Dense::Matrix<T> A;
  for (const auto &dir_idx: bundle_template) {
    if (dir_idx >= directions.size()) {
      SAPO_ERROR("templates must contains as values "
                 "indices of the directions vectors",
                 std::domain_error);
    }

    A.emplace_back(directions[dir_idx]);
  }

  return Dense::rank(A) == A.size();
}

template<typename T>
void validate_templates(
    const std::vector<LinearAlgebra::Vector<T>> &directions,
    const std::set<std::vector<unsigned int>> &templates)
{
  if (directions.size() == 0) {
    SAPO_ERROR("directions must be non empty", std::domain_error);
  }

  const size_t dim = directions[0].size();

  for (const auto &bundle_template: templates) {
    if (bundle_template.size() != dim) {
      SAPO_ERROR("templates must have as many columns "
                 "as directions",
                 std::domain_error);
    }

    if (!template_is_a_complete_basis(directions, bundle_template)) {
      std::ostringstream oss;

      oss << bundle_template << " does not contain linearly "
          << "independent directions";
      SAPO_ERROR(oss.str(), std::domain_error);
    }
  }
}

template<typename T>
void validate_directions(
    const std::vector<LinearAlgebra::Vector<T>> &directions,
    const LinearAlgebra::Vector<T> &lower_bounds,
    const LinearAlgebra::Vector<T> &upper_bounds)
{
  if (directions.size() == 0) {
    SAPO_ERROR("direction vector must be non empty", std::domain_error);
  }

  const size_t dim = directions[0].size();
  for (const auto &dir: directions) {
    if (dir.size() != dim) {
      SAPO_ERROR("all the directions must have the same dimension",
                 std::domain_error);
    }

    if (LinearAlgebra::norm_infinity(dir) == 0) {
      SAPO_ERROR("all the directions must be non-null", std::domain_error);
    }
  }

  if (directions.size() != upper_bounds.size()) {
    SAPO_ERROR("directions and upper_bounds must have the same size",
               std::domain_error);
  }
  if (directions.size() != lower_bounds.size()) {
    SAPO_ERROR("directions and lower_bounds must have the same size",
               std::domain_error);
  }
}

void duplicate_dynamic_directions(
    std::vector<LinearAlgebra::Vector<double>> &directions,
    LinearAlgebra::Vector<double> &lower_bounds,
    LinearAlgebra::Vector<double> &upper_bounds,
    const std::set<std::vector<unsigned int>> &static_templates,
    std::set<std::vector<unsigned int>> &dynamic_templates)
{
  // mark as to-be-duplicated all the static directions
  std::vector<bool> require_duplicate(directions.size(), false);
  for (const auto &bundle_template: static_templates) {
    for (const auto &idx: bundle_template) {
      require_duplicate[idx] = true;
    }
  }

  // build a new dynamic template set
  std::set<std::vector<unsigned int>> new_dynamic_templates;

  // for each template in the dynamic template set
  for (auto bundle_template: dynamic_templates) {
    bool changed = false;
    for (auto &idx: bundle_template) {

      // if one of template directions has been already mentioned
      if (require_duplicate[idx]) {

        // the template must be changed
        changed = true;

        // copy the mentioned direction as a new direction
        directions.emplace_back(directions[idx]);
        lower_bounds.push_back(lower_bounds[idx]);
        upper_bounds.push_back(upper_bounds[idx]);

        // change the direction id in the template
        idx = directions.size() - 1;
      } else { // if not yet mentioned, mark as to-be-duplicated
        require_duplicate[idx] = true;
      }
    }

    // if the bundle changed
    if (changed) {
      canonize_template(bundle_template);
    }
    new_dynamic_templates.insert(bundle_template);
  }

  std::swap(dynamic_templates, new_dynamic_templates);
}

Bundle::Bundle(const std::vector<LinearAlgebra::Vector<double>> &directions,
               const LinearAlgebra::Vector<double> &lower_bounds,
               const LinearAlgebra::Vector<double> &upper_bounds,
               const std::set<BundleTemplate> &templates):
    _directions(directions),
    _lower_bounds(lower_bounds), _upper_bounds(upper_bounds),
    _templates(templates)
{
}

Bundle::Bundle(std::vector<LinearAlgebra::Vector<double>> &&directions,
               LinearAlgebra::Vector<double> &&lower_bounds,
               LinearAlgebra::Vector<double> &&upper_bounds,
               const std::set<BundleTemplate> &templates):
    _directions(std::move(directions)),
    _lower_bounds(std::move(lower_bounds)),
    _upper_bounds(std::move(upper_bounds)), _templates(templates)
{
}

Bundle::Bundle(const std::vector<LinearAlgebra::Vector<double>> &directions,
               const LinearAlgebra::Vector<double> &lower_bounds,
               const LinearAlgebra::Vector<double> &upper_bounds,
               std::set<std::vector<unsigned int>> static_templates,
               std::set<std::vector<unsigned int>> dynamic_templates,
               const template_less_fate_type template_less_fate):
    _directions(directions),
    _lower_bounds(lower_bounds), _upper_bounds(upper_bounds), _templates()
{
  validate_directions(_directions, _lower_bounds, _upper_bounds);
  validate_templates(_directions, static_templates);
  validate_templates(_directions, dynamic_templates);

  // filter duplicated direction, update templates, and bring them in canonical
  // form
  auto new_pos = filter_duplicated_directions(_directions, _lower_bounds,
                                              _upper_bounds);

  static_templates = update_templates_directions(static_templates, new_pos);
  dynamic_templates = update_templates_directions(dynamic_templates, new_pos);

  auto used_dirs = collect_directions_in_templates(static_templates);
  collect_directions_in_templates(used_dirs, dynamic_templates);

  if (template_less_fate != REMOVE_DIRECTION) {
    // add missing directions in templates
    std::set<size_t> missing_dirs;
    for (size_t idx = 0; idx < _directions.size(); ++idx) {
      if (used_dirs.count(idx) == 0) {
        missing_dirs.insert(idx);
      }
    }
    if (template_less_fate == STATIC_TEMPLATES) {
      add_missing_templates(static_templates, _directions, missing_dirs);
    } else {
      add_missing_templates(dynamic_templates, _directions, missing_dirs);
    }
  } else {
    if (static_templates.size() == 0 && dynamic_templates.size() == 0) {
      SAPO_ERROR("template vector must be non empty", std::domain_error);
    }

    // if the number of directions appearing in the templates is smaller
    // than the number of directions
    if (used_dirs.size() != directions.size()) {

      // find a new position for the used directions
      auto new_pos = find_new_position_for(used_dirs);

      resort_directions(_directions, _lower_bounds, _upper_bounds, new_pos);
      static_templates
          = update_templates_directions(static_templates, new_pos);
      dynamic_templates
          = update_templates_directions(dynamic_templates, new_pos);
    }
  }
  duplicate_dynamic_directions(_directions, _lower_bounds, _upper_bounds,
                               static_templates, dynamic_templates);

  _templates.insert(std::begin(static_templates), std::end(static_templates));

  for (auto &dynamic_template: dynamic_templates) {
    _templates.emplace(std::move(dynamic_template), true);
  }
}

Bundle::Bundle(const std::vector<LinearAlgebra::Vector<double>> &directions,
               const LinearAlgebra::Vector<double> &lower_bounds,
               const LinearAlgebra::Vector<double> &upper_bounds,
               const std::set<std::vector<unsigned int>> &static_templates,
               const template_less_fate_type template_less_fate):
    Bundle(directions, lower_bounds, upper_bounds, static_templates, {},
           template_less_fate)
{
}

Bundle::Bundle(std::vector<LinearAlgebra::Vector<double>> &&directions,
               LinearAlgebra::Vector<double> &&lower_bounds,
               LinearAlgebra::Vector<double> &&upper_bounds,
               const std::set<std::vector<unsigned int>> &static_templates,
               const template_less_fate_type template_less_fate):
    Bundle(std::move(directions), std::move(lower_bounds),
           std::move(upper_bounds), static_templates, {}, template_less_fate)
{
}

Bundle::Bundle(const std::vector<LinearAlgebra::Vector<double>> &directions,
               const LinearAlgebra::Vector<double> &lower_bounds,
               const LinearAlgebra::Vector<double> &upper_bounds):
    Bundle(directions, lower_bounds, upper_bounds, {}, {}, STATIC_TEMPLATES)
{
}

Bundle::Bundle(std::vector<LinearAlgebra::Vector<double>> &&directions,
               LinearAlgebra::Vector<double> &&lower_bounds,
               LinearAlgebra::Vector<double> &&upper_bounds):
    Bundle(std::move(directions), std::move(lower_bounds),
           std::move(upper_bounds), {}, {}, STATIC_TEMPLATES)
{
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

  std::set<size_t> missing_dirs;
  for (size_t i = 0; i < _directions.size(); ++i) {
    missing_dirs.insert(i);
  }

  add_missing_templates(_templates, _directions, missing_dirs);
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
  this->_templates = std::set<BundleTemplate>();

  this->_directions.reserve(orig._directions.size());

  for (const auto &dir: orig._directions) {
    this->_directions.emplace_back(dir);
  }

  for (const auto &bundle_template: orig._templates) {
    this->_templates.emplace(bundle_template);
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

Parallelotope
Bundle::get_parallelotope(const BundleTemplate &bundle_template) const
{
  using namespace std;

  vector<double> lbound, ubound;
  vector<LinearAlgebra::Vector<double>> Lambda;

  auto it = std::begin(bundle_template.get_direction_indices());
  // upper facets
  for (unsigned int j = 0; j < this->dim(); j++) {
    const unsigned int idx = *it;
    if (idx >= this->_directions.size()) {
      SAPO_ERROR("the parameter is not a template for the bundle",
                 std::domain_error);
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
    SAPO_ERROR("the two bundles must have the same dimension",
               std::domain_error);
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
std::list<Bundle> &Bundle::split_bundle(
    std::list<Bundle> &res, LinearAlgebra::Vector<double> &lower_bounds,
    LinearAlgebra::Vector<double> &upper_bounds, const size_t idx,
    const double &max_magnitude, const double &split_ratio) const
{
  if (idx == this->size()) {
    Bundle new_bundle(this->_directions, this->_lower_bounds,
                      this->_upper_bounds, this->_templates);
    res.push_back(std::move(new_bundle));

    return res;
  }

  if (std::abs(get_upper_bound(idx) - get_lower_bound(idx)) > max_magnitude) {
    double lower_bound = get_lower_bound(idx);

    do {
      const double upper_bound = std::min(
          lower_bound + split_ratio * max_magnitude, get_upper_bound(idx));

      lower_bounds[idx] = lower_bound;
      upper_bounds[idx] = upper_bound;
      this->split_bundle(res, lower_bounds, upper_bounds, idx + 1,
                         max_magnitude, split_ratio);

      lower_bound = upper_bound;
    } while (get_upper_bound(idx) != lower_bound);
  } else {
    lower_bounds[idx] = get_lower_bound(idx);
    upper_bounds[idx] = get_upper_bound(idx);
    this->split_bundle(res, lower_bounds, upper_bounds, idx + 1, max_magnitude,
                       split_ratio);
  }
  return res;
}

std::list<Bundle> Bundle::split(const double max_magnitude,
                                const double split_ratio) const
{
  std::list<Bundle> split_list;

  LinearAlgebra::Vector<double> upper_bounds(this->size());
  LinearAlgebra::Vector<double> lower_bounds(this->size());

  this->split_bundle(split_list, lower_bounds, upper_bounds, 0, max_magnitude,
                     split_ratio);

  return split_list;
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
void add_mapped_templates(
    std::set<BundleTemplate> &dest_templates,
    const std::set<BundleTemplate> &source_templates,
    const std::vector<unsigned int> &new_direction_indices)
{
  // check whether some of the templates must be copied
  for (const BundleTemplate &temp: source_templates) {
    // if this is the case, copy and update the template indices
    std::vector<unsigned int> t_copy(temp.get_direction_indices());
    for (unsigned int j = 0; j < t_copy.size(); ++j) {
      t_copy[j] = new_direction_indices[t_copy[j]];
    }

    // add the new template to the intersected bundle
    dest_templates.emplace(std::move(t_copy), temp.have_dynamic_directions());
  }
}

Bundle &Bundle::intersect_with(const Bundle &A)
{
  using namespace LinearAlgebra;

  if (dim() != A.dim()) {
    SAPO_ERROR("the two bundles differ in dimensions", std::domain_error);
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
      double scaled_A_upper_bound, scaled_A_lower_bound;

      if (dep_coeff > 0) {
        scaled_A_lower_bound = A._lower_bounds[i] * dep_coeff;
        scaled_A_upper_bound = A._upper_bounds[i] * dep_coeff;
      } else {
        scaled_A_lower_bound = A._upper_bounds[i] * dep_coeff;
        scaled_A_upper_bound = A._lower_bounds[i] * dep_coeff;
      }

      // increase the lower bound if necessary
      if (scaled_A_lower_bound > this->_lower_bounds[new_i]) {
        this->_lower_bounds[new_i] = scaled_A_lower_bound;
      }

      // decrease the upper bound if necessary
      if (scaled_A_upper_bound < this->_upper_bounds[new_i]) {
        this->_upper_bounds[new_i] = scaled_A_upper_bound;
      }
    }
  }

  add_mapped_templates(_templates, A.templates(), new_ids);

  return *this;
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
    SAPO_ERROR("the bundle and the linear system differ in "
               "dimensions",
               std::domain_error);
  }

  std::set<size_t> outside_templates;
  const Matrix<double> &A = ls.A();

  std::vector<unsigned int> new_ids(A.size());
  // for each row in the linear system
  for (size_t i = 0; i < A.size(); ++i) {
    const LinearAlgebra::Vector<double> &A_row = A[i];

    new_ids[i] = get_a_linearly_dependent_in(A_row, this->_directions);

    // if the direction is not present in this object
    if (new_ids[i] == this->size()) {

      // add the direction and the corresponding boundaries
      outside_templates.insert(new_ids[i]);
      this->_directions.push_back(A_row);
      this->_lower_bounds.push_back(-std::numeric_limits<double>::infinity());
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
  add_missing_templates(_templates, _directions, outside_templates);

  return this->canonize();
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
    SAPO_ERROR("the two bundles differ in dimensions", std::domain_error);
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

  add_mapped_templates(res._templates, b2.templates(), new_ids);

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
