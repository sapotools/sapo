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
    _dir_indices(orig._dir_indices), _adaptive_indices(orig._adaptive_indices)
{
}

BundleTemplate::BundleTemplate(BundleTemplate &&orig):
    _dir_indices(std::move(orig._dir_indices)),
    _adaptive_indices(std::move(orig._adaptive_indices))
{
}

BundleTemplate::BundleTemplate(const std::vector<unsigned int> &dir_indices,
                               const std::set<size_t> &adaptive_indices):
    _dir_indices(dir_indices.size())
{
  std::vector<size_t> position(dir_indices.size());
  std::iota(position.begin(), position.end(), 0);

  std::sort(std::begin(position), std::end(position),
            [&dir_indices](const size_t &a, const size_t &b) {
              return dir_indices[a] < dir_indices[b];
            });

  for (size_t i = 0; i < dir_indices.size(); ++i) {
    _dir_indices[position[i]] = dir_indices[i];
  }

  for (size_t i = 1; i < _dir_indices.size(); ++i) {
    if (_dir_indices[i] == _dir_indices[i - 1]) {
      SAPO_ERROR("Template contains twice the same "
                 "direction index",
                 std::domain_error);
    }
  }

  for (size_t i = 0; i < _dir_indices.size(); ++i) {
    if (adaptive_indices.count(_dir_indices[i]) > 0) {
      _adaptive_indices.insert(i);
    }
  }
}

std::ostream &operator<<(std::ostream &os,
                         const BundleTemplate &bundle_template)
{
  os << "{direction_indices: [";
  std::string sep = "";
  for (const auto &index: bundle_template.direction_indices()) {
    os << sep << index;
    sep = ",";
  }
  os << "], adaptive_indices: {";
  sep = "";
  for (const auto &index: bundle_template.adaptive_indices()) {
    os << sep << index;
    sep = ",";
  }
  os << "}}";

  return os;
}

bool std::less<BundleTemplate>::operator()(const BundleTemplate &a,
                                           const BundleTemplate &b) const
{
  std::less<std::vector<unsigned int>> cmp;

  return cmp(a.direction_indices(), b.direction_indices());
}

Bundle::Bundle() {}

Bundle::Bundle(const Bundle &orig):
    _directions(orig._directions),
    _adaptive_directions(orig._adaptive_directions),
    _lower_bounds(orig._lower_bounds), _upper_bounds(orig._upper_bounds),
    _templates(orig._templates)
{
}

Bundle::Bundle(Bundle &&orig)
{
  swap(*this, orig);
}

void swap(Bundle &A, Bundle &B)
{
  std::swap(A._directions, B._directions);
  std::swap(A._adaptive_directions, B._adaptive_directions);
  std::swap(A._upper_bounds, B._upper_bounds);
  std::swap(A._lower_bounds, B._lower_bounds);
  std::swap(A._templates, B._templates);
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
 * @param templates is the set of templates whose directions is aimed for
 * @return the set of directions mentioned in the templates
 */
std::set<size_t> collect_directions_in_templates(
    const std::set<std::vector<unsigned int>> &templates)
{
  std::set<size_t> template_directions;
  for (const auto &bundle_template: templates) {
    template_directions.insert(std::begin(bundle_template),
                               std::end(bundle_template));
  }

  return template_directions;
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
    const std::set<size_t> &adaptive_directions, const size_t &dir_idx)
{
  using namespace LinearAlgebra;

  bool dir_is_dyn = (adaptive_directions.count(dir_idx) > 0);

  const size_t new_idx = new_pos[dir_idx];

  const auto &dir = directions[dir_idx];
  for (size_t j = dir_idx + 1; j < directions.size(); ++j) {
    if (new_pos[j] == j) {
      bool j_is_dyn = (adaptive_directions.count(j) > 0);
      if (j_is_dyn == dir_is_dyn
          && are_linearly_dependent(dir, directions[j])) {
        new_pos[j] = new_idx;
      }
    }
  }
}

/**
 * @brief Filter duplicate directions
 *
 * @tparam T is the scalar type of the directions
 * @param directions is the direction vector
 * @param adaptive_directions is the set of adaptive direction indices
 * @param lower_bounds is the vector of the direction lower bounds
 * @param upper_bounds is the vector of the direction upper bounds
 * @return a vector mapping each index of the input direction vector in the
 *       corresponding index in the final direction vector
 */
template<typename T>
std::vector<size_t>
filter_duplicated_directions(std::vector<LinearAlgebra::Vector<T>> &directions,
                             std::set<size_t> &adaptive_directions,
                             LinearAlgebra::Vector<T> &lower_bounds,
                             LinearAlgebra::Vector<T> &upper_bounds)
{
  std::vector<LinearAlgebra::Vector<T>> new_directions;
  std::set<size_t> new_adaptive_dirs;
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

      // add directions[i], lower_bounds[i], and
      // upper_bounds[i] in new_directions, new_lower_bounds,
      // and new_upper_bounds, respectively
      new_directions.push_back(directions[i]);
      new_lower_bounds.push_back(lower_bounds[i]);
      new_upper_bounds.push_back(upper_bounds[i]);

      // if the i-th direction is dynamic add among the
      // adaptive directions of the new direction vector
      if (adaptive_directions.count(i) > 0) {
        new_adaptive_dirs.insert(new_idx);
      }

      // mark in new_pos directions all direction vector that are
      // linearly dependent to directions[i] and are among dynamic
      // directions as direction[i] is/is not.
      mark_linearly_dep_dirs(new_pos, directions, adaptive_directions, i);
    } else {
      // directions[i] has been already inserted in new_directions
      // in position new_pos[i]
      using namespace LinearAlgebra;

      // get the index of directions[i] in new_directions
      const size_t new_idx = new_pos[i];

      const T coeff = new_directions[new_idx] / directions[i];

      // if coeff < 0, then coeff*lower_bounds[i] and 
      // coeff*upper_bounds[i] are the new potential 
      // upper and lower bounds, respectively 
      if (coeff < 0) {
        std::swap(lower_bounds[i], upper_bounds[i]);
      }

      // if necessary, increase the lower bound
      new_lower_bounds[new_idx]
          = std::max(new_lower_bounds[new_idx], lower_bounds[i] * coeff);

      // if necessary, decrease the upper bound
      new_upper_bounds[new_idx]
          = std::min(new_upper_bounds[new_idx], upper_bounds[i] * coeff);
    }
  }

  // replace old directions, adaptive_directions, lower_bounds, and
  // upper_bounds
  std::swap(directions, new_directions);
  std::swap(adaptive_directions, new_adaptive_dirs);
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
 * @param directions is the vector of the directions
 * @param adaptive_directions is the set of adaptive direction indices
 * @param lower_bounds is the vector of the direction lower bounds
 * @param upper_bounds is the vector of the direction upper bounds
 * @param new_positions is the map of the new positions
 */
template<typename T>
void resort_directions(std::vector<LinearAlgebra::Vector<T>> &directions,
                       std::set<size_t> &adaptive_directions,
                       LinearAlgebra::Vector<T> &lower_bounds,
                       LinearAlgebra::Vector<T> &upper_bounds,
                       const std::map<size_t, size_t> &new_positions)
{
  std::vector<LinearAlgebra::Vector<T>> new_directions(new_positions.size());
  std::set<size_t> new_adaptive_dirs;
  LinearAlgebra::Vector<T> new_lower_bounds(new_positions.size()),
      new_upper_bounds(new_positions.size());

  for (const auto &new_pos: new_positions) {
    std::swap(directions[new_pos.first], new_directions[new_pos.second]);
    std::swap(lower_bounds[new_pos.first], new_lower_bounds[new_pos.second]);
    std::swap(upper_bounds[new_pos.first], new_upper_bounds[new_pos.second]);
  }

  for (const auto &adaptive_direction: adaptive_directions) {
    new_adaptive_dirs.insert(new_positions.at(adaptive_direction));
  }

  std::swap(directions, new_directions);
  std::swap(adaptive_directions, new_adaptive_dirs);
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
    const std::set<size_t> &adaptive_directions, std::set<size_t> missing_dirs)
{
  std::set<std::vector<unsigned int>> raw_templates;

  add_missing_templates(raw_templates, directions, missing_dirs);

  for (auto &raw_template: raw_templates) {
    templates.emplace(raw_template, adaptive_directions);
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

  return rank(A) == A.size();
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

void duplicate_adaptive_directions(
    std::vector<LinearAlgebra::Vector<double>> &directions,
    std::set<size_t> &adaptive_directions,
    LinearAlgebra::Vector<double> &lower_bounds,
    LinearAlgebra::Vector<double> &upper_bounds,
    std::set<std::vector<unsigned int>> &templates)
{
  std::vector<bool> require_duplicate(directions.size(), false);

  // build a new dynamic template set
  std::set<std::vector<unsigned int>> new_templates;

  // for each template in the dynamic template set
  for (auto bundle_template: templates) {
    bool changed = false;
    for (auto &idx: bundle_template) {

      // if one of template directions has been already mentioned
      if (adaptive_directions.count(idx) > 0 && require_duplicate[idx]) {

        // the template must be changed
        changed = true;

        // add the index of the direction that we are going to
        // insert at the end of the direction vector among
        // the adaptive direction indices
        adaptive_directions.insert(directions.size());

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
    new_templates.insert(bundle_template);
  }

  std::swap(templates, new_templates);
}

Bundle::Bundle(const std::set<size_t> &adaptive_directions,
               const std::vector<LinearAlgebra::Vector<double>> &directions,
               const LinearAlgebra::Vector<double> &lower_bounds,
               const LinearAlgebra::Vector<double> &upper_bounds,
               const std::set<BundleTemplate> &templates):
    _directions(directions),
    _adaptive_directions(adaptive_directions), _lower_bounds(lower_bounds),
    _upper_bounds(upper_bounds), _templates(templates)
{
}

Bundle::Bundle(const std::set<size_t> &adaptive_directions,
               std::vector<LinearAlgebra::Vector<double>> &&directions,
               LinearAlgebra::Vector<double> &&lower_bounds,
               LinearAlgebra::Vector<double> &&upper_bounds,
               const std::set<BundleTemplate> &templates):
    _directions(std::move(directions)),
    _adaptive_directions(adaptive_directions),
    _lower_bounds(std::move(lower_bounds)),
    _upper_bounds(std::move(upper_bounds)), _templates(templates)
{
}

Bundle::Bundle(const std::vector<LinearAlgebra::Vector<double>> &directions,
               const LinearAlgebra::Vector<double> &lower_bounds,
               const LinearAlgebra::Vector<double> &upper_bounds,
               std::set<std::vector<unsigned int>> templates,
               const std::set<size_t> &adaptive_directions,
               const bool remove_unused_directions):
    _directions(directions),
    _adaptive_directions(adaptive_directions), _lower_bounds(lower_bounds),
    _upper_bounds(upper_bounds), _templates()
{
  validate_directions(_directions, _lower_bounds, _upper_bounds);
  validate_templates(_directions, templates);

  for (const auto &dir_index: adaptive_directions) {
    if (dir_index >= _directions.size()) {
      SAPO_ERROR("the adaptive directions must be valid "
                 "indices for the direction vector",
                 std::domain_error);
    }
  }

  // filter duplicated direction, update templates, and bring them in canonical
  // form
  auto new_pos = filter_duplicated_directions(
      _directions, _adaptive_directions, _lower_bounds, _upper_bounds);

  templates = update_templates_directions(templates, new_pos);

  auto used_dirs = collect_directions_in_templates(templates);

  if (!remove_unused_directions) {
    // add missing directions in templates
    std::set<size_t> missing_dirs;
    for (size_t idx = 0; idx < _directions.size(); ++idx) {
      if (used_dirs.count(idx) == 0) {
        missing_dirs.insert(idx);
      }
    }
    add_missing_templates(templates, _directions, missing_dirs);
  } else {
    if (templates.size() == 0) {
      SAPO_ERROR("template vector must be non empty", std::domain_error);
    }

    // if the number of directions appearing in the templates is smaller
    // than the number of directions
    if (used_dirs.size() != directions.size()) {

      // find a new position for the used directions
      auto new_pos = find_new_position_for(used_dirs);

      resort_directions(_directions, _adaptive_directions, _lower_bounds,
                        _upper_bounds, new_pos);
      templates = update_templates_directions(templates, new_pos);
    }
  }
  duplicate_adaptive_directions(_directions, _adaptive_directions,
                                _lower_bounds, _upper_bounds, templates);

  for (auto &bundle_template: templates) {
    _templates.emplace(bundle_template, _adaptive_directions);
  }
}

Bundle::Bundle(const std::vector<LinearAlgebra::Vector<double>> &directions,
               const LinearAlgebra::Vector<double> &lower_bounds,
               const LinearAlgebra::Vector<double> &upper_bounds,
               const std::set<std::vector<unsigned int>> &templates,
               const bool remove_unused_directions):
    Bundle(directions, lower_bounds, upper_bounds, templates, {},
           remove_unused_directions)
{
}

Bundle::Bundle(const std::vector<LinearAlgebra::Vector<double>> &directions,
               const LinearAlgebra::Vector<double> &lower_bounds,
               const LinearAlgebra::Vector<double> &upper_bounds,
               const bool remove_unused_directions):
    Bundle(directions, lower_bounds, upper_bounds, {}, {},
           remove_unused_directions)
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
      _upper_bounds.push_back(P.maximize(P.A(i)).objective_value());
      _lower_bounds.push_back(P.minimize(P.A(i)).objective_value());
    }
  }

  std::set<size_t> missing_dirs;
  for (size_t i = 0; i < _directions.size(); ++i) {
    missing_dirs.insert(i);
  }

  add_missing_templates(_templates, _directions, {}, missing_dirs);
}

/**
 * @brief Copy operator
 *
 * @param orig is the original object to be copied
 * @return a reference to the updated object
 */
Bundle &Bundle::operator=(const Bundle &orig)
{
  this->_directions = orig._directions;
  this->_adaptive_directions = orig._adaptive_directions;
  this->_templates = orig._templates;

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

  auto it = std::begin(bundle_template.direction_indices());
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
    _lower_bounds[i] = bund.minimize(this->_directions[i]).objective_value();
    _upper_bounds[i] = bund.maximize(this->_directions[i]).objective_value();
  }
  return *this;
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

  using namespace LinearAlgebra;

  Dense::Matrix<double> ls_A;
  Vector<double> ls_b;

  auto add_bundle_constraints = [&ls_A, &ls_b](const Bundle &bundle) {
    for (size_t i = 0; i < bundle.size(); ++i) {
      ls_A.push_back(bundle.get_direction(i));
      ls_b.push_back(bundle.get_upper_bound(i));
      ls_A.push_back(-bundle.get_direction(i));
      ls_b.push_back(-bundle.get_lower_bound(i));
    }
  };

  add_bundle_constraints(A);
  add_bundle_constraints(B);

  SimplexMethodOptimizer optimizer;

  auto result = optimizer(ls_A, ls_b, Vector<double>(A.dim()));

  return result.feasible_set_is_empty();
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
    if (P_this.minimize(bundle.get_direction(dir_idx)).objective_value()
        < bundle.get_lower_bound(dir_idx)) {
      return false;
    }

    // if the maximum of this object on that direction is greater than
    // the bundle maximum, this object is not a subset of the bundle
    if (P_this.maximize(bundle.get_direction(dir_idx)).objective_value()
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
    if (P_this.maximize(ls.A(dir_idx)).objective_value() > ls.b(dir_idx)) {
      return false;
    }
  }

  // if none of the previous conditions hold, the bundle is a
  // subset for this object
  return true;
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
    Bundle new_bundle(this->_adaptive_directions, this->_directions,
                      this->_lower_bounds, this->_upper_bounds,
                      this->_templates);
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
 * vector, the templates involving non-shared must be added.
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
    std::vector<unsigned int> t_copy(temp.direction_indices());
    for (unsigned int j = 0; j < t_copy.size(); ++j) {
      t_copy[j] = new_direction_indices[t_copy[j]];
    }

    // add the new template to the intersected bundle
    dest_templates.emplace(t_copy, std::set<size_t>{});
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
  add_missing_templates(_templates, _directions, _adaptive_directions,
                        outside_templates);

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
    res._lower_bounds[i] = std::min(p2.minimize(res_dir).objective_value(),
                                    res.get_lower_bound(i));
    res._upper_bounds[i] = std::max(p2.maximize(res_dir).objective_value(),
                                    res.get_upper_bound(i));
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
      double lower_bound = std::min(p1.minimize(b2_dir).objective_value(),
                                    b2.get_lower_bound(i));
      double upper_bound = std::max(p1.maximize(b2_dir).objective_value(),
                                    b2.get_upper_bound(i));

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
    auto new_bound = p1.maximize(b2.get_direction(i)).objective_value();
    if (new_bound > b2.get_upper_bound(i)) {
      Bundle new_b1 = b1;
      new_b1._directions.push_back(b2.get_direction(i));
      new_b1._lower_bounds.push_back(b2.get_upper_bound(i));
      new_b1._upper_bounds.push_back(new_bound);

      su.add(std::move(new_b1.canonize()));
    }

    new_bound = p1.minimize(b2.get_direction(i)).objective_value();
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
