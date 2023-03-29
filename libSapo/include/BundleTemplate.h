/**
 * @file BundleTemplate.h
 * @author Alberto Casagrande <acasagrande@units.it>
 * @brief Represent bundle templates
 * @version 0.1
 * @date 2023-03-29
 *
 * @copyright Copyright (c) 2023
 */

#ifndef _BUNDLE_TEMPLATE_H_
#define _BUNDLE_TEMPLATE_H_

#include <set>
#include <vector>

#include "LinearAlgebra.h"

template<typename T, typename APPROX_TYPE>
class Bundle;

/**
 * @brief Bundle template
 *
 * This class represents bundle templates. A template of a $n$-dimensional
 * bundle is a vector of the indices of $n$ linearly independent directions of
 * the bundle itself. A bundle consists in the intersection of parallelotopes.
 * A template symbolically represents one of these parallelotopes.
 */
class BundleTemplate
{
  std::vector<unsigned int>
      _dir_indices; //!< indices of the template directions
  std::set<size_t>
      _adaptive_indices; //!< the template adaptive direction indices

  /**
   * @private
   * @brief Canonize an index vector representing a template
   *
   * A template is in canonical form if it is sorted. This function
   * sorts the direction indices of a template and brings it in
   * canonical form.
   *
   * @param bundle_template is the template to be brought in canonical form
   * @return a reference to the updated template
   */
  static std::vector<unsigned int> &
  canonize(std::vector<unsigned int> &bundle_template);

  /**
   * @private
   * @brief Update the direction positions in a set of templates
   *
   * This method returns a set of vectors representing templates. The
   * vectors themselves are built by updating the indices of those in
   * the input set according to the vector of the new direction
   * positions. The returned templates are also brought to the
   * canonical form.
   *
   * @param templates is the set of the templates whose direction position
   *               must be updated
   * @param new_positions is a vector mapping each old direction index in the
   * new direction index
   * @return the updated set of templates
   */
  static std::set<std::vector<unsigned int>>
  update_directions(const std::set<std::vector<unsigned int>> &templates,
                    const std::vector<size_t> &new_positions);

  /**
   * @brief Update the direction positions in a set of templates
   *
   * This method returns a set of vectors representing templates. The
   * vectors themselves are built by updating the indices of those in
   * the input set according to the map of the new direction positions.
   * The returned templates are also brought to the canonical form.
   *
   * @param templates is the set of the templates whose direction position
   *               must be updated
   * @param new_positions is the map of new direction positions
   * @return the updated set of templates
   */
  static std::set<std::vector<unsigned int>>
  update_directions(const std::set<std::vector<unsigned int>> &templates,
                    const std::map<unsigned int, unsigned int> &new_positions);

  /**
   * @brief Collect directions in templates represented as vectors
   *
   * @param templates is the set of templates whose directions is aimed for
   * @return the set of directions mentioned in the templates
   */
  static std::set<unsigned int>
  collect_directions(const std::set<std::vector<unsigned int>> &templates);

  /**
   * @brief Add templates including non-mentioned directions
   *
   * @tparam T is the direction scalar type
   * @param[out] templates is a set of templates represented as index vectors
   * @param[in] directions is the direction vector
   * @param[in,out] missing_dirs is the set of directions that are not
   * mentioned in the templates
   */
  template<typename T>
  static void add_missing_templates(
      std::set<std::vector<unsigned int>> &templates,
      const std::vector<LinearAlgebra::Vector<T>> &directions,
      std::set<size_t> missing_dirs);

  /**
   * @brief Add templates including non-mentioned directions
   *
   * @tparam T is the direction scalar type
   * @param[out] templates is a set of templates represented as index vectors
   * @param[in] directions is the direction vector
   * @param[in] adaptive_directions is the set of adaptive directions
   * @param[in,out] missing_dirs is the set of directions that are not
   * mentioned in the templates
   */
  template<typename T>
  static void add_missing_templates(
      std::set<BundleTemplate> &templates,
      const std::vector<LinearAlgebra::Vector<T>> &directions,
      const std::set<size_t> &adaptive_directions,
      std::set<size_t> missing_dirs);

public:
  /**
   * @brief A copy constructor
   *
   * @param[in] orig is the copied bundle template
   */
  BundleTemplate(const BundleTemplate &orig);

  /**
   * @brief A move constructor
   *
   * @param[in] orig is the moved bundle template
   */
  BundleTemplate(BundleTemplate &&orig);

  /**
   * @brief A constructor for bundle templates
   *
   * @param dir_indices is a vector of bundle direction indices
   * @param adaptive_indices is the set of bundle adaptive direction indices
   */
  BundleTemplate(const std::vector<unsigned int> &dir_indices,
                 const std::set<size_t> &adaptive_indices);

  /**
   * @brief A move operator
   *
   * @param orig is the original template
   * @return a reference to the updated template
   */
  BundleTemplate &operator=(BundleTemplate &&orig);

  /**
   * @brief A copy operator
   *
   * @param orig is the original template
   * @return a reference to the updated template
   */
  BundleTemplate &operator=(const BundleTemplate &orig);

  /**
   * @brief Get the vector of the template direction indices
   *
   * @return the vector of the template direction indices
   */
  const std::vector<unsigned int> &direction_indices() const;

  /**
   * @brief Get the set of the template adaptive direction indices
   *
   * @return the set of the template adaptive direction indices
   */
  const std::set<size_t> &adaptive_indices() const;

  /**
   * @brief Test whether the template includes adaptive direction indices
   *
   * @return `true` if and only if the template includes adaptive
   *      direction indices
   */
  bool is_adaptive() const;

  /**
   * @brief Get the number of direction indices in the template
   *
   * @return the number of direction indices in the template
   */
  size_t dim() const;

  /**
   * @brief Get the i-th direction index in the template
   *
   * @param i is the position of the aimed direction index
   * @return a reference to the i-th direction index in the
   *        template
   */
  const unsigned int &operator[](const size_t &i) const;

  template<typename T, typename APPROX_TYPE>
  friend class Bundle;
};

/**
 * @brief Print a bundle template
 *
 * @param os is the output stream
 * @param bundle_template is the bundle template to be printed
 * @return a reference to the output stream
 */
std::ostream &operator<<(std::ostream &os,
                         const BundleTemplate &bundle_template);

template<>
struct std::less<BundleTemplate> {
  bool operator()(const BundleTemplate &a, const BundleTemplate &b) const;
};

// Implementation


template<typename T>
void BundleTemplate::add_missing_templates(
    std::set<std::vector<unsigned int>> &templates,
    const std::vector<LinearAlgebra::Vector<T>> &directions,
    std::set<size_t> missing_dirs)
{
  using namespace LinearAlgebra;
  using namespace LinearAlgebra::Dense;

  while (!missing_dirs.empty()) {
    // fill T with all the bundle directions starting
    // from those not included in a template
    Matrix<T> D = directions;
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

    BundleTemplate::canonize(new_template);

    templates.insert(new_template);
    for (const auto &idx: new_template) {
      missing_dirs.erase(idx);
    }
  }
}

template<typename T>
void BundleTemplate::add_missing_templates(
    std::set<BundleTemplate> &templates,
    const std::vector<LinearAlgebra::Vector<T>> &directions,
    const std::set<size_t> &adaptive_directions, std::set<size_t> missing_dirs)
{
  std::set<std::vector<unsigned int>> raw_templates;

  BundleTemplate::add_missing_templates(raw_templates, directions,
                                        missing_dirs);

  for (auto &raw_template: raw_templates) {
    templates.emplace(raw_template, adaptive_directions);
  }
}

inline BundleTemplate &BundleTemplate::operator=(BundleTemplate &&orig)
{
  std::swap(_dir_indices, orig._dir_indices);

  return *this;
}

inline BundleTemplate &BundleTemplate::operator=(const BundleTemplate &orig)
{
  _dir_indices = orig._dir_indices;

  return *this;
}

inline const std::vector<unsigned int> &
BundleTemplate::direction_indices() const
{
  return _dir_indices;
}

inline const std::set<size_t> &BundleTemplate::adaptive_indices() const
{
  return _adaptive_indices;
}

inline bool BundleTemplate::is_adaptive() const
{
  return _adaptive_indices.size() > 0;
}

inline size_t BundleTemplate::dim() const
{
  return _dir_indices.size();
}

inline const unsigned int &BundleTemplate::operator[](const size_t &i) const
{
  return _dir_indices[i];
}

#endif /* _BUNDLE_TEMPLATE_H_ */
