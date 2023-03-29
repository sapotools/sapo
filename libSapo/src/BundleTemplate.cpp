/**
 * @file BundleTemplate.cpp
 * @author Alberto Casagrande <acasagrande@units.it>
 * @brief Represent bundle template
 * @version 0.1
 * @date 2023-03-29
 *
 * @copyright Copyright (c) 2023
 */

#include "BundleTemplate.h"

#include <string>
#include <algorithm> // std::sort, std::is_sorted

#include "ErrorHandling.h"

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

std::vector<unsigned int> &
BundleTemplate::canonize(std::vector<unsigned int> &bundle_template)
{
  if (!std::is_sorted(std::begin(bundle_template),
                      std::end(bundle_template))) {
    std::sort(std::begin(bundle_template), std::end(bundle_template));
  }

  return bundle_template;
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

std::set<unsigned int> BundleTemplate::collect_directions(
    const std::set<std::vector<unsigned int>> &templates)
{
  std::set<unsigned int> template_directions;
  for (const auto &bundle_template: templates) {
    template_directions.insert(std::begin(bundle_template),
                               std::end(bundle_template));
  }

  return template_directions;
}

std::set<std::vector<unsigned int>> BundleTemplate::update_directions(
    const std::set<std::vector<unsigned int>> &templates,
    const std::vector<size_t> &new_positions)
{
  std::set<std::vector<unsigned int>> new_templates;

  for (const auto &bundle_template: templates) {
    std::vector<unsigned int> new_template(bundle_template);
    for (auto &dir_id: new_template) {
      dir_id = new_positions[dir_id];
    }
    BundleTemplate::canonize(new_template);

    new_templates.insert(std::move(new_template));
  }

  return new_templates;
}

std::set<std::vector<unsigned int>> BundleTemplate::update_directions(
    const std::set<std::vector<unsigned int>> &templates,
    const std::map<unsigned int, unsigned int> &new_positions)
{
  std::set<std::vector<unsigned int>> new_templates;

  for (const auto &bundle_template: templates) {
    std::vector<unsigned int> new_template(bundle_template);
    for (auto &dir_id: new_template) {
      dir_id = new_positions.at(dir_id);
    }
    BundleTemplate::canonize(new_template);

    new_templates.insert(std::move(new_template));
  }

  return new_templates;
}
