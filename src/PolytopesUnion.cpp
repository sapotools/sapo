/**
 * @file PolytopesUnion.cpp
 * Represent and manipulate a set of polytopes
 * It can be used to represent a symbolic union of polytopes
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#include "PolytopesUnion.h"

using namespace std;

/**
 * Constructor that instantiates an empty set
 */
PolytopesUnion::PolytopesUnion() {}

/**
 * Constructor that instantiates a singleton
 *
 * @param[in] P is a polytope
 */
PolytopesUnion::PolytopesUnion(const Polytope &P)
{
  if (!P.is_empty()) {
    this->push_back(P);

    this->back().simplify();
  }
}

/**
 * Constructor that instantiates a singleton
 *
 * @param[in] P is a polytope
 */

PolytopesUnion::PolytopesUnion(Polytope &&P)
{
  if (!P.is_empty()) {
    P.simplify();
    this->push_back(P);
  }
}

/**
 * A copy constructor for a polytopes union
 *
 * @param[in] orig a polytopes union
 */
PolytopesUnion::PolytopesUnion(const PolytopesUnion &orig):
    std::vector<Polytope>(orig)
{
}

PolytopesUnion::PolytopesUnion(PolytopesUnion &&orig)
{
  std::swap(*(static_cast<std::vector<Polytope> *>(this)),
            *(static_cast<std::vector<Polytope> *>(&orig)));
}

PolytopesUnion &PolytopesUnion::operator=(const PolytopesUnion &orig)
{
  resize(0);

  for (auto it = std::cbegin(orig); it != std::cend(orig); ++it) {
    this->push_back(*it);
  }

  return *this;
}

PolytopesUnion &PolytopesUnion::operator=(PolytopesUnion &&orig)
{
  std::swap(*(static_cast<std::vector<Polytope> *>(this)),
            *(static_cast<std::vector<Polytope> *>(&orig)));

  return *this;
}

bool PolytopesUnion::contains(const Polytope &P)
{
  for (auto it = std::cbegin(*this); it != std::cend(*this); ++it) {
    if (it->contains(P)) {
      return true;
    }
  }

  return false;
}

/**
 * Add a polytope to the set
 *
 * @param[in] P is the polytope to be added
 */
PolytopesUnion &PolytopesUnion::add(const Polytope &P)
{
  if (size() != 0 && (front().dim() != P.dim())) {
    std::cerr << "Adding to a polytopes union a "
              << "polytope having a different dimension" << std::endl;
  }

  if (!(P.is_empty() || this->contains(P))) {
    this->push_back(P);
  }

  return *this;
}

/**
 * Add a polytope to the set
 *
 * @param[in] ls polytope to be added
 */
PolytopesUnion &PolytopesUnion::add(Polytope &&P)
{
  if (size() != 0 && (front().dim() != P.dim())) {
    std::cerr << "Adding to a polytopes union a "
              << "polytope having a different dimension" << std::endl;
  }

  if (!(P.is_empty() || this->contains(P))) {
    this->push_back(P);
  }

  return *this;
}

PolytopesUnion &PolytopesUnion::simplify()
{
  for (auto it = std::begin(*this); it != std::end(*this); ++it) {
    it->simplify();
  }

  return *this;
}

/**
 * Compute the intersection of two polytopes
 *
 * @param[in] A is a polytope
 * @param[in] B is a polytope
 * @return the polytopes union representing the intersection
 *     of the two parameters.
 */
PolytopesUnion intersect(const PolytopesUnion &A, const PolytopesUnion &B)
{
  PolytopesUnion result;

  for (auto t_it = std::begin(A); t_it != std::end(A); ++t_it) {
    for (auto s_it = std::begin(B); s_it != std::end(B); ++s_it) {
      Polytope I = intersect(*t_it, *s_it);

      if (!I.is_empty()) {
        result.add(I);
      }
    }
  }

  return result;
}

/**
 * Compute the intersection between a polytopes unions and a polytope
 *
 * @param[in] A is a polytopes union
 * @param[in] B is a polytope
 * @return the polytopes union that represents the set of values satisfying
 * both the parameters
 */
PolytopesUnion intersect(const PolytopesUnion &A, const Polytope &B)
{
  PolytopesUnion result;

  for (auto t_it = std::begin(A); t_it != std::end(A); ++t_it) {
    Polytope I = intersect(*t_it, B);

    if (!I.is_empty()) {
      result.add(I);
    }
  }

  return result;
}

/**
 * Union of polytopes unions
 *
 * @param[in] Pu the polytopes union to be adjoint
 * @returns a reference to the updated union.
 */
PolytopesUnion &PolytopesUnion::add(const PolytopesUnion &Pu)
{
  for (auto it = std::cbegin(Pu); it != std::cend(Pu); ++it) {
    this->add(*it);
  }

  return *this;
}

PolytopesUnion unite(const PolytopesUnion &A, const PolytopesUnion &B)
{
  return PolytopesUnion(A).add(B);
}

/**
 * Union of polytopes unions
 *
 * @param[in] Pu the polytopes union to be adjoint
 * @returns a reference to the updated union.
 */
PolytopesUnion &PolytopesUnion::add(PolytopesUnion &&Pu)
{
  for (auto it = std::begin(Pu); it != std::end(Pu); ++it) {
    this->add(*it);
  }

  return *this;
}

/**
 * Compute the union of two polytopes
 *
 * @param[in] A is a polytope
 * @param[in] B is a polytope
 * @return the union of the two polytopes.
 */
PolytopesUnion unite(const Polytope &A, const Polytope &B)
{
  PolytopesUnion result(A);

  result.add(B);

  return result;
}

/**
 * Sum of volumes of boxes containing the sets
 *
 * @returns sum of bounding boxes
 */
double PolytopesUnion::volume_of_bounding_boxes() const
{

  double vol = 0;
  for (unsigned int i = 0; i < this->size(); i++) {
    vol = vol + (*this)[i].bounding_box_volume();
  }
  return vol;
}

unsigned int PolytopesUnion::dim() const
{
  if (this->empty()) {
    return 0;
  }

  return front().dim();
}

/**
 * Check if the current set is empty
 *
 * @returns true if the set is empty
 */
bool PolytopesUnion::is_empty() const
{
  if (this->empty()) {
    return true;
  }

  for (auto it = std::begin(*this); it != std::end(*this); ++it) {
    if (!it->is_empty()) {
      return false;
    }
  }

  return true;
}

/**
 * Print the polytope in Matlab format (for plotregion script)
 *
 * @param[in] os is the output stream
 * @param[in] color color of the polytope to plot
 */
void PolytopesUnion::plotRegion(std::ostream &os, const char color) const
{
  for (auto it = std::begin(*this); it != std::end(*this); ++it) {
    it->plotRegion(os, color);
  }
}

PolytopesUnion::~PolytopesUnion()
{
  // TODO Auto-generated destructor stub
}

bool every_set_is_empty(const std::list<PolytopesUnion> &ps_list)
{
  for (auto ps_it = std::begin(ps_list); ps_it != std::end(ps_list); ++ps_it) {
    if (!ps_it->is_empty()) {
      return false;
    }
  }

  return true;
}
