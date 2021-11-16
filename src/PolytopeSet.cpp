/**
 * @file PolytopeSet.cpp
 * Represent and manipulate a set of polytopes
 * It can be used to represent a symbolic union of polytopes
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#include "PolytopeSet.h"

using namespace std;
using namespace GiNaC;

/**
 * Constructor that instantiates an empty set
 */
PolytopeSet::PolytopeSet() {}

/**
 * Constructor that instantiates a singleton
 *
 * @param[in] P is a polytope
 */
PolytopeSet::PolytopeSet(const Polytope &P)
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

PolytopeSet::PolytopeSet(Polytope &&P)
{
  if (!P.is_empty()) {
    P.simplify();
    this->push_back(P);
  }
}

/**
 * A copy constructor for a polytope set
 *
 * @param[in] orig a polytope set
 */
PolytopeSet::PolytopeSet(const PolytopeSet &orig):
    std::vector<Polytope>(orig)
{
}

PolytopeSet::PolytopeSet(PolytopeSet &&orig)
{
  std::swap(*(static_cast<std::vector<Polytope> *>(this)),
            *(static_cast<std::vector<Polytope> *>(&orig)));
}

PolytopeSet &PolytopeSet::operator=(const PolytopeSet &orig)
{
  resize(0);

  for (auto it = std::cbegin(orig); it != std::cend(orig); ++it) {
    this->push_back(*it);
  }

  return *this;
}

PolytopeSet &PolytopeSet::operator=(PolytopeSet &&orig)
{
  std::swap(*(static_cast<std::vector<Polytope> *>(this)),
            *(static_cast<std::vector<Polytope> *>(&orig)));

  return *this;
}

/**
 * Get the set of polytopes
 *
 * @returns the current collection of polytopes
 */

bool satisfiesOneIn(const Polytope &set, const PolytopeSet &S)
{

#if MINIMIZE_LS_SET_REPRESENTATION
  for (PolytopeSet::const_iterator it = std::cbegin(S); it != std::cend(S); ++it) {
    if (it->contains(set)) {
      return true;
    }
  }
#endif

  return false;
}

/**
 * Add a polytope to the set
 *
 * @param[in] P is the polytope to be added
 */
PolytopeSet &PolytopeSet::add(const Polytope &P)
{
  if (size() != 0 && (front().dim() != P.dim())) {
    std::cerr << "Adding to a polytope set a "
              << "polytope having a different dimension" << std::endl;
  }

  if ((!P.is_empty()) && (!satisfiesOneIn(P, *this))) {
    this->push_back(P);
  }

  return *this;
}

/**
 * Add a polytope to the set
 *
 * @param[in] ls polytope to be added
 */
PolytopeSet &PolytopeSet::add(Polytope &&ls)
{
  if (size() != 0 && (front().dim() != ls.dim())) {
    std::cerr << "Adding to a polytope set a "
              << "polytope having a different dimension" << std::endl;
  }

  if ((!ls.is_empty()) && (!satisfiesOneIn(ls, *this))) {
    this->push_back(ls);
  }

  return *this;
}

PolytopeSet &PolytopeSet::simplify()
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
 * @return the polytope set representing the intersection
 *     of the two parameters.
 */
PolytopeSet intersection(const PolytopeSet &A,
                             const PolytopeSet &B)
{
  PolytopeSet result;

  for (auto t_it = std::begin(A); t_it != std::end(A); ++t_it) {
    for (auto s_it = std::begin(B); s_it != std::end(B); ++s_it) {
      Polytope I = intersection(*t_it, *s_it);

      if (!I.is_empty()) {
        result.add(I);
      }
    }
  }

  return result;
}

/**
 * Compute the intersection between a polytope sets and a polytope
 *
 * @param[in] A is a polytope set
 * @param[in] B is a polytope
 * @return the polytope set that represents the set of values satisfying
 * both the parameters
 */
PolytopeSet intersection(const PolytopeSet &A, const Polytope &B)
{
  PolytopeSet result;

  for (auto t_it = std::begin(A); t_it != std::end(A); ++t_it) {
    Polytope I = intersection(*t_it, B);

    if (!I.is_empty()) {
      result.add(I);
    }
  }

  return result;
}

/**
 * Union of sets
 *
 * @param[in] LSset set to union with
 * @returns merged sets
 */
PolytopeSet &PolytopeSet::unionWith(const PolytopeSet &LSset)
{
  for (auto it = std::cbegin(LSset); it != std::cend(LSset); ++it) {
    this->add(*it);
  }

  return *this;
}

PolytopeSet unionset(const PolytopeSet &A, const PolytopeSet &B)
{
  return PolytopeSet(A).unionWith(B);
}

/**
 * Union of sets
 *
 * @param[in] LSset set to union with
 * @returns merged sets
 */
PolytopeSet &PolytopeSet::unionWith(PolytopeSet &&LSset)
{
  for (auto it = std::begin(LSset); it != std::end(LSset); ++it) {
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
PolytopeSet unionset(const Polytope &A, const Polytope &B)
{
  PolytopeSet result(A);

  result.add(B);

  return result;
}

/**
 * Union of two sets of polytopes up to bounded cardinality
 *
 * @param[in] LSset set to union with
 * @param[in] bound set size bound
 * @returns merged sets
 */
/* // TODO: remove this code
PolytopeSet& PolytopeSet::boundedUnionWith(PolytopeSet *LSset,
const unsigned int bound) {

        if (this->size() > bound) {
                std::cerr << "PolytopeSet::boundedUnionWith : size of
actual box larger than bound" << std::endl; exit (EXIT_FAILURE);
        }

        int iters = min(bound-this->size(), (unsigned int)LSset->size());

        container::iterator it = std::begin(set);
        for (int i=0; i<iters&&it!=std::end(set); i++) {
                this->set.push_back(*it);
        }

        return *this;
}
*/

/**
 * Sum of volumes of boxes containing the sets
 *
 * @returns sum of bounding boxes
 */
double PolytopeSet::volume_of_bounding_boxes() const
{

  double vol = 0;
  for (unsigned int i = 0; i < this->size(); i++) {
    vol = vol + (*this)[i].bounding_box_volume();
  }
  return vol;
}

unsigned int PolytopeSet::dim() const
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
bool PolytopeSet::is_empty() const
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
 * Print the set of polytopes
 */
void PolytopeSet::print() const
{
  std::cout << *this << std::endl;
}

/**
 * Print the polytope in Matlab format (for plotregion script)
 *
 * @param[in] os is the output stream
 * @param[in] color color of the polytope to plot
 */
void PolytopeSet::plotRegion(std::ostream &os, const char color) const
{
  for (auto it = std::begin(*this); it != std::end(*this); ++it) {
    it->plotRegion(os, color);
  }
}

PolytopeSet::~PolytopeSet()
{
  // TODO Auto-generated destructor stub
}

bool every_set_is_empty(const std::list<PolytopeSet> &ps_list)
{
  for (auto ps_it = std::begin(ps_list); ps_it != std::end(ps_list);
       ++ps_it) {
    if (!ps_it->is_empty()) {
      return false;
    }
  }

  return true;
}
