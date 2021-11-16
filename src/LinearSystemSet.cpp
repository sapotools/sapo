/**
 * @file LinearSystemSet.cpp
 * Represent and manipulate a set of linear systems
 * It can be used to represent a symbolic union of polytopes
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#include "LinearSystemSet.h"

using namespace std;
using namespace GiNaC;

/**
 * Constructor that instantiates an empty set
 */
LinearSystemSet::LinearSystemSet() {}

/**
 * Constructor that instantiates a singleton set
 *
 * @param[in] ls element of the set
 */
LinearSystemSet::LinearSystemSet(const LinearSystem &ls)
{
  if (!ls.isEmpty()) {
    this->push_back(ls.get_simplified());
  }
}

/**
 * Constructor that instantiates a singleton set
 *
 * @param[in] orig a linear system
 */

LinearSystemSet::LinearSystemSet(LinearSystem &&ls)
{
  if (!ls.isEmpty()) {
    ls.simplify();
    this->push_back(ls);
  }
}

/**
 * A copy constructor for a linear system set
 *
 * @param[in] orig a linear system set
 */
LinearSystemSet::LinearSystemSet(const LinearSystemSet &orig):
    std::vector<LinearSystem>(orig)
{
}

LinearSystemSet::LinearSystemSet(LinearSystemSet &&orig)
{
  std::swap(*(static_cast<std::vector<LinearSystem> *>(this)),
            *(static_cast<std::vector<LinearSystem> *>(&orig)));
}

LinearSystemSet &LinearSystemSet::operator=(const LinearSystemSet &orig)
{
  resize(0);

  for (auto it = std::cbegin(orig); it != std::cend(orig); ++it) {
    this->push_back(*it);
  }

  return *this;
}

LinearSystemSet &LinearSystemSet::operator=(LinearSystemSet &&orig)
{
  std::swap(*(static_cast<std::vector<LinearSystem> *>(this)),
            *(static_cast<std::vector<LinearSystem> *>(&orig)));

  return *this;
}

/**
 * Get the set of linear systems
 *
 * @returns the current collection of linear systems
 */

bool satisfiesOneIn(const LinearSystem &set, const LinearSystemSet &S)
{

#if MINIMIZE_LS_SET_REPRESENTATION
  for (LinearSystemSet::const_iterator it = S.begin(); it != S.end(); ++it) {
    if (set.satisfies(*it)) {
      return true;
    }
  }
#endif

  return false;
}

/**
 * Add a linear system to the set
 *
 * @param[in] ls linear system to be added
 */
LinearSystemSet &LinearSystemSet::add(const LinearSystem &ls)
{
  if (size() != 0 && (front().dim() != ls.dim())) {
    std::cerr << "Adding to a linear system set a "
              << "linear system having a different dimension" << std::endl;
  }

  if ((!ls.isEmpty()) && (!satisfiesOneIn(ls, *this))) {
    this->push_back(ls);
  }

  return *this;
}

/**
 * Add a linear system to the set
 *
 * @param[in] ls linear system to be added
 */
LinearSystemSet &LinearSystemSet::add(LinearSystem &&ls)
{
  if (size() != 0 && (front().dim() != ls.dim())) {
    std::cerr << "Adding to a linear system set a "
              << "linear system having a different dimension" << std::endl;
  }

  if ((!ls.isEmpty()) && (!satisfiesOneIn(ls, *this))) {
    this->push_back(ls);
  }

  return *this;
}

LinearSystemSet &LinearSystemSet::simplify()
{
  for (auto it = std::begin(*this); it != std::end(*this); ++it) {
    it->simplify();
  }

  return *this;
}

/**
 * Compute the intersection of two linear systems
 *
 * @param[in] A is a linear system
 * @param[in] B is a linear system
 * @return the linear system that represents the set of values satisfying both
 *          the parameters.
 */
LinearSystemSet intersection(const LinearSystemSet &A,
                             const LinearSystemSet &B)
{
  LinearSystemSet result;

  for (auto t_it = std::begin(A); t_it != std::end(A); ++t_it) {
    for (auto s_it = std::begin(B); s_it != std::end(B); ++s_it) {
      LinearSystem I = intersection(*t_it, *s_it);

      if (!I.isEmpty()) {
        result.add(I);
      }
    }
  }

  return result;
}

/**
 * Compute the intersection between a linear system sets and a linear system
 *
 * @param[in] A is a linear system set
 * @param[in] B is a linear system
 * @return the linear system set that represents the set of values satisfying
 * both the parameters
 */
LinearSystemSet intersection(const LinearSystemSet &A, const LinearSystem &B)
{
  LinearSystemSet result;

  for (auto t_it = std::begin(A); t_it != std::end(A); ++t_it) {
    LinearSystem I = intersection(*t_it, B);

    if (!I.isEmpty()) {
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
LinearSystemSet &LinearSystemSet::unionWith(const LinearSystemSet &LSset)
{
  for (auto it = std::cbegin(LSset); it != std::cend(LSset); ++it) {
    this->add(*it);
  }

  return *this;
}

LinearSystemSet unionset(const LinearSystemSet &A, const LinearSystemSet &B)
{
  return LinearSystemSet(A).unionWith(B);
}

/**
 * Union of sets
 *
 * @param[in] LSset set to union with
 * @returns merged sets
 */
LinearSystemSet &LinearSystemSet::unionWith(LinearSystemSet &&LSset)
{
  for (auto it = std::begin(LSset); it != std::end(LSset); ++it) {
    this->add(*it);
  }

  return *this;
}

/**
 * Compute the union of two linear systems
 *
 * @param[in] A is a linear system
 * @param[in] B is a linear system
 * @return the linear system set that represents the set of values satisfying
 *          at least one of the parameters.
 */
LinearSystemSet unionset(const LinearSystem &A, const LinearSystem &B)
{
  LinearSystemSet result(A);

  result.add(B);

  return result;
}

/**
 * Union of two sets of linear systems up to bounded cardinality
 *
 * @param[in] LSset set to union with
 * @param[in] bound set size bound
 * @returns merged sets
 */
/* // TODO: remove this code
LinearSystemSet& LinearSystemSet::boundedUnionWith(LinearSystemSet *LSset,
const unsigned int bound) {

        if (this->size() > bound) {
                std::cerr << "LinearSystemSet::boundedUnionWith : size of
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
double LinearSystemSet::boundingVol() const
{

  double vol = 0;
  for (unsigned int i = 0; i < this->size(); i++) {
    vol = vol + (*this)[i].volBoundingBox();
  }
  return vol;
}

unsigned int LinearSystemSet::dim() const
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
bool LinearSystemSet::isEmpty() const
{
  if (this->empty()) {
    return true;
  }

  for (auto it = std::begin(*this); it != std::end(*this); ++it) {
    if (!it->isEmpty()) {
      return false;
    }
  }

  return true;
}

/**
 * Print the set of linear systems
 */
void LinearSystemSet::print() const
{
  std::cout << *this << std::endl;
}

/**
 * Print the linear system in Matlab format (for plotregion script)
 *
 * @param[in] os is the output stream
 * @param[in] color color of the polytope to plot
 */
void LinearSystemSet::plotRegion(std::ostream &os, const char color) const
{
  for (auto it = std::begin(*this); it != std::end(*this); ++it) {
    it->plotRegion(os, color);
  }
}

LinearSystemSet::~LinearSystemSet()
{
  // TODO Auto-generated destructor stub
}

bool every_set_is_empty(const std::list<LinearSystemSet> &lss_list)
{
  for (auto lss_it = std::begin(lss_list); lss_it != std::end(lss_list);
       ++lss_it) {
    if (!lss_it->isEmpty()) {
      return false;
    }
  }

  return true;
}
