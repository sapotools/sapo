/**
 * @file LinearSystemSet.cpp
 * Represent and manipulate a set of linear systems
 * It can be used to represent a symbolic union of polytopes
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#include "LinearSystemSet.h"

#include <thread>
#include <shared_mutex>

/**
 * Constructor that instantiates an empty set
 */
LinearSystemSet::LinearSystemSet(): set() {}

/**
 * Constructor that instantiates a singleton set
 *
 * @param[in] ls element of the set
 */
LinearSystemSet::LinearSystemSet(const LinearSystem &ls)
{
  if (!ls.isEmpty()) {
    this->set.push_back(std::make_shared<LinearSystem>(ls.get_simplified()));
  }
}

/**
 * Constructor that instantiates a singleton set
 *
 * @param[in] ls pointer to the element of the set
 */
LinearSystemSet::LinearSystemSet(pointer ls)
{
  if (!ls->isEmpty()) {
    ls->simplify();
    this->set.push_back(ls);
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
    this->set.push_back(std::make_shared<LinearSystem>(ls));
  }
}

/**
 * Constructor that instantiates a set from a vector of sets
 *
 * @param[in] set vector of linear systems
 */
LinearSystemSet::LinearSystemSet(const container &set)
{
  for (unsigned int i = 0; i < set.size(); i++) {
    if (!set[i]->isEmpty()) {
      this->set.push_back(set[i]);
    }
  }
}

/**
 * A copy constructor for a linear system set
 *
 * @param[in] orig a linear system set
 */
LinearSystemSet::LinearSystemSet(const LinearSystemSet &orig)
{
  for (auto it = orig.set.cbegin(); it != orig.set.cend(); ++it) {
    this->set.push_back(std::make_shared<LinearSystem>(*(*it)));
  }
}


/**
 * Get the set of linear systems
 *
 * @returns the current collection of linear systems
 */
bool satisfiesOneIn(const LinearSystem &set, const LinearSystemSet &S)
{

#if MINIMIZE_LS_SET_REPRESENTATION
  class ThreadResult {
    mutable std::shared_timed_mutex mutex;
    bool value;
  public:
    ThreadResult(): value(false)
    {}

    bool get() const 
    {
      std::shared_lock<std::shared_timed_mutex> rlock(mutex, std::defer_lock);

      return value;
    }

    void set(const bool& value) 
    {
      std::unique_lock<std::shared_timed_mutex> wlock(mutex, std::defer_lock);

      this->value = value;
    }
  };

  ThreadResult result;
  
  auto check_and_update = [&result, &set](const LinearSystem& ls) { 
    if (!result.get() && set.satisfies(ls)) {
      result.set(true);
    }
  };

  std::vector<std::thread> threads;
  for (LinearSystemSet::const_iterator it = S.cbegin(); it != S.cend(); ++it) {
    threads.push_back(std::thread(check_and_update, std::ref(*it)));
  }

  for (std::thread & th : threads) {
    if (th.joinable())
        th.join();
  }

  return result.get();
#endif

  return false;
}

/**
 * Add a linear system to the set
 *
 * @param[in] ls linear system to add
 */
LinearSystemSet &LinearSystemSet::add(const LinearSystem &ls)
{
  return this->add(std::make_shared<LinearSystem>(ls));
}

/**
 * Add a linear system to the set
 *
 * @param[in] ls linear system to be added
 */
LinearSystemSet &LinearSystemSet::add(LinearSystem &&ls)
{
  return this->add(std::make_shared<LinearSystem>(ls));
}

/**
 * Add a linear system to the set
 *
 * @param[in] ls linear system to be added
 */
LinearSystemSet &LinearSystemSet::add(pointer ls)
{
  if (size() != 0 && set[0]->dim() != ls->dim()) {
    std::cerr << "Adding to a linear system set a "
              << "linear system having a different dimension" << std::endl;
  }

  if ((!ls->isEmpty()) && (!satisfiesOneIn(*ls, *this))) {
    this->set.push_back(ls);
  }

  return *this;
}

// TODO: parallelize the following method
LinearSystemSet &LinearSystemSet::simplify()
{
  class ThreadResult {
    mutable std::shared_timed_mutex mutex;
    unsigned int non_empty;
    std::map<unsigned int, unsigned int> new_position;
  public:
    ThreadResult(): non_empty(0)
    {}

    unsigned int get_non_empty() const 
    {
      std::shared_lock<std::shared_timed_mutex> read_lock(mutex, std::defer_lock);

      return non_empty;
    }

    void set_non_empty(const unsigned int& index) 
    {
      std::unique_lock<std::shared_timed_mutex> write_lock(mutex, std::defer_lock);

      this->new_position[non_empty++] = index;
    }

    const unsigned int& old_pos(const unsigned int& new_index) const  
    {
      std::shared_lock<std::shared_timed_mutex> read_lock(mutex, std::defer_lock);

      return this->new_position.at(new_index);
    }
  };

  ThreadResult result;

  auto test_emptiness_and_simplify = [&result](LinearSystem& ls, const unsigned int& i) {
    if (!ls.isEmpty()) {
      ls.simplify();
      result.set_non_empty(i);
    } 
  };

  std::vector<bool> non_empty_vect(set.size());
  std::vector<std::thread> threads;
  for (unsigned int i=0; i<set.size(); ++i) {
    /*
    threads.push_back(std::thread(test_emptiness_and_simplify, std::ref(*set[i]), std::ref(i)));
    /*/
    test_emptiness_and_simplify(std::ref(*set[i]), std::ref(i));
    //*/
  }

  for (std::thread & th : threads) {
    if (th.joinable())
        th.join();
  }

  container new_set(result.get_non_empty());
  for (unsigned int i=0; i<new_set.size(); ++i) {
    new_set[i] = set[result.old_pos(i)];
  }

  swap(set, new_set);

  return *this;
}

LinearSystemSet LinearSystemSet::get_a_finer_covering() const
{
  LinearSystemSet covering;

  for (auto it = std::begin(set); it != std::end(set); ++it) {
    covering.unionWith((*it)->get_a_finer_covering());
  }

  return covering;
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

  for (auto t_it = std::begin(A.set); t_it != std::end(A.set); ++t_it) {
    for (auto s_it = std::begin(B.set); s_it != std::end(B.set); ++s_it) {
      LinearSystem I = intersection(*(*t_it), *(*s_it));

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

  for (auto t_it = std::begin(A.set); t_it != std::end(A.set); ++t_it) {
    LinearSystem I = intersection(*(*t_it), B);

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
  for (container::const_iterator it = std::cbegin(LSset.set);
       it != std::cend(LSset.set); ++it) {
    this->add(std::make_shared<LinearSystem>(*(*it)));
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
  for (container::iterator it = std::begin(LSset.set);
       it != std::end(LSset.set); ++it) {
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
    vol = vol + this->set[i]->volBoundingBox();
  }
  return vol;
}

unsigned int LinearSystemSet::dim() const
{
  if (this->set.empty()) {
    return 0;
  }

  return (this->set[0])->dim();
}

/**
 * Check if the current set is empty
 *
 * @returns true if the set is empty
 */
bool LinearSystemSet::isEmpty() const
{
  if (this->set.empty()) {
    return true;
  }

  for (auto it = std::begin(set); it != std::end(set); ++it) {
    if (!(*it)->isEmpty()) {
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
  for (auto it = std::begin(set); it != std::end(set); ++it) {
    (*it)->plotRegion(os, color);
  }
}

LinearSystemSet::~LinearSystemSet()
{
  // TODO Auto-generated destructor stub
}

std::ostream &operator<<(std::ostream &out, const LinearSystemSet &ls)
{
  using namespace std;

  if (ls.size() == 0) {
    out << "---- empty set ----" << endl;
  } else {
    out << "--------------";
    const LinearSystemSet::const_iterator last(ls.cend() - 1);
    for (auto it = ls.cbegin(); it != last; ++it) {
      out << endl << *it << endl;
    }

    out << endl << *last;
  }

  return out;
}
