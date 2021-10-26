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
LinearSystemSet::LinearSystemSet(): set() {}

/**
 * Constructor that instantiates a singleton set
 *
 * @param[in] ls element of the set
 */
LinearSystemSet::LinearSystemSet(const LinearSystem& ls)
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
 * Constructor that instantiates a set from a vector of sets
 *
 * @param[in] set vector of linear systems
 */
LinearSystemSet::LinearSystemSet(const container& set)
{
	for (unsigned int i=0; i<set.size(); i++) {
		if (!set[i]->isEmpty()) {
			this->set.push_back(set[i]);
		}
	}
}

/**
 * Get the set of linear systems
 *
 * @returns the current collection of linear systems
 */

bool satisfiesOneIn(const LinearSystem &set, const LinearSystemSet& S)
{

#if MINIMIZE_LS_SET_REPRESENTATION
	for (LinearSystemSet::const_iterator it=S.cbegin(); it!=S.cend(); ++it) {
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
 * @param[in] ls linear system to add
 */
LinearSystemSet& LinearSystemSet::add(const LinearSystem& ls)
{
	return this->add(std::make_shared<LinearSystem>(ls));
}

/**
 * Add a linear system to the set
 *
 * @param[in] ls linear system to be added
 */
LinearSystemSet& LinearSystemSet::add(LinearSystem&& ls)
{
	return this->add(std::make_shared<LinearSystem>(ls));
}

/**
 * Add a linear system to the set
 *
 * @param[in] ls linear system to be added
 */
LinearSystemSet& LinearSystemSet::add(pointer ls)
{
	if (size()!=0 && set[0]->dim()!=ls->dim()) {
		std::cerr << "Adding to a linear system set a "
				  << "linear system having a different dimension" << std::endl;
	}

	if ((!ls->isEmpty()) && (!satisfiesOneIn(*ls, *this))) {
		this->set.push_back(ls);
	}

	return *this;
}

LinearSystemSet& LinearSystemSet::simplify()
{
	for (auto it=std::begin(set); it!=std::end(set); ++it) {
		(*it)->simplify();
	}

	return *this;
}

LinearSystemSet* LinearSystemSet::get_a_finer_covering() const
{
	LinearSystemSet* covering = new LinearSystemSet();

	for (auto it=std::begin(set); it!=std::end(set); ++it) {
		covering->unionWith((*it)->get_a_finer_covering());
	}

	return covering;
}

/**
 * Intersect to sets of linear systems
 *
 * @param[in] LSset set to intersect with
 * @returns intersected sets
 */
LinearSystemSet* LinearSystemSet::getIntersectionWith(const LinearSystemSet *LSset) const {

	LinearSystemSet* result = new LinearSystemSet();
	for (auto t_it=std::begin(set); t_it!=std::end(set); ++t_it) {
		for (auto s_it=std::begin(LSset->set); s_it!=std::end(LSset->set); ++s_it) {
			result->add(std::make_shared<LinearSystem>((*t_it)->getIntersectionWith(*(*s_it))));
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
LinearSystemSet& LinearSystemSet::unionWith(LinearSystemSet *LSset) {
	for (container::iterator it=std::begin(LSset->set);
	                         it!=std::end(LSset->set); ++it) {
		this->add(*it);
	}

	return *this;
}

/**
 * Union of two sets of linear systems up to bounded cardinality
 *
 * @param[in] LSset set to union with
 * @param[in] bound set size bound
 * @returns merged sets
 */
LinearSystemSet& LinearSystemSet::boundedUnionWith(LinearSystemSet *LSset, const unsigned int bound) {

	if (this->size() > bound) {
		std::cerr << "LinearSystemSet::boundedUnionWith : size of actual box larger than bound" << std::endl;
		exit (EXIT_FAILURE);
	}

	int iters = min(bound-this->size(), (unsigned int)LSset->size());

	container::iterator it = std::begin(set);
	for (int i=0; i<iters&&it!=std::end(set); i++) {
		this->set.push_back(*it);
	}

	return *this;
}

/**
 * Sum of volumes of boxes containing the sets
 *
 * @returns sum of bounding boxes
 */
double LinearSystemSet::boundingVol() const{

	double vol = 0;
	for (unsigned int i=0; i<this->size(); i++) {
		vol = vol + this->set[i]->volBoundingBox();
	}
	return vol;

}

/**
 * Get the i-th linear system
 *
 * @param[in] index of the linear system to fetch
 * @returns i-th linear system
 */
/*
LinearSystem* LinearSystemSet::at(int i) {
	return this->set[i];
}
*/

/**
 * Get the i-th linear system
 *
 * @param[in] index of the linear system to fetch
 * @returns i-th linear system
 */
/*
const LinearSystem* LinearSystemSet::at(int i) const {
	return this->set[i];
}
*/

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
bool LinearSystemSet::isEmpty() const{
	if (this->set.empty()) {
		return true;
	}

	for (auto it=std::begin(set); it!=std::end(set); ++it) {
		if (!(*it)->isEmpty()) {
			return false;
		}
	}

	return true;
}

/**
 * Print the set of linear systems
 */
void LinearSystemSet::print() const {
	std::cout << *this << std::endl;
}

/**
 * Print the linear system in Matlab format (for plotregion script)
 * 
 * @param[in] os is the output stream
 * @param[in] color color of the polytope to plot
 */
void LinearSystemSet::plotRegion(std::ostream& os, const char color) const
{
	for (auto it=std::begin(set); it!=std::end(set); ++it) {
		(*it)->plotRegion(os, color);
	}
}

LinearSystemSet::~LinearSystemSet() {
	// TODO Auto-generated destructor stub
}

std::ostream& operator<<(std::ostream& out, const LinearSystemSet& ls)
{
	if ( ls.size() == 0) {
		out<<"---- empty set ----" << endl;
	}else{
		out << "--------------";
		const LinearSystemSet::const_iterator last(ls.cend()-1);
		for (auto it=ls.cbegin(); it!=last; ++it) {
			out << endl << *it  << endl;
		}

		out << endl << *last;
	}

	return out;
}