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
LinearSystemSet::LinearSystemSet(const LinearSystem& ls){
	if(!ls.isEmpty()){
		this->set.push_back(new LinearSystem(ls.get_simplified()));
	}
}

/**
 * Constructor that instantiates a singleton set
 *
 * @param[in] ls element of the set
 */
LinearSystemSet::LinearSystemSet(LinearSystem* ls){
	if(!ls->isEmpty()){
		ls->simplify();
		this->set.push_back(ls);
	}
}

/**
 * Constructor that instantiates a set from a vector of sets
 *
 * @param[in] set vector of linear systems
 */
LinearSystemSet::LinearSystemSet(const vector<LinearSystem*>& set){

	for(unsigned int i=0; i<set.size(); i++){
		if(!set[i]->isEmpty()){
			this->set.push_back(set[i]);
		}
	}
}

/**
 * Get the set of linear systems
 *
 * @returns the current collection of linear systems
 */

bool satisfiesOneOf(const LinearSystem &set, const std::vector<LinearSystem*>& S) {

#if MINIMIZE_LS_SET_REPRESENTATION
	for (std::vector<LinearSystem*>::const_iterator it=std::begin(S); it!=std::end(S); ++it) {
		if (set.satisfies(*(*it))) {
			return true;
		}
	}
#endif

	return false;
}

/**
 * Add a linear system to the set
 *
 * @param[in] LS linear system to add
 */
void LinearSystemSet::add(LinearSystem *ls)
{
	if ((!ls->isEmpty()) && (!satisfiesOneOf(*ls, this->set))) {
		this->set.push_back(ls);
	}
}

LinearSystemSet* LinearSystemSet::get_a_finer_covering() const
{
	LinearSystemSet* covering = new LinearSystemSet();

	for (std::vector<LinearSystem *>::const_iterator it=std::begin(set);
													 it!=std::end(set); ++it) {
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
	vector<LinearSystem*> set = LSset->get_set();

	for(int i=0; i<(signed)this->set.size(); i++){
		for(int j=0; j<(signed)set.size(); j++){
			LinearSystem *intLS = new LinearSystem(this->set[i]->getIntersectionWith(*set[j])); // intersect
			result->add(intLS);
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
LinearSystemSet& LinearSystemSet::unionWith(LinearSystemSet *LSset){
	for (std::vector<LinearSystem*>::iterator it=std::begin(LSset->set); 
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
LinearSystemSet& LinearSystemSet::boundedUnionWith(LinearSystemSet *LSset, const int bound){

	if(this->size() > bound){
		cout<<"LinearSystemSet::boundedUnionWith : size of actual box larger than bound";
		exit (EXIT_FAILURE);
	}

	vector<LinearSystem*> set = LSset->get_set();
	int set_card = set.size();
	int iters = min(bound-this->size(),set_card);

	for(int i=0; i<iters; i++){
		this->set.push_back(set[i]);
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
	for(int i=0; i<this->size(); i++){
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
LinearSystem* LinearSystemSet::at(int i) {
	return this->set[i];
}

/**
 * Get the i-th linear system
 *
 * @param[in] index of the linear system to fetch
 * @returns i-th linear system
 */
const LinearSystem* LinearSystemSet::at(int i) const {
	return this->set[i];
}

/**
 * Check if the current set is empty
 *
 * @returns true if the set is empty
 */
bool LinearSystemSet::isEmpty() const{
	return this->set.empty();
}

/**
 * Print the set of linear systems
 */
void LinearSystemSet::print() const {

	if( this->set.size() <= 0){
		cout<<"--- empty set ----" << endl;
	}else{
		for(int i=0; i<(signed)this->set.size(); i++){
			cout << "--------------" << endl 
			     << *(this->set[i]) << endl << endl;
		}
	}

}

LinearSystemSet::~LinearSystemSet() {
	// TODO Auto-generated destructor stub
}

