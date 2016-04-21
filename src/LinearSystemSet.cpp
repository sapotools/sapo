/**
 * @file LinearSystemSet.cpp
 * Represent and manipulate a set of linear systems
 * It can be used to represent a symbolic union of polytopes
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#include "LinearSystemSet.h"

/**
 * Constructor that instantiates an empty set
 */
LinearSystemSet::LinearSystemSet(){
}

/**
 * Constructor that instantiates a singleton set
 *
 * @param[in] LS element of the set
 */
LinearSystemSet::LinearSystemSet(LinearSystem *LS){
	if(!LS->isEmpty()){
		this->set.push_back(LS);
	}
}

/**
 * Constructor that instantiates a set from a vector of sets
 *
 * @param[in] set vector of linear systems
 */
LinearSystemSet::LinearSystemSet(vector<LinearSystem*> set){

	for(int i=0; i<(signed)set.size(); i++){
		if(!set[i]->isEmpty()){
			this->set.push_back(set[i]);
		}
	}
}

/**
 * Get the set of linear systems
 *
 * @returns actual collection of linear systems
 */
vector<LinearSystem*> LinearSystemSet::getSet(){
	return this->set;
}

/**
 * Add a linear system to the set
 *
 * @param[in] LS linear system to add
 */
void LinearSystemSet::add(LinearSystem *LS){
	if(!LS->isEmpty()){
		this->set.push_back(LS);
	}
}

/**
 * Intersect to sets of linear systems
 *
 * @param[in] LSset set to intersect with
 * @returns intersected sets
 */
LinearSystemSet* LinearSystemSet::intersectWith(LinearSystemSet *LSset){

	vector<LinearSystem*> intSet; // new intersection set
	vector<LinearSystem*> set = LSset->getSet();

	for(int i=0; i<(signed)this->set.size(); i++){
		for(int j=0; j<(signed)set.size(); j++){
			LinearSystem *intLS = this->set[i]->appendLinearSystem(set[j]); // intersect
			if(!intLS->isEmpty()){
				intSet.push_back(intLS);	// add inteserction
			}
		}
	}
	return new LinearSystemSet(intSet);
}

/**
 * Union of sets
 *
 * @param[in] LSset set to union with
 * @returns merged sets
 */
LinearSystemSet* LinearSystemSet::unionWith(LinearSystemSet *LSset){

	vector<LinearSystem*> uniSet = this->set; 		// new union set
	vector<LinearSystem*> set = LSset->getSet();

	for(int i=0; i<(signed)set.size(); i++){
		uniSet.push_back(set[i]);
	}

	return new LinearSystemSet(uniSet);

}

/**
 * Union of two sets of linear systems up to bounded cardinality
 *
 * @param[in] LSset set to union with
 * @param[in] bound set size bound
 * @returns merged sets
 */
LinearSystemSet* LinearSystemSet::boundedUnionWith(LinearSystemSet *LSset, int bound){

	if(this->size() > bound){
		cout<<"LinearSystemSet::boundedUnionWith : size of actual box larger than bound";
		exit (EXIT_FAILURE);
	}

	vector<LinearSystem*> uniSet = this->set; 		// new union set
	vector<LinearSystem*> set = LSset->getSet();
	int set_card = set.size();
	int iters = min(bound-this->size(),set_card);

	for(int i=0; i<iters; i++){
		uniSet.push_back(set[i]);
	}

	return new LinearSystemSet(uniSet);

}

/**
 * Sum of volumes of boxes containing the sets
 *
 * @returns sum of bounding boxes
 */
double LinearSystemSet::boundingVol(){

	double vol = 0;
	for(int i=0; i<this->size(); i++){
		vol = vol + this->at(i)->volBoundingBox();
	}
	return vol;

}

/**
 * Get the size of this set, i.e,
 * the number of linear systems
 *
 * @returns number of linear systems in the set
 */
int LinearSystemSet::size() {
	return this->set.size();
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
 * Check if the current set is empty
 *
 * @returns true if the set is empty
 */
bool LinearSystemSet::isEmpty() {
	return this->set.empty();
}

/**
 * Print the set of linear systems
 */
void LinearSystemSet::print() {

	if( this->set.size() <= 0){
		cout<<"--- empty set ----\n";
	}else{
		for(int i=0; i<(signed)this->set.size(); i++){
			cout<<"--------------\n";
			this->set[i]->print();
		}
	}

}

LinearSystemSet::~LinearSystemSet() {
	// TODO Auto-generated destructor stub
}

