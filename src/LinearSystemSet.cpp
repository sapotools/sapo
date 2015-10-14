/*
 * LinearSystemSet.cpp
 *
 *  Created on: Nov 17, 2014
 *      Author: Tommaso Dreossi
 */

#include "LinearSystemSet.h"

LinearSystemSet::LinearSystemSet(){
}

LinearSystemSet::LinearSystemSet(LinearSystem *LS){
	if(!LS->isEmpty()){
		this->set.push_back(LS);
	}
}

LinearSystemSet::LinearSystemSet(vector<LinearSystem*> set){

	for(int i=0; i<(signed)set.size(); i++){
		if(!set[i]->isEmpty()){
			this->set.push_back(set[i]);
		}
	}
}

vector<LinearSystem*> LinearSystemSet::getSet(){
	return this->set;
}

// Add a linear system to the set
void LinearSystemSet::add(LinearSystem *LS){
	if(!LS->isEmpty()){
		this->set.push_back(LS);
	}
}

// Intersect two sets of linear systems
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

// Union of two sets of linear systems
LinearSystemSet* LinearSystemSet::unionWith(LinearSystemSet *LSset){

	vector<LinearSystem*> uniSet = this->set; 		// new union set
	vector<LinearSystem*> set = LSset->getSet();

	for(int i=0; i<(signed)set.size(); i++){
		uniSet.push_back(set[i]);
	}

	return new LinearSystemSet(uniSet);

}

// Union of two sets of linear systems up to bounded cardinality cardinality
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

// Volume of boxes containing the sets
double LinearSystemSet::boundingVol(){

	double vol = 0;
	for(int i=0; i<this->size(); i++){
		vol = vol + this->at(i)->volBoundingBox();
	}
	return vol;

}

// Return the size of this set, i.e, the number of linear systems
int LinearSystemSet::size() {
	return this->set.size();
}

// Return the i-th linear system
LinearSystem* LinearSystemSet::at(int i) {
	return this->set[i];
}

// Check if the current set is empty
bool LinearSystemSet::isEmpty() {
	return this->set.empty();
}

// Print the linear system
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

