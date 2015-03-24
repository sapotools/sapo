/*
 * LinearSystemSet.cpp
 *
 *  Created on: Nov 17, 2014
 *      Author: dreossi
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

	for(int i=0; i<set.size(); i++){
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

	for(int i=0; i<this->set.size(); i++){
		for(int j=0; j<set.size(); j++){
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

