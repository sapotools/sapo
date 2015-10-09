/*
 * CoeffMap.cpp
 *
 *  Created on: Sep 23, 2015
 *      Author: dreossi
 */

#include "CoeffMap.h"

CoeffMap::CoeffMap() {
	// TODO Auto-generated constructor stub

}

void CoeffMap::insert(vector<int> key, lst element){
	this->keys.push_back(key);
	this->elements.push_back(element);
}

int CoeffMap::find(vector<int> key){
	for(int i=0; i<this->keys.size(); i++){
		if( this->keys[i] == key ){
			return i;
		}
	}
	return -1;
}

lst CoeffMap::get(int i){
	return this->elements[i];
}

bool CoeffMap::isIn(vector<int> key){

	for(int i=0; i<this->keys.size(); i++){
		if( this->keys[i] == key ){
			return true;
		}
	}
	return false;
}

CoeffMap::~CoeffMap() {
	// TODO Auto-generated destructor stub
}

