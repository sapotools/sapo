/**
 * @file Flowpipe.cpp
 * Represent and manipulate flowpipes of bundle
 * Used to represent the reachable set of a dynamical system
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#include "Flowpipe.h"
#include <string>

/**
 * Constructor that instantiates Flowpipe
 */
Flowpipe::Flowpipe() {
}

/**
 * Constructor that instantiates Flowpipe
 *
 * @param[in] flowpipe vector of bundles
 */
Flowpipe::Flowpipe(vector< Bundle* > flowpipe){
	this->flowpipe = flowpipe;
}

/**
 * Return the i-th bundle
 *
 * @param[in] i index
 * @return i-th bundle
 */
Bundle* Flowpipe::get(int i){
	if(( 0<= i ) && (i < (signed)this->size())){
				return this->flowpipe[i];
		}
		cout<<"Flowpipe::get : i must be between 0 and the flowpipe size";
		exit (EXIT_FAILURE);
}

/**
 * Append a bundle to the flowpipe
 *
 * @param[in] bundle bundle to append
 */
void Flowpipe::append( Bundle* bundle ){
	this->flowpipe.push_back(bundle);
}

/**
 * Display the flowpipe
 */
void Flowpipe::print(){
	for(int i=0; i<this->size(); i++){
		this->flowpipe[i]->getBundle()->print();
	}
}

/**
 * Display the flowpipe in Matlab formal
 */
void Flowpipe::plotRegion(){
	for(int i=0; i<this->size(); i++){
		this->flowpipe[i]->getBundle()->plotRegion();
	}
}

/**
 * Print the linear system in Matlab format (for plotregion script) into a file
 *
 * @param[in] file_name name of the file
 * @param[in] color color of the polytope to plot
 */
void Flowpipe::plotRegionToFile(char *file_name, char color){
	for(int i=0; i<this->size(); i++){
		this->flowpipe[i]->getBundle()->plotRegionToFile(file_name,color);
	}
}

Flowpipe::~Flowpipe() {
	// TODO Auto-generated destructor stub
}

