/**
 * @file Flowpipe.cpp
 * Represent and manipulate flowpipes of bundle
 * Used to represent the reachable set of a dynamical system
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#include <string>
#include <fstream>

#include "Flowpipe.h"

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
const Bundle* Flowpipe::get(const unsigned int i) const{
	if(( 0<= i ) && (i < this->size())){
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
void Flowpipe::print() const {
	for (auto it=std::begin(flowpipe); it!=std::end(flowpipe); ++it) {
		std::cout << (*it)->getLinearSystem() << std::endl << std::endl;
	}
}

/**
 * Display the flowpipe in Matlab formal
 */
void Flowpipe::plotRegion() const{
	for (auto it=std::begin(flowpipe); it!=std::end(flowpipe); ++it) {
		(*it)->getLinearSystem().plotRegion();
	}
}

/**
 * Print the linear system in Matlab format (for plotregion script) into a file
 *
 * @param[in] file_name name of the file
 * @param[in] color color of the polytope to plot
 */
void Flowpipe::plotRegionToFile(const char *file_name, const char color) const {
	std::ofstream outstream(file_name, ios_base::app);
	for (auto it=std::begin(flowpipe); it!=std::end(flowpipe); ++it) {
		(*it)->getLinearSystem().plotRegion(outstream, color);
	}
}

/**
 * Print the projection of the variable in time using Matlab format into a file
 *
 * @param[in] var variable to be plotted
 * @param[in] time_step time step for time projection
 * @param[in] file_name name of the file
 * @param[in] color color of the polytope to plot
 */
void Flowpipe::plotProjToFile(const unsigned int var, const double time_step, const char *file_name, const char color) const {

	if( var < 0 || var >= this->flowpipe[0]->getDim() ){
		cout<<"Flowpipe::plotProjToFile : i must be between 0 and the system dimension";
		exit (EXIT_FAILURE);
	}

	ofstream matlab_script(file_name, ios_base::app);

	// select figure
	matlab_script<<"figure("<<var+1<<")\n";

	// print time
	matlab_script<<"t = [ ";
	for(unsigned int i=0; i<this->size(); i++){
		matlab_script<<i*time_step<<" ";
	}
	matlab_script<<" ];\n";

	// print lower offsets
	matlab_script<<"varm = [ ";
	for(auto it=std::begin(this->flowpipe); it!=std::end(this->flowpipe); ++it){
		matlab_script<<-(*it)->getOffm(var)<<" ";
	}
	matlab_script<<" ];\n";

	// print upper offsets
	matlab_script<<"varp = [ ";
	for(auto it=std::begin(this->flowpipe); it!=std::end(this->flowpipe); ++it){
		matlab_script<<(*it)->getOffp(var)<<" ";
	}
	matlab_script<<" ];\n";

	matlab_script<<"T = [t,fliplr(t)];\n";
	matlab_script<<"X = [varm,fliplr(varp)];\n";
	matlab_script<<"fill(T,X,'"<<color<<"');\n";
	matlab_script.close();

}

Flowpipe::~Flowpipe() {
	// TODO Auto-generated destructor stub
}

