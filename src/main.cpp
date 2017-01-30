/**
 * @file demo_sir_reach.cpp
 * Demo: Reachability analysis of SIR epidemic model
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#include <stdio.h>
#include <iostream>

#include "Common.h"
#include "Bundle.h"
#include "Sapo.h"

#include "SIR.h"

using namespace std;

/**
 * Main function
 *
 */
int main(int argc,char** argv){

    // Init model
    Model *model = new SIR();

    ///// SAPO core /////

    // Sapo's options
    sapo_opt options;
    options.trans = 1;			 // Set transformation (0=OFO, 1=AFO)
    options.decomp = 0;			  // Template decomposition (0=no, 1=yes)
    //options.alpha = 0.5;		// Weight for bundle size/orthgonal proximity
    options.verbose = false;

    Sapo *sapo = new Sapo(model,options);
    int reach_steps = 300;
    Flowpipe* flowpipe = sapo->reach(model->getReachSet(),reach_steps);	// reachability analysis

    // Store the constructed flowpipe in file sir.m (in Matlab format)
    char file_name[] = "sir_flowpipe.m";
    flowpipe->plotRegionToFile(file_name,'w');

}
