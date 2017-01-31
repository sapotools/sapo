/**
 * @file main.cpp
 * Demo: Reachability analysis demos
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#include <stdio.h>
#include <iostream>

#include "Common.h"
#include "Bundle.h"
#include "Sapo.h"

#include "VanDerPol.h"
#include "Rossler.h"
#include "SIR.h"
#include "LotkaVolterra.h"
#include "Phosphorelay.h"
#include "Quadcopter.h"

#include "SIRp.h"
#include "Influenza.h"
#include "Ebola.h"

using namespace std;

/**
 * Main function
 *
 */
int main(int argc,char** argv){

    if(argc != 2){
      cerr<<"One argument expected (demo model ID)\n";
      exit(EXIT_FAILURE);
    }

    // Init model
    int reach_steps;
    Model *model;

    switch(atoi(argv[1])){
      case 0:
        model = new VanDerPol(); reach_steps = 300;
      break;
      case 1:
        model = new Rossler(); reach_steps = 250;
      break;
      case 2:
        model = new SIR(); reach_steps = 300;
      break;
      case 3:
        model = new LotkaVolterra(); reach_steps = 500;
      break;
      case 4:
        model = new Phosphorelay(); reach_steps = 200;
      break;
      case 5:
        model = new Quadcopter(); reach_steps = 300;
      break;
      case 10:
        model = new SIRp();
      break;
      case 11:
        model = new Influenza();
      break;
      case 12:
        model = new Ebola();
      break;
      default:
        cerr<<"Unknown demo ID\n";
        exit(EXIT_FAILURE);
    }


    // Sapo's options
    sapo_opt options;
    options.trans = 1;			 // Set transformation (0=OFO, 1=AFO)
    options.decomp = 0;			  // Template decomposition (0=no, 1=yes)
    //options.alpha = 0.5;		// Weight for bundle size/orthgonal proximity
    options.verbose = false;

    if(atoi(argv[1]) < 10){
      Sapo *sapo = new Sapo(model,options);
      Flowpipe* flowpipe = sapo->reach(model->getReachSet(),reach_steps);	// reachability analysis
    }else{
      Sapo *sapo = new Sapo(model,options);
    	LinearSystemSet *synth_parameter_set = sapo->synthesize(model->getReachSet(),model->getParaSet(),model->getSpec());	// parameter synthesis
    }

    exit(EXIT_SUCCESS);
}
