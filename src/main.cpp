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
      cerr<<"One argument expected (case study ID)\n";
      exit(EXIT_FAILURE);
    }


    int caseID = atoi(argv[1]);

    // Init model
    int reach_steps;
    Model *model;

    switch(caseID){
      case 11:
        model = new VanDerPol(); reach_steps = 300;
      break;
      case 12:
        model = new Rossler(); reach_steps = 250;
      break;
      case 13:
      case 32:
        model = new SIR(false); reach_steps = 300;
      break;
      case 14:
        model = new LotkaVolterra(); reach_steps = 500;
      break;
      case 15:
        model = new Phosphorelay(); reach_steps = 200;
      break;
      case 16:
        model = new Quadcopter(); reach_steps = 300;
      break;
      case 21:
      case 41:
      case 42:
        model = new SIRp(); reach_steps = 300;
      break;
      case 22:
        model = new Influenza();
      break;
      case 23:
        model = new Ebola();
      break;
      case 31:
        model = new SIR(true); reach_steps = 300;
      break;
      default:
        cerr<<"Unknown case study ID\n";
        exit(EXIT_FAILURE);
    }


    // Sapo's options
    sapo_opt options;
    options.trans = 1;			 // Set transformation (0=OFO, 1=AFO)
    options.decomp = 0;			  // Template decomposition (0=no, 1=yes)
    //options.alpha = 0.5;		// Weight for bundle size/orthgonal proximity
    options.verbose = false;

    if(caseID < 20){
      Sapo *sapo = new Sapo(model,options);
      Flowpipe* flowpipe = sapo->reach(model->getReachSet(),reach_steps);	// reachability analysis
      char file_name[] = "chiurlo.m";
      flowpipe->plotRegionToFile(file_name,'w');

      exit(EXIT_SUCCESS);
    }

    if(caseID < 30){
      Sapo *sapo = new Sapo(model,options);
    	LinearSystemSet *synth_parameter_set = sapo->synthesize(model->getReachSet(),model->getParaSet(),model->getSpec());	// parameter synthesis
      exit(EXIT_SUCCESS);
    }

    if(caseID < 40){ // Figures 3a 3b
      Sapo *sapo = new Sapo(model,options);
      Flowpipe* flowpipe = sapo->reach(model->getReachSet(),reach_steps);	// reachability analysis

      //	// Store the constructed flowpipe in file sir.m (in Matlab format)
      char file_name[] = "plotFigure.m";
      flowpipe->plotRegionToFile(file_name,'w');

      // Set picture appearence
      ofstream matlab_script;
    	matlab_script.open (file_name, ios_base::app);
      matlab_script<<"xlabel('s');\n";
      matlab_script<<"ylabel('i');\n";
      matlab_script<<"zlabel('r');\n";
      matlab_script<<"axis([0 1 0 0.7 0 0.8]);\n";
      matlab_script<<"view([74 23]);\n";
      matlab_script<<"grid on;";
      matlab_script.close();

      exit(EXIT_SUCCESS);
    }

    if(caseID < 50){ // Figures 4a 4b
      Sapo *sapo = new Sapo(model,options);
    	LinearSystemSet *synth_parameter_set = sapo->synthesize(model->getReachSet(),model->getParaSet(),model->getSpec());	// parameter synthesis

      // Store the first synthesized parameter set in file sir_synth.m (in Matlab format)
      char file_name[] = "plotFigure.m";

      if(caseID == 41){
        model->getParaSet()->at(0)->plotRegionToFile(file_name,'w');
        synth_parameter_set->at(0)->plotRegionToFile(file_name,'k');
        // Set picture appearence
        ofstream matlab_script;
    	   matlab_script.open (file_name, ios_base::app);
         matlab_script<<"xlabel('\\beta');\n";
         matlab_script<<"ylabel('\\gamma');\n";
         matlab_script<<"axis([0.17 0.21 0.045 0.065]);\n";
         matlab_script<<"grid on;";
         matlab_script.close();
         exit(EXIT_SUCCESS);
       }

       if(caseID == 42){
         //Rechability computation under the synthesized parameter set
         Flowpipe* flowpipe = sapo->reach(model->getReachSet(),synth_parameter_set->at(0),reach_steps);
         flowpipe->plotRegionToFile(file_name,'w');
         // Set picture appearence
         ofstream matlab_script;
       	 matlab_script.open (file_name, ios_base::app);
         matlab_script<<"xlabel('s');\n";
         matlab_script<<"ylabel('i');\n";
         matlab_script<<"zlabel('r');\n";
         matlab_script<<"axis([0 1 0 0.7 0 0.8]);\n";
         matlab_script<<"view([120 40]);\n";
         matlab_script<<"grid on;";
         matlab_script.close();

         exit(EXIT_SUCCESS);
       }
    }


    exit(EXIT_SUCCESS);
}
