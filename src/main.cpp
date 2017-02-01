/**
 * @file main.cpp
 * main: This main file reproduces the experiments reported in "Sapo: Reachability Computation and Parameter Synthesis of Polynomial Dynamical Systems"
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

int main(int argc,char** argv){

  // Sapo's options
  sapo_opt options;
  options.trans = 1;			 // Set transformation (0=OFO, 1=AFO)
  options.decomp = 0;			  // Template decomposition (0=no, 1=yes)
  //options.alpha = 0.5;		// Weight for bundle size/orthgonal proximity
  options.verbose = false;


  cout<<"TABLE 1"<<endl;
  // Load modles
  vector< Model* > reach_models;
  vector< int > reach_steps;
  reach_models.push_back(new VanDerPol());      reach_steps.push_back(300);
  reach_models.push_back(new Rossler());        reach_steps.push_back(250);
  reach_models.push_back(new SIR(false));       reach_steps.push_back(300);
  reach_models.push_back(new LotkaVolterra());  reach_steps.push_back(500);
  reach_models.push_back(new Phosphorelay());   reach_steps.push_back(200);
  reach_models.push_back(new Quadcopter());     reach_steps.push_back(300);

  // Compute reach sets
  for(int i=0; i<reach_models.size(); i++){

    cout<<"Model: "<<reach_models[i]->getName()<<"\tReach steps: "<<reach_steps[i]<<"\t";

    Sapo *sapo = new Sapo(reach_models[i],options);
    Flowpipe* flowpipe = sapo->reach(reach_models[i]->getReachSet(),reach_steps[i]);	// reachability analysis
  }
  cout<<"\n";


  cout<<"TABLE 2"<<endl;
  // Load models
  vector< Model* > synth_models;
  synth_models.push_back(new SIRp());
  synth_models.push_back(new Influenza());
  synth_models.push_back(new Ebola());

  // Synthesize parameters
  for(int i=0; i<synth_models.size(); i++){

    cout<<"Model: "<<synth_models[i]->getName()<<"\t";

    Sapo *sapo = new Sapo(synth_models[i],options);
    LinearSystemSet *synth_parameter_set = sapo->synthesize(synth_models[i]->getReachSet(),synth_models[i]->getParaSet(),synth_models[i]->getSpec());	// parameter synthesis
  }
  cout<<"\n";


  cout<<"FIGURE 3a"<<endl;
  Model* sir3a = new SIR(true);
  Sapo *sapo3a = new Sapo(sir3a,options);

  // Compute reach set with box template
  cout<<"Model: "<<sir3a->getName()<<"\tReach steps: 300\t";
  Flowpipe* flowpipe3a = sapo3a->reach(sir3a->getReachSet(),300);

  // Generate matlab script to plot flowpipe
  char fig3a[] = "plotFigure3a.m";
  flowpipe3a->plotRegionToFile(fig3a,'w');
  // Set picture appearence
  ofstream matlab_script;
  matlab_script.open (fig3a, ios_base::app);
  matlab_script<<"xlabel('s');\n";
  matlab_script<<"ylabel('i');\n";
  matlab_script<<"zlabel('r');\n";
  matlab_script<<"axis([0 1 0 0.7 0 0.8]);\n";
  matlab_script<<"view([74 23]);\n";
  matlab_script<<"grid on;";
  matlab_script.close();
  cout<<fig3a<<" generated\n"<<endl;


  cout<<"FIGURE 3b"<<endl;
  Model* sir3b = new SIR(false);
  Sapo *sapo3b = new Sapo(sir3b,options);

  // Compute reach set with bundle template
  cout<<"Model: "<<sir3b->getName()<<"\tReach steps: 300\t";
  Flowpipe* flowpipe3b = sapo3b->reach(sir3b->getReachSet(),300);

  // Generate matlab script to plot flowpipe
  char fig3b[] = "plotFigure3b.m";
  flowpipe3b->plotRegionToFile(fig3b,'w');
  // Set picture appearence
  matlab_script.open (fig3b, ios_base::app);
  matlab_script<<"xlabel('s');\n";
  matlab_script<<"ylabel('i');\n";
  matlab_script<<"zlabel('r');\n";
  matlab_script<<"axis([0 1 0 0.7 0 0.8]);\n";
  matlab_script<<"view([74 23]);\n";
  matlab_script<<"grid on;";
  matlab_script.close();
  cout<<fig3b<<" generated\n"<<endl;


  cout<<"FIGURE 4a"<<endl;
  Model *sirp = new SIRp();
  Sapo *sapo_synth = new Sapo(sirp,options);

  // Synthesize parameters
  cout<<"Model: "<<sirp->getName()<<"\t";
  LinearSystemSet *synth_parameter_set = sapo_synth->synthesize(sirp->getReachSet(),sirp->getParaSet(),sirp->getSpec());

  // Generate matlab script to plot the parameter set
  char fig4a[] = "plotFigure4a.m";
  sirp->getParaSet()->at(0)->plotRegionToFile(fig4a,'w');
  synth_parameter_set->at(0)->plotRegionToFile(fig4a,'k');
  // Set picture appearence
  matlab_script.open (fig4a, ios_base::app);
  matlab_script<<"xlabel('\\beta');\n";
  matlab_script<<"ylabel('\\gamma');\n";
  matlab_script<<"axis([0.17 0.21 0.045 0.065]);\n";
  matlab_script<<"grid on;";
  matlab_script.close();
  cout<<fig4a<<" generated\n"<<endl;


  cout<<"FIGURE 4b"<<endl;
  cout<<"Model: "<<sirp->getName()<<"\t";

  // Compute reach set with contrained parameters
  Sapo *sapo_reach = new Sapo(sirp,options);
  Flowpipe* flowpipe = sapo_reach->reach(sirp->getReachSet(),synth_parameter_set->at(0),300);

  // Generate matlab script to plot the reach set
  char fig4b[] = "plotFigure4b.m";
  flowpipe->plotRegionToFile(fig4b,'w');
  // Set picture appearence
  matlab_script.open (fig4b, ios_base::app);
  matlab_script<<"xlabel('s');\n";
  matlab_script<<"ylabel('i');\n";
  matlab_script<<"zlabel('r');\n";
  matlab_script<<"axis([0 1 0 1 0 1]);\n";
  matlab_script<<"view([126 35]);\n";
  matlab_script<<"grid on;";
  matlab_script.close();
  cout<<fig4b<<" generated"<<endl;

  exit(EXIT_SUCCESS);
}
