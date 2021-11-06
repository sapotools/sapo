/**
 * @file main.cpp
 * main: This main file reproduces the experiments reported in "Sapo:
 * Reachability Computation and Parameter Synthesis of Polynomial Dynamical
 * Systems"
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#include <fstream>
#include <iostream>
#include <stdio.h>

#include "AutoGenerated.h"
#include "Bundle.h"
#include "Common.h"
#include "Sapo.h"
#include "driver.h"

using namespace std;

Sapo init_sapo(Model *model, const AbsSyn::InputData& data, const bool verbose=false)
{
  Sapo sapo(model);

  sapo.trans = data.getTransValue();
  sapo.decomp = data.getDecomposition() ? 1 : 0;
  sapo.alpha = data.getAlpha();
  sapo.time_horizon = data.getIterations();
  sapo.max_param_splits = data.getMaxParameterSplits();
  sapo.verbose = verbose;

  return sapo;
}

void reach_analysis(Sapo &sapo, const Model *model)
{
    Flowpipe flowpipe;
    // if the model does not specify any parameter set
    if (model->getParams().nops() == 0) {

      // perform the reachability analysis
      flowpipe = sapo.reach(*(model->getReachSet()), sapo.time_horizon);
    } else {

      // perform the parametric reachability analysis
      flowpipe = sapo.reach(*(model->getReachSet()), *(model->getParaSet()), sapo.time_horizon);
    }

    cout << "FLOWPIPES" << endl
         << "=================" << endl
         << flowpipe << endl;
}

void apply_parameters_and_print(std::ostream& os, Sapo& sapo,
                                const std::list<LinearSystemSet>& synth_params, 
                                std::function<void(Sapo &, const LinearSystemSet&)> funct)
{
  if (every_set_is_empty(synth_params)) {
    os << "=================" << endl
          << "----empty set----" << endl
          << endl;

    return;
  } 

  for (auto p_it = std::cbegin(synth_params); p_it!= std::cend(synth_params); ++p_it) {
    if (p_it->size() != 0) {
      os << "=================" << endl;
      funct(sapo, *p_it);
    }
  }
}

void synthesis(Sapo &sapo, const Model *model)
{
    // Synthesize parameters
  std::list<LinearSystemSet> synth_params
      = sapo.synthesize(*(model->getReachSet()), *(model->getParaSet()),
                        model->getSpec(), sapo.max_param_splits);

  cout << "FLOWPIPES" << endl;
  auto reachability = [model](Sapo& sapo, const LinearSystemSet& pSet) {
    cout << sapo.reach(*(model->getReachSet()), pSet,
                        sapo.time_horizon) << endl;
  };
  apply_parameters_and_print(cout, sapo, synth_params, reachability);

  cout << "PARAMETER SETS" << endl;

  auto identity = [](Sapo& sapo, const LinearSystemSet& pSet) {
    (void)sapo; // this is just to avoid the warning

    cout << pSet << endl << endl;
  };
  apply_parameters_and_print(cout, sapo, synth_params, identity); 
}

int main(int argc, char **argv)
{

  driver drv;
  string file;

  if (argc == 1) {
    file = "-"; // causes lexer to use stdin
  } else if (argc == 2) {
    file = argv[1]; // provided file used
  } else {
    cerr << "usage: " << argv[0] << " <filename>" << endl;
    cerr << "       " << argv[0] << endl;
    return 1;
  }

  if (drv.parse(file) != 0)
    return 2;

  Model *model = new AutoGenerated(drv.data);
  Sapo sapo = init_sapo(model, drv.data);

  switch (drv.data.getProblem()){
    case AbsSyn::problemType::REACH:
      reach_analysis(sapo, model);
      break;
    case AbsSyn::problemType::SYNTH:
      synthesis(sapo, model);
      break;
    default:
      std::cerr << "Unsupported problem type" << std::endl;
      exit(EXIT_FAILURE);
  }

  delete model;

  exit(EXIT_SUCCESS);
}
