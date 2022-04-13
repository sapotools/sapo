# Sapo

Sapo is a tool for the formal analysis of discrete-time polynomial dynamical systems.

The problems treated by Sapo are:

- Reachability computation, i.e., the calculation of the set of states reachable by the system from a set of initial conditions

- Parameter synthesis, i.e., the refinement of a set of parameters so that the system satisfies a given specification.
  
  For reachability analysis Sapo produces a flowpipe that over-approximates the set of states reachable by the system from a set of initial conditions.

For parameter synthesis Sapo computes a refinement of the given set of parameters such that the system satisfies a given specification. The specification is formalized as a Signal Temporal Logic (STL) formula.

In both cases, the analysis can be done on bounded time.

Sapo consists in a C++ library, named `libSapo`, and an optional command line application, named [`sapo`](sapo-Application), that is meant to ease Sapo usability.

A web user interface for Sapo is available [here](https://github.com/LucaDorigo/webSapo).



Please, refer to the wiki pages to have more pieces of information about:

* (Handled models and set representation)[https://github.com/dreossi/sapo/wiki]
* (Downloading, Building, and Installing Sapo)[https://github.com/dreossi/sapo/wiki/Building-Sapo]
* (the `sapo` command line application)[https://github.com/dreossi/sapo/wiki/sapo-Application]
