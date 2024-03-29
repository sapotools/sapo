# Sapo

Sapo is a tool for the formal analysis of discrete-time polynomial dynamical systems.

Sapo can:

- compute the reachable set, i.e., the set of states reachable by the system from a set of initial conditions

- synthesize a parameter set satisfying a specification.

For reachability analysis, Sapo produces a flowpipe that over-approximates the set of states reachable by the system from a set of initial conditions.

For parameter synthesis, Sapo computes a refinement of the given set of parameters such that the system satisfies a given specification. The specification is formalized as a Signal Temporal Logic (STL) formula.

In both cases, the analysis can be done on bounded time.

Sapo consists in a C++ library, named `libSapo`, and an optional command line application, named `sapo`, that is meant to ease Sapo usability.

A web user interface for Sapo is available [here](https://github.com/LucaDorigo/webSapo).



Please, refer to the wiki pages to have more pieces of information about:
* [handled models and set representation](https://github.com/dreossi/sapo/wiki)
* [downloading, building, and installing Sapo](https://github.com/dreossi/sapo/wiki/Building-Sapo)
* [the `sapo` command line application](https://github.com/dreossi/sapo/wiki/sapo-Application)
