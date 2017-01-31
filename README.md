# Sapo
## Description
Sapo is a C++ tool for the formal analysis of discrete-time polynomial dynamical systems.

The problems treated by Sapo are:

- Reachability computation, i.e., the calculation of the set of states reachable by the system from a set of initial conditions
- Parameter synthesis, i.e., the refinement of a set of parameters so that the system satisfies a given specification.
For reachability analysis Sapo produces a flowpipe that over-approximates the set of states reachable by the system from a set of initial conditions.

For parameter synthesis Sapo computes a refinement of the given set of parameters such that the system satisfies a given specification. The specification is formalized as a Signal Temporal Logic (STL) formula.

In both cases, the analysis can be done on bounded time.

###Models
The dynamical systems supported by Sapo are discrete-time polynomial dynamical systems, i.e., dynamical systems whose evolutions can be described by difference equations of the form x_{k+1} = f(x_k,p)

Reachability computation can be carried out also on systems without parameters whose dynamics look like x_{k+1} = f(x_k) with f : R^n to R^n polynomial.

###Set representation
The flowpipe representing the reachable set consists in a series of sets. The sets supported by Sapo are:

- Boxes (or hyperrectangles), i.e., n-dimensional rectangles
- Parallelotopes, i.e., n-dimensional parallelograms
- Parallelotopes bundles, i.e., finite sets of parallelotopes whose intersections generate polytopes

The parameter synthesis produces a refined set of parameters represented by:
- Polytopes, i.e., n-dimensional polygon

## Required Packages

- C++11-compatible compiler, <a href="https://cmake.org/">cmake</a>, <a href="https://www.gnu.org/software/make/">make</a>, <a href="https://www.freedesktop.org/wiki/Software/pkg-config/">pkg-config</a>
- <a href="https://www.gnu.org/software/glpk/">GLPK</a>, <a href="http://www.ginac.de/CLN/">CLN</a>, <a href="http://www.ginac.de/">GiNaC</a>

## Build Sapo
``` sh
git https://github.com/tommasodreossi/sapo
cd sapo
cmake .
make
```

## Visualization
2D/3D or projections of higher dimensional reachable and parameter sets computed by Sapo can be visualized using the Matlab package plotregion.
