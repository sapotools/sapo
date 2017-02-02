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

## Build Sapo
To compile the source code, the following packages are required:

- C++11-compatible compiler, <a href="https://cmake.org/">cmake</a>, <a href="https://www.gnu.org/software/make/">make</a>, <a href="https://www.freedesktop.org/wiki/Software/pkg-config/">pkg-config</a>
- <a href="http://www.ginac.de/CLN/">CLN</a>,  <a href="http://www.ginac.de/">GiNaC</a>, <a href="https://www.gnu.org/software/glpk/">GLPK</a> libraries

### Install CLN

1. Download latest <a href="http://www.ginac.de/CLN/">CLN</a> and untar
2. In CLN folder:
``` sh
./configure
make
make check
sudo make install
```

### Install GiNaC

1. Download latest <a href="http://www.ginac.de/">GiNaC</a> and untar
2. In GiNaC folder:
``` sh
./configure
make
make check
sudo make install
```

### Install GLPK

1. Download latest <a href="https://www.gnu.org/software/glpk/">GLPK</a> and untar
2. In glpk folder:
``` sh
./configure
make
make check
sudo make install
```

### Install Sapo

Once that the required packages are installed, download, build and install Sapo:
``` sh
git clone https://github.com/tommasodreossi/sapo
cd sapo
cmake .
make
```

