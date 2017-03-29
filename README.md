# Sapo
## Description
Sapo is a C++ tool for the formal analysis of discrete-time polynomial dynamical systems.

The problems treated by Sapo are:

- Reachability computation, i.e., the calculation of the set of states reachable by the system from a set of initial conditions
- Parameter synthesis, i.e., the refinement of a set of parameters so that the system satisfies a given specification.
For reachability analysis Sapo produces a flowpipe that over-approximates the set of states reachable by the system from a set of initial conditions.

For parameter synthesis Sapo computes a refinement of the given set of parameters such that the system satisfies a given specification. The specification is formalized as a Signal Temporal Logic (STL) formula.

In both cases, the analysis can be done on bounded time.

### Models
The dynamical systems supported by Sapo are discrete-time polynomial dynamical systems, i.e., dynamical systems whose evolutions can be described by difference equations of the form x_{k+1} = f(x_k,p)

Reachability computation can be carried out also on systems without parameters whose dynamics look like x_{k+1} = f(x_k) with f : R^n to R^n polynomial.

### Set representation
The flowpipe representing the reachable set consists in a series of sets. The sets supported by Sapo are:

- Boxes (or hyperrectangles), i.e., n-dimensional rectangles
- Parallelotopes, i.e., n-dimensional parallelograms
- Parallelotopes bundles, i.e., finite sets of parallelotopes whose intersections generate polytopes

The parameter synthesis produces a refined set of parameters represented by:
- Polytopes, i.e., n-dimensional polygon

## <a name="buildsapo">Build Sapo</a>

To compile the source code, the following packages are required:

- C++11-compatible compiler, <a href="https://cmake.org/">cmake</a>, <a href="https://www.gnu.org/software/make/">make</a>, <a href="https://www.freedesktop.org/wiki/Software/pkg-config/">pkg-config</a>
- <a href="http://www.ginac.de/CLN/">CLN</a>,  <a href="http://www.ginac.de/">GiNaC</a>, <a href="https://www.gnu.org/software/glpk/">GLPK</a> libraries

### Install CLN

Download latest <a href="http://www.ginac.de/CLN/">CLN</a> and install:
``` sh
curl http://www.ginac.de/CLN/cln-1.3.4.tar.bz2 | tar -xj;
cd cln-1.3.4/;
./configure;
make;
make check;
sudo make install;
```

### Install GiNaC

Download latest <a href="http://www.ginac.de/">GiNaC</a> and install:
``` sh
curl http://www.ginac.de/ginac-1.7.2.tar.bz2 | tar -xj;
cd ginac-1.7.2/;
./configure;
make;
make check;
sudo make install;
```

### Install GLPK

Download latest <a href="https://www.gnu.org/software/glpk/">GLPK</a> and install:
``` sh
curl http://ftp.gnu.org/gnu/glpk/glpk-4.61.tar.gz | tar -xz;
cd glpk-4.61/;
./configure;
make;
make check;
sudo make install;
```

### Install Sapo

Once that the required packages are installed, download, build and install Sapo:
``` sh
git clone https://github.com/tommasodreossi/sapo
cd sapo
cmake .
make
```

This generates the executable ``./bin/sapo``.
To reproduce the case studies:
``` sh
cd bin
./sapo
```

To visualize the figures go to the [Visualize Figures](#visfigs) section.

## <a name="visfigs">Visualize Figures</a>

The executable ``./bin/sapo`` produces the scripts
``plotFigure3a.m``, ``plotFigure3b.m``, ``plotFigure4a.m``, and``plotFigure4b.m`` that can be used to generate the figures
of the paper. The scripts can be run in both Octave and Matlab
and require the ``plotregion`` package available
<a href="https://www.mathworks.com/matlabcentral/fileexchange/9261-plot-2d-3d-region">here</a> (MatWorks account required) or <a href="https://www.dropbox.com/sh/4ex9yqc3y0p1618/AACnl43b9knKovYaHVTwlkxVa?dl=0">here</a>.


For instance, from the virtual machine, launch Octave:
``` sh
octave
```
and then from Octave command window, include the ``plotregion`` package and run the scripts. For instance:
``` sh
cd ~/sapo/bin
addpath("~/Downloads")
figure(1); plotFigure3a;
figure(2); plotFigure3b;
figure(3); plotFigure4a;
figure(4); plotFigure4b;
```
