# Sapo

Sapo is a tool for the formal analysis of discrete-time polynomial dynamical systems.

The problems treated by Sapo are:

- Reachability computation, i.e., the calculation of the set of states reachable by the system from a set of initial conditions
- Parameter synthesis, i.e., the refinement of a set of parameters so that the system satisfies a given specification.
  For reachability analysis Sapo produces a flowpipe that over-approximates the set of states reachable by the system from a set of initial conditions.

For parameter synthesis Sapo computes a refinement of the given set of parameters such that the system satisfies a given specification. The specification is formalized as a Signal Temporal Logic (STL) formula.

In both cases, the analysis can be done on bounded time.

Sapo is provided as both a C++ library, to promote integration by other projects, and a stand-alone tool, to ease its usability.

A web user interface for Sapo is avaiable [here](https://github.com/LucaDorigo/webSapo).

### Models

The dynamical systems supported by Sapo are discrete-time polynomial dynamical systems, i.e., dynamical systems whose evolutions can be described by difference equations of the form <!-- $x_{k+1} = f(x_k,p)$ --> <img style="transform: translateY(0.1em); background: white;" src="svg/kukKgGlU7t.svg"/>.

Reachability computation can be carried out also on systems without parameters whose dynamics look like <!-- $x_{k+1} = f(x_k)$ --> <img style="transform: translateY(0.1em); background: white;" src="svg/O2tXhyFGXU.svg"/>  with <!-- $f : R^n \rightarrow R^n$ --> <img style="transform: translateY(0.1em); background: white;" src="svg/6nV9NfgdJZ.svg"/> polynomial.

### Set representation

The flowpipe representing the reachable set consists in a series of sets. The sets supported by Sapo are:

- Boxes (or hyperrectangles), i.e., n-dimensional rectangles
- Parallelotopes, i.e., n-dimensional parallelograms
- Parallelotopes bundles, i.e., finite sets of parallelotopes whose intersections generate polytopes

The parameter synthesis produces a refined set of parameters represented by:

- Polytopes, i.e., n-dimensional polygon

## <a name="buildsapo">Build Sapo</a>

To compile the source code, the following packages are required:

- C++11-compatible compiler, <a href="https://cmake.org/">cmake</a> (version>=3.6), <a href="https://www.gnu.org/software/make/">make</a>, <a href="https://www.freedesktop.org/wiki/Software/pkg-config/">pkg-config</a>
- <a href="http://www.ginac.de/CLN/">CLN</a> (version>=1.3.6), <a href="http://www.ginac.de/">GiNaC</a> (version>=1.7.8), <a href="https://www.gnu.org/software/glpk/">GLPK</a> (version >=5.0) libraries
- <a href="https://github.com/westes/flex">Flex</a> (version >=2.6.3) and <a href="https://www.gnu.org/software/bison/manual">Bison</a>

### Install CLN

Download latest <a href="http://www.ginac.de/CLN/">CLN</a> and install:

```sh
curl http://www.ginac.de/CLN/cln-1.3.6.tar.bz2 | tar -xj
cd cln-1.3.6/
./configure
make
make check
sudo make install
```

### Install GiNaC

Download latest <a href="http://www.ginac.de/">GiNaC</a> and install:

```sh
curl http://www.ginac.de/ginac-1.8.0.tar.bz2 | tar -xj
cd ginac-1.8.0/
./configure
make
make check
sudo make install
```

### Install GLPK

Download latest <a href="https://www.gnu.org/software/glpk/">GLPK</a> and install:

```sh
curl http://ftp.gnu.org/gnu/glpk/glpk-4.61.tar.gz | tar -xz
cd glpk-4.61/
./configure
make
make check
sudo make install
```

### Install Flex

Download latest <a href="https://github.com/westes/flex">Flex</a> and install:

```sh
git clone https://github.com/westes/flex
cd flex/
./autogen.sh
./configure
make
make install
```

### Install Bison
Download latest <a href="https://www.gnu.org/software/bison/manual">Bison</a> and install:

```sh
curl http://ftp.gnu.org/gnu/bison/bison-3.5.1.tar.gz | tar -xz
cd bison-3.5.1/
./configure
make
make check
sudo make install
```

### Install Sapo

Once that the required packages are installed, download, build and install Sapo:

```sh
git clone https://github.com/LucaDorigo/sapo
cd sapo
cmake .
make
```

This generates the executable `./bin/sapo`.

## Using Sapo stand-alone application
The Sapo stand-alone application `sapo` expects as input a SIL file.
The SIL language is defined in the [SIL.md](SIL.md) file in this repository.

`sapo` can be invoked passing the path of the input file

```sh
./bin/sapo path/to/file.sil
```

or without arguments. In that case, `sapo` reads its input from standard input.

```sh
./bin/sapo
```

The outputs are written on standard output.
