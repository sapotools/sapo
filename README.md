# Sapo

Sapo is a tool for the formal analysis of discrete-time polynomial dynamical systems.

The problems treated by Sapo are:

- Reachability computation, i.e., the calculation of the set of states reachable by the system from a set of initial conditions
- Parameter synthesis, i.e., the refinement of a set of parameters so that the system satisfies a given specification.
  For reachability analysis Sapo produces a flowpipe that over-approximates the set of states reachable by the system from a set of initial conditions.

For parameter synthesis Sapo computes a refinement of the given set of parameters such that the system satisfies a given specification. The specification is formalized as a Signal Temporal Logic (STL) formula.

In both cases, the analysis can be done on bounded time.

Sapo is provided as both a C++ library, to promote integration by other projects, and a stand-alone tool, to ease its usability.

A web user interface for Sapo is available [here](https://github.com/LucaDorigo/webSapo).

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

## <a name="buildsapo">Building Sapo</a>

In order to compile the source code, the following packages are required:

- a C++11-compatible compiler
- <a href="https://cmake.org/">cmake</a> (version>=3.6)
- <a href="https://www.gnu.org/software/make/">make</a>
- <a href="https://www.freedesktop.org/wiki/Software/pkg-config/">pkg-config</a>
- <a href="https://github.com/westes/flex">Flex</a> (version >=2.6.3)
- <a href="https://www.gnu.org/software/bison/manual">Bison</a>
- <a href="http://www.ginac.de/CLN/">CLN</a> (version>=1.3.6)
- <a href="http://www.ginac.de/">GiNaC</a> (version>=1.7.8)
- <a href="https://www.gnu.org/software/glpk/">GLPK</a> (version >=5.0)

### Downloading and Compiling Sapo

Once all the required packages have been installed, download and build Sapo by using the following commands:

```sh
git clone https://github.com/dreossi/sapo
cd sapo
cmake .
make
```

This generates the executable `./bin/sapo`.

## Testing Sapo

In order to test whether Sapo is properly working in your system, please call from 
the command line:

```
make test
```

## Using Sapo stand-alone application
The Sapo stand-alone application `sapo` accepts input in the
[SIL file format](SIL.md) from either the standard input or a file.

Users can execute `sapo` by either passing the path of the input file, as in

```sh
./bin/sapo path/to/file.sil
```

or without any argument. In the latter case, the tool reads its input from 
standard input.

```sh
cat path/to/file.sil | ./bin/sapo
```

The outputs are always written on standard output.

Some examples of SIL files are provided for your convenience in the directory [examples](examples) in this repository.
