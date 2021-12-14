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

- a C++14-compatible compiler
- <a href="https://cmake.org/">cmake</a> (version>=3.6)
- <a href="https://www.gnu.org/software/make/">make</a>
- <a href="https://github.com/westes/flex">Flex</a> (version >=2.6.3)
- <a href="https://www.gnu.org/software/bison/manual">Bison</a>
- <a href="https://gmplib.org">GNU Multi-Precision Library</a>
- <a href="https://www.gnu.org/software/glpk/">GLPK</a> (version >=5.0)

### Downloading and Compiling Sapo<a id="compile-multithreading"></a>

Once all the required packages have been installed, download and build Sapo by using the following commands:

```sh
git clone https://github.com/dreossi/sapo
cd sapo
cmake .
make
```

This generates the executable `./bin/sapo` which supports multi-threading even though 
[it only runs one thread by default](#multithreading).

If you prefer to compile plain sequential code, then replace the last two lines of
above instructions by

```sh
cmake -DTHREADED_VERSION:BOOL=FALSE .
make
```

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

A complete list of `sapo` command line options can be obtained by using the `-h` option.


### Multi-threading<a id="multithreading"></a>

Even if `sapo` was [compiled with multi-threading support](#compile-multithreading), it only runs one thread 
by default. In order to increase the number of threads, use the `sapo` command line option `-t`. 

This option can take as parameter either nothing or one natural number greater than 0. By issuing the command
```sh
./bin/sapo -t path/to/file.sil
```
`sapo` uses the maximum number of concurrent threads for the executing achitecture, while the line
```sh
./bin/sapo -t 5 path/to/file.sil
```
makes `sapo` running 5 threads. 
