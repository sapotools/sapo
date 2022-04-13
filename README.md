# Sapo

Sapo is a tool for the formal analysis of discrete-time polynomial dynamical systems.

The problems treated by Sapo are:

- Reachability computation, i.e., the calculation of the set of states reachable by the system from a set of initial conditions
- Parameter synthesis, i.e., the refinement of a set of parameters so that the system satisfies a given specification.
  For reachability analysis Sapo produces a flowpipe that over-approximates the set of states reachable by the system from a set of initial conditions.

For parameter synthesis Sapo computes a refinement of the given set of parameters such that the system satisfies a given specification. The specification is formalized as a Signal Temporal Logic (STL) formula.

In both cases, the analysis can be done on bounded time.

Sapo consists in a C++ library, named `libSapo`, and an optional standalone application, named `sapo`, that is meant to ease Sapo usability.

A web user interface for Sapo is available [here](https://github.com/LucaDorigo/webSapo).

### Models

The dynamical systems supported by Sapo are discrete-time polynomial dynamical systems, i.e., dynamical systems whose evolutions can be described by difference equations of the form <!-- $x_{k+1} = f(x_k,p)$ --> <img style="transform: translateY(0.1em); background: white;" src="https://github.com/dreossi/sapo/blob/assets/svg/kukKgGlU7t.svg"/>.

Reachability computation can be carried out also on systems without parameters whose dynamics look like <!-- $x_{k+1} = f(x_k)$ --> <img style="transform: translateY(0.1em); background: white;" src="https://github.com/dreossi/sapo/blob/assets/svg/O2tXhyFGXU.svg"/>  with <!-- $f : R^n \rightarrow R^n$ --> <img style="transform: translateY(0.1em); background: white;" src="https://github.com/dreossi/sapo/blob/assets/svg/6nV9NfgdJZ.svg"/> polynomial.

### Set representation

The flowpipe representing the reachable set consists in a series of sets. 
The sets supported by Sapo are:

- Boxes (or hyperrectangles), i.e., n-dimensional rectangles
- Parallelotopes, i.e., n-dimensional parallelograms
- Parallelotopes bundles, i.e., finite sets of parallelotopes whose intersections generate polytopes

The parameter synthesis produces a refined set of parameters represented by:

- Polytopes, i.e., n-dimensional polygon

## Building Sapo

<a id="requirements"></a>
In order to compile *libSapo*, the following packages are mandatory:

- a C++14-compatible compiler
- <a href="https://cmake.org/">cmake</a> (version >= 3.6)
- <a href="https://www.gnu.org/software/make/">make</a>
- <a href="https://gmplib.org">GNU Multi-Precision Library</a>
- <a href="https://www.gnu.org/software/glpk/">GLPK</a> (version >= 5.0)
 
In addition to the above-mentioned packages, compiling the standalone application *sapo* requires the following softwares:

- <a href="https://github.com/westes/flex">Flex</a> (version >= 2.6.3)
- <a href="https://www.gnu.org/software/bison/manual">Bison</a> (version >= 3.6)

### Compiling `libSapo` and `sapo`<a id="compile-multithreading"></a>

The Sapo source code can be obtained by cloning its GitHub repository:
```shell
git clone https://github.com/dreossi/sapo
```
At the end of the command execution the current directory will include a subdirectory, named `sapo`, containing the Sapo source code.

Once all [the required packages](#requirements) have been installed, build both `libSapo` and `sapo` 
by using the following commands:

```sh
cd sapo
cmake .
make
```

This command generates both the dynamic and the static versions of `libSapo` and, 
whenever the requirements are met, produced the executable `./sapo/bin/sapo`. 
The standalone application supports 
multi-threading even though [it only uses one thread by default](#multithreading).

If you prefer to compile plain sequential code, then replace the last two lines of
above instructions by

```sh
cmake -DTHREADED_VERSION:BOOL=FALSE .
make
```

### Testing Sapo

In order to test whether Sapo is properly working in your system, please call from 
the command line:

```
make test
```

### Installing Sapo

In order to install `libSapo` and, possibly, `sapo` system-wise use the command:
```sh
make install
```
The default install directory is `/usr/local/`. The standalone application will be 
placed in `<INSTALL_DIRECTORY}>/bin`, the library in `<INSTALL_DIRECTORY>/lib`, and its header files 
in `<INSTALL_DIRECTORY>/include`. 

To change the default install directory execute:
```sh
cmake -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIRECTORY> .
make install
```
where `<INSTALL_DIRECTORY>` must be replaced by the aimed installed directory, e.g., `/usr` or `./`.

## Using `sapo`
The Sapo standalone application `sapo` accepts input in the
[SIL file format](sapo/doc/SIL.md) from either the standard input or a file.

If the tool has been installed in a directory in the user path, 
the `sapo` application can be executed simply by typing its name on a shell. 
Otherwise, the name of the directory containing the executable must be 
prepended to the executable name itself (e.g., `./sapo/bin/sapo`). In the next examples, 
we will assume that `sapo` has been installed in a directory included among 
those in the user path.

`sapo` reads the problem specification from an input file as in 
```sh
sapo path/to/file.sil
```
or from the standard input 
```sh
cat path/to/file.sil | sapo
```

The outputs are always written on standard output.

Some examples of SIL files are provided for your convenience in the directory [examples](sapo/examples) in this repository.

A complete list of `sapo` command line options can be obtained by using the `-h` option.


### Multi-threading<a id="multithreading"></a>

Even if `sapo` was [compiled with multi-threading support](#compile-multithreading), it only runs one thread 
by default. In order to increase the number of threads, use the `sapo` command line option `-t`. 

This option can take as parameter either nothing or one natural number greater than 0. By issuing the command
```sh
sapo -t path/to/file.sil
```
`sapo` uses the maximum number of concurrent threads for the executing architecture, while the line
```sh
sapo -t 5 path/to/file.sil
```
makes `sapo` running 5 threads. 
