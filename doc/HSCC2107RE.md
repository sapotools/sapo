# Repeatability Evaluation
## Description

This document explains how to reproduce the results presented in
Section 4 of the submission “Sapo: A Tool for the Reachability Computation and Parameter Synthesis of Polynomial Dynamical Systems” by T. Dreossi.

The executable ``./bin/sapo`` reproduces the reachbility analysis of Table 1, the parameter synthesis of Table 2, and generates the scripts ``plotFigure3a.m``, ``plotFigure3b.m``, ``plotFigure4a.m``, and ``plotFigure4b.m`` that can be used to plot the Figures 3a, 3b, 4a, 4b.


There are two ways to reproduce the data:

1. Using the [Virtual Machine](#virtualmachine)
2. [Building Sapo](#buildsapo) from source

## <a name="virtualmachine">Virtual Machine</a>

1. Download the virtual machine <a href="https://www.dropbox.com/sh/4ex9yqc3y0p1618/AACnl43b9knKovYaHVTwlkxVa?dl=0">here</a>
2. Launch the virtual machine using <a href="">Virtual Box</a>
3. If required, login with usersame ``sapo`` and password ``sapo``
4. To reproduce the case studies:
``` sh
cd ~/sapo/bin
./sapo
```

To visualize the figures go to the [Visualize Figures](#visfigs) section.


## <a name="buildsapo">Build Sapo</a>

To compile the source code, the following packages are required:

- C++11-compatible compiler, <a href="https://cmake.org/">cmake</a>, <a href="https://www.gnu.org/software/make/">make</a>, <a href="https://www.freedesktop.org/wiki/Software/pkg-config/">pkg-config</a>
- <a href="http://www.ginac.de/CLN/">CLN</a>, <a href="https://www.gnu.org/software/glpk/">GLPK</a>, <a href="http://www.ginac.de/">GiNaC</a> libraries

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

Once that the required packages are installed, donwload, build and install Sapo:
``` sh
git clone https://github.com/tommasodreossi/sapo
cd sapo
cmake .
make
```

This generates the executable ``./bin/sapo``. To reproduce the
case studies:
``` sh
cd bin
./sapo
```

To visualize the figures go to the [Visualize Figures](#visfigs) section.

## <a name="visfigs">Visualize Figures</a>

The executable ``./bin/sapo`` produces the scripts
``plotFigure3a.m``, ``plotFigure3b.m``, ``plotFigure4a.m``, and``plotFigure4b.m`` that can be used to generate the figures
of the paper. The scripts can be run in both Octave and Matlab
and require the ``plotregion`` package available at
<a href="https://www.mathworks.com/matlabcentral/fileexchange/9261-plot-2d-3d-region">here</a> (MatWorks account required) or <a href="https://www.dropbox.com/sh/4ex9yqc3y0p1618/AACnl43b9knKovYaHVTwlkxVa?dl=0">here</a>.


For instance, from the virtual machine, launch Octave:
``` sh
octave
```
and then from Octave command window, include the ``plotregion`` package and run the scripts. For instance:
``` sh
cd ~/sapo/bin
addpath("~/Downloads/plotregion")
figure(1); plotFigure3a;
figure(2); plotFigure3b;
figure(3); plotFigure4a;
figure(4); plotFigure4b;
```

Notes:

1. The scripts might take some time (around 20 seconds) to plot the reachable sets
2. Depending on Matlab or Octave, the color of the plots might differ from those from the paper
