# Xyce Configure, Build and Installation Guide Using CMake

This guide describes the basic process for compiling and installing a Xyce
binary using the CMake build system (it is easiest to view these instructions
with full formatting on the [Xyce GitHub
website](https://github.com/Xyce/Xyce/blob/master/INSTALL)). For instructions
on building Xyce with the autotools system, see the [Xyce Building
Guide](https://xyce.sandia.gov/documentation/BuildingGuide.html). These
directions are appropriate for Unix-like systems (like Linux, BSD, MacOS or
Cygwin, e.g.). For building natively on Windows, the primary differences
involve syntax (DOS vs. \*nix), and the lack of a Windows package manager.

If you do not want to build from source, binary installers for Windows, Mac and
Red Hat Linux are made available for every release of Xyce on the
[Xyce website](https://xyce.sandia.gov).

Xyce can be built in two variants: serial or with MPI parallelism. The variant
is determined by whether Trilinos is built in serial or with MPI. While the
parallel and serial variants can co-exist on a system, they must be installed
in different directories (for both Trilinos and Xyce).

If you have not yet done so, clone the [Xyce
repository](https://github.com/Xyce/Xyce) from GitHub, or download one of the
release tarballs.

Much of the difficulty in building Xyce involves properly obtaining and/or
building the third-party libraries (TPLs), particularly
[Trilinos](https://trilinos.github.io/). We are in the process of developing a
CMake "superbuild" capability that can automatically download and build many of
the TPLs. The [Using the Superbuild](#using-the-superbuild) section covers that
approach. The [Standard Build Approach](#standard-build-approach) section
covers the more traditional method. For either method, a certain minimal set of
dependencies are required, which is covered first.

If many flags are being given to CMake on the command line (for Trilinos or
Xyce), you might want to create a configuration script. See the [Configuration
Scripts](#configuration-scripts) section for more detail.

## Obtaining the Dependencies

Most of the software needed to build Xyce is available via package management
systems. On MacOS we recommend using [MacPorts](https://www.macports.org/) or
[Homebrew](https://brew.sh/). You will need to obtain the following tools:
- C++11-compatible compiler suite (e.g., gcc 4.9 or later, clang 3.3 or later)
- Fortran compiler (e.g., gfortran)\
  (This is for Trilinos, and is technically optional; but without it, the build
  will not be optimized for performance.)
- [CMake](https://cmake.org)
- [Bison](https://www.gnu.org/software/bison) (3.0.4 or later)
- [flex](https://github.com/westes/flex) (2.5.34 or later)

On Windows, flex and Bison are available via the
[WinFlexBison](https://github.com/lexxmark/winflexbison) package.

You will also need the following libraries:
- [BLAS](http://www.netlib.org/blas)
- [LAPACK](http://www.netlib.org/lapack)

On MacOS, BLAS and LAPACK can be provided by Xcode or the package managers.

Note that, on some systems, you may also need to install the "development"
versions of certain libraries.

## Using the Superbuild

While easy, this approach has not been thoroughly tested, so should be
considered a "beta" capability. Also, it currently installs everything into an
"install" directory in the build directory, which may not be ideal for most
users. Therefore, we encourage people to use the [standard build
approach](#standard-build-approach), which involves only a few more steps.

Assuming the dependencies have been installed, CMake will automatically build
the following components:
- ADMS
- FFTW (or uses MKL?)
- SuiteSparse
- Trilinos

Then CMake will compile a serial build of Xyce that enables the (optional)
[Xyce/ADMS](https://xyce.sandia.gov/documentation/XyceADMSGuide.html) model
plugin capability. Note that, since CMake is building several packages, the
process could take a long time.

To perform a superbuild, follow this procedure:
- Create a "build" directory somewhere on your system (the location and
  directory name are not important)
- From the build directory, run the following command:
  ```sh
  cmake -D Xyce_USE_SUPERBUILD=ON <path/to/Xyce>
  ```
- Then run
  ```sh
  cmake --build . -j 2
  ```
  The "-j 2" designates the number of processors being used for compiling.
  Choose an appropriate number for your system.

## Standard Build Approach

The standard build approach requires all the third-party software to be in
place prior to invoking CMake to configure Xyce. Be sure you have followed the
instructions in the [Obtaining the Dependencies](#obtaining-the-dependencies)
section, above. Then proceed through the following tasks, which are broken into
three steps:
- Obtaining the TPLs (except Trilinos)
- Building [Trilinos](https://trilinos.github.io/)
- Building Xyce

### Obtaining the TPLs

The minimal third-party libraries needed for Xyce are:
- [SuiteSparse](http://faculty.cse.tamu.edu/davis/suitesparse.html) (2.1.1 or
  later)
- [FFTW](http://www.fftw.org) (3.x)

For parallel builds you will also need [Open MPI](http://www.open-mpi.org)
([MPICH](https://www.mpich.org/) should also work, but the Xyce team does not
regularly build with that library).

The above are available in most package managers. Again, note that, on some
systems, you may also need to install the "development" versions of the
libraries. If the libraries are not available via your package manager, you
will need to download and build them from source.

### Building Trilinos

While [Trilinos](https://trilinos.github.io/) is available in some package
managers, it may not have all the features required by Xyce. The following
process will provide a serial Trilinos installation that will contain only the
libraries needed by Xyce. For the parallel version, be sure to understand the
serial build process prior to reading the [Building Trilinos with MPI
Parallelism](#building-trilinos-with-mpi-parallelism) section, below.

First, download Trilinos version 12.12.1 from the [Trilinos GitHub
Page](https://github.com/trilinos/trilinos), or use this [direct
link](https://github.com/trilinos/Trilinos/archive/refs/tags/trilinos-release-12-12-1.tar.gz).
Be sure to compile version 12.12.1. Other versions will not work with this
process, as they require slightly different build options.

An "initial cache" file is included in the Xyce repository with a typical set
of settings for a Xyce-oriented serial Trilinos build. The file,
`trilinos-config.cmake` is located in the `cmake/trilinos` directory.

Prior to building Trilinos, create a "build" directory in a convenient
location. (The name does not matter, and can be something simple like
`trilinos_build`). On Unix-like systems, the default Trilinos installation
location is `/usr/local`. This can be changed by adding the following flag to
the CMake invocation.
```
-D CMAKE_INSTALL_PREFIX=</path/to/install>
```
As with Xyce, both a parallel and serial build of Trilinos can exist on the
same system, but they must be in different directories. (We recommend
specifying unique sub-directories in `/usr/local`, such as
`/usr/local/trilinos_serial`.) If you plan to build Xyce again, such as when a
new version comes out, you will not need to repeat the Trilinos build/install
process (unless the recommended version of Trilinos is changed, which does not
happen often).

If you have compilers or libraries in non-standard locations, see the [Other
Trilinos Options](#other-trilinos-options) section, below.

To build Trilinos (and install in `/usr/local`), enter the build directory and
run:
```sh
cmake -C path/to/Xyce/cmake/trilinos/trilinos-config.cmake <path/to/Trilinos>
```
Once the configuration step is done, run the following in the build directory:
```sh
cmake --build . -j 2 -t install
```
The "-j 2" designates the number of processors being used for compiling
Trilinos. Choose an appropriate number for your system. After the Xyce
installation is successful, you may delete the Trilinos build directory if you
do not plan on building Trilinos again (see the [Uninstalling
Xyce](uninstalling-xyce) section below first, however).

#### Other Trilinos Options

The following are various flags that might be needed when building Trilinos.
These can be added to the command line, or put in a [configuration
script](#configuration-scripts).

CMake will use the first compiler set it finds on your system. Alternatively,
you can specify the compilers by adding the following flags to the CMake
invocation:
```sh
-D CMAKE_C_COMPILER=<C-compiler> \
-D CMAKE_CXX_COMPILER=<C++-compiler> \
-D CMAKE_Fortran_COMPILER=<Fortran-compiler> \
```
You may need to use a full path, if they are not in your `$PATH`.

Similarly, if the third-party libraries (AMD, BLAS and LAPACK) are not in your
`$PATH`, use the following flags to help CMake find the libraries:
```sh
-D AMD_LIBRARY_DIRS=</path/to/AMD/lib> \
-D AMD_INCLUDE_DIRS=</path/to/AMD/include> \
-D BLAS_LIBRARY_DIRS=</path/to/BLAS/lib> \
-D LAPACK_LIBRARY_DIRS=</path/to/LAPACK/lib> \
```
See the [Trilinos Build
Reference](https://docs.trilinos.org/files/TrilinosBuildReference.html) for
more details on the Trilinos build options.

#### Building Trilinos with MPI Parallelism

To enable MPI parallelism in Xyce, Trilinos must be built with MPI enabled. A
second "initial cache" file, called `trilinos-config-MPI.cmake`, is provided
for MPI builds. In addition, the MPI compilers must be explicitly specified to
CMake. The following CMake invocation should work on most systems:
```sh
cmake\
-C path/to/Xyce/cmake/trilinos/trilinos-config-MPI.cmake \
-D CMAKE_C_COMPILER=mpicc \
-D CMAKE_CXX_COMPILER=mpicxx \
-D CMAKE_Fortran_COMPILER=mpifort \
<path/to/Trilinos>
```
Note that the exact compiler names may be different on your system. Also, the
above invocation will install Trilinos in `/usr/local`.

As with the serial build, once the configuration step is done, run the
following in the build directory to compile and install Trilinos:
```sh
cmake --build . -j 2 -t install
```

### Building Xyce

Once Trilinos is built, you can build and install Xyce. By default, Xyce will
be installed in the `/usr/local/` directory. To specify a different
installation location, add the following flag to the CMake invocation:
```sh
-D CMAKE_INSTALL_PREFIX=</path/to/install>
```
Again, if you plan to have both a parallel and serial build of Xyce on your
system, they must be in different directories. (We recommend specifying unique
sub-directories in `/usr/local`, such as `/usr/local/xyce_serial`.)

If Trilinos is not located in your path, add the following flag to the CMake
invocation:
```sh
-D Trilinos_ROOT=</path/to/Trilinos_install>
```
Create a build directory for Xyce (such as `xyce_build`) and go into that
directory. Then run CMake using:
```
cmake path/to/Xyce
```
Then, to build and install Xyce, run:
```
cmake --build . -j 2 -t install
```
The "-j 2" designates the number of processors being used for compiling Xyce.
Choose an appropriate number for your system.

#### Adding the Xyce/ADMS capability

To compile Verilog-A models using the "Xyce/ADMS" capability,
[ADMS](https://github.com/Qucs/ADMS) must be installed first. See the
[Xyce/ADMS Users Guide](https://xyce.sandia.gov/documentation/XyceADMSGuide.html)
for more information. To enable the capability in the CMake build, add the
following flag to the Xyce CMake invocation:
```
-D Xyce_PLUGIN_SUPPORT=ON
```

## Uninstalling Xyce

The Xyce CMake does not create uninstall script. However, on installation it
does produce an "install manifest" file, which lists the path to every
installed file. To remove a Xyce installation from a Unix-like system, simply
run:
```
xargs rm < install_manifest.txt
```
You may also want to uninstall Trilinos. As with Xyce, an
`install_manifest.txt` file is produced when Trilinos is installed. If you do
not want to keep the Trilinos build directory, simply copy that file to a safe
location (the file name can be changed). Then Trilinos can be uninstalled from
a unix-like system at any time using the `xarg` command, above.

## Running the Test Suite

After Xyce is installed, and you wish to test the installation, run the
[Xyce Regression Suite](https://github.com/Xyce/Xyce_Regression). See the
[Running the Xyce Regression Suite](https://xyce.sandia.gov/documentation/RunningTheTests.html)
documentation on the Xyce home page.

## Configuration Scripts

If many flags are being applied to CMake (for Trilinos or Xyce), you might want
to create a configuration script. This is essentially a shell script with the
CMake invocation. For example, when building Trilinos one might create a file
containing:
```sh
#!/bin/sh

cmake \
-C path/to/Xyce/cmake/trilinos/trilinos-config.cmake \
-D CMAKE_C_COMPILER=mpicc \
-D CMAKE_CXX_COMPILER=mpicxx \
-D CMAKE_Fortran_COMPILER=mpifort \
-D CMAKE_INSTALL_PREFIX=$HOME/install/Trilinos/mpi/opt \
<path/to/Trilinos>
```
More options can be added as needed. The file can be made executable by running
the following on the command line:
```sh
chmod a+x <filename>
```
Running the configuration script in the build directory is the came as running
the CMake configure line directly.

