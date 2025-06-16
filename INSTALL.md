# Xyce Configure, Build and Installation Guide Using CMake

> **Note:** _The Autotools configuration and build system was DEPRECATED in Xyce 7.10_

------------------------------------------------------------------------
------------------------------------------------------------------------

This guide describes the basic process for compiling and installing a Xyce
binary using the CMake build system. It is easiest to view these instructions
with full formatting on the [Xyce GitHub website](https://github.com/Xyce/Xyce/blob/master/INSTALL.md).
For instructions on building Xyce with the deprecated Autotools system, see the
[Xyce Building Guide - DEPRECATED](https://xyce.sandia.gov/documentation-tutorials/building-guide/).

Xyce can be built for serial execution or with distributed-memory (MPI)
parallelism. This is determined by the parallelism enabled in the build of
Trilinos used when configuring Xyce. Different variants of Xyce and Tilinos can
co-exist on a system, but the variants must be installed in different
directories.

__Binary installers for Windows, macOS and Red Hat Linux are made available for
the latest release of Xyce on the [Xyce website](https://xyce.sandia.gov)__.
Installers for previous versions of Xyce may be made available upon request at
the Xyce [contact-us](https://xyce.sandia.gov/contact-us/) page.

## Overview

Building and installing Xyce from source involves the following steps:
1. Install the [prerequisite](#prerequisites) libraries and tools
2. [Build and install](#building-trilinos) Trilinos from the
   [source code](https://github.com/trilinos/Trilinos)
3. [Build and install](#building-xyce) Xyce from the Xyce
   [website](https://xyce.sandia.gov/downloads/source-code) or
   [GitHub](https://github.com/Xyce/Xyce)

> **Note:** Some install/uninstall commands may require administrative
> privileges on Unix-like systems.

## Prerequisites

Xyce depends on external libraries to supply enhanced parsing, mathematical
algorithms and parallel communication. Most Linux distributions will have the
required components already installed, or they can be added via the package
manager. For Mac and Windows, refer to the
[System-Specific Modifications](#system-specific-modifications) section for
guidance. If needed, the prerequisites can be built from source.

Be sure all the prerequisites are installed prior to [building Trilinos](#building-trilinos).

| Dependency         | Version               | Required | Notes |
| :----------------  | :-------------------- | :------: | :---  |
| CMake              | 3.22 or later         |   \*No   | CMake automates the Xyce build process. (\*[CMake](https://cmake.org) will be required in release 7.11.) |
| Build system       | various               |   Yes    | CMake requires a build system, such as [Make](https://www.gnu.org/software/make/), [Ninja](https://ninja-build.org/), or [Jom](https://wiki.qt.io/Jom). |
| C/C++ compiler     | various               |   Yes    | Any modern C++17 compliant compiler—such as gcc, Clang, or Intel—should work. |
| Fortran compiler   | various               |   No     | Trilinos has Fortran code, but Xyce does not. It is technically [optional](https://docs.trilinos.org/files/TrilinosBuildReference.html#disabling-the-fortran-compiler-and-all-fortran-code), but not using the Fortran code could result in slower performance. Both the gcc and Intel compiler suites include a Fortran compiler. |
| MPI compiler       | various               |   No     | This only needed for a distributed-memory parallel build of Xyce. [Open MPI](http://www.open-mpi.org) is recommended. ([MPICH](https://www.mpich.org/) should also work, but the Xyce project does not regularly build with that library.) |
| bison              | 3.3 or later          |   Yes    | [Bison](http://www.gnu.org/software/bison/) is used for various parsing/lexing in Xyce. |
| flex               | 2.6 or later          |   Yes    | [flex](https://github.com/westes/flex) is used for various parsing/lexing in Xyce. |
| BLAS               | any                   |   Yes    | The [BLAS](http://www.netlib.org/blas) numerical package is used via Trilinos. |
| LAPACK             | any                   |   Yes    | The [LAPACK](http://www.netlib.org/lapack) numerical package is used via Trilinos. |
| FFT                | FFTW 3.x or Intel MKL |   No     | This is required for Harmonic Balance analysis, which will be disabled without an FFT library. Xyce can use [FFTW](http://www.fftw.org) or the [Intel MKL](http://software.intel.com). |
| SuiteSparse        | 7.8.3 or later        |   Yes    | Xyce uses AMD, which is part of the SuiteSparse library. It enabled via Trilinos. If needed, [download](https://github.com/DrTimothyAldenDavis/SuiteSparse), compile, and install the library according to the [Building SuiteSparse](#building-suitesparse) section, below. |

### Building SuiteSparse

<details>
<summary>Click here for build instructions</summary>

If SuiteSparse is not provided by a package manager, the more recent versions
are easy to build with CMake. SuiteSparse version 7.8.3 or later is
recommended. Xyce depends only on the AMD and SuiteSparse_config packages in
SuiteSparse.

To build and install AMD, first, download SuiteSparse from
[GitHub](https://github.com/DrTimothyAldenDavis/SuiteSparse).
Next, create a build directory and go into it. Finally enter the following,
at the command line, replacing the items in angle brackets with the appropriate
values. We recommend installing AMD to the same location as Trilinos. 
```sh
cmake \
-D CMAKE_INSTALL_PREFIX=<path/to/where-you-will-install-Trilinos> \
-D SUITESPARSE_ENABLE_PROJECTS="suitesparse_config;amd" \
<path/to/SuiteSparse> 

cmake --build -t install
```
</details>

### Broken Open MPI on Older Linux Distros

<details>
<summary>Click here for more information</summary>

Some Linux distributions from between 2017 and 2020 have broken versions of
Open MPI in their package repositories. We are not aware of continuing
problems in newer releases of these Linux systems. In the problem releases,
the version of Open MPI in the repositories is compiled with the
`--enable-heterogeneous` option, which breaks MPI's standard compliance and
causes Xyce to fail in some situations.

We have encountered with this issue with Ubuntu, beginning with release 17.10
and continuing up to, but not including, 20.04 LTS. See comment 11 of the
[Launchpad bug report](https://bugs.launchpad.net/ubuntu/+source/openmpi/+bug/1731938/comments/11).

To check whether your system's install of Open MPI has this issue, you can run
a small test program from the Xyce source tree. To test your Open MPI package,
copy the file, [Xyce/src/test/MPITest/testBUG967.c](src/test/MPITest/testBUG967.c),
into temporary directory. Then compile and run it using the following commands:
```sh
mpicc -o testBUG967 testBUG967.c
mpirun -np 2 ./testBUG967
```
If the run produces any output with the word __BAD__, your Open MPI install is
broken and cannot be used. 

</details>

## Building Trilinos

> **Note:** [Trilinos](https://trilinos.github.io/) is available in some
> package managers, but the particular installation may not have all the
> features required by Xyce. If Trilinos is available, you can try to use it to
> configure and build Xyce. If any issue is encountered, build Trilinos using
> these instructions.

> **Note:** To build Trilinos on Windows, see the Windows section under
> [System-Specific Modifications](#system-specific-modifications).

Trilinos is an extensive set of numerical packages for a wide range of
computational problems. This section will provide a general overview of
building Trilinos for Xyce. For detailed questions on Trilinos or its build
system see the Trilinos
[getting started page](https://trilinos.github.io/getting_started.html).

The current minimum version of [Trilinos](https://github.com/trilinos/Trilinos)
usable in a CMake-configured Xyce is 14.4. __Versions after 16.1 have not been
rigorously tested with Xyce and may not work properly.__ A zip file of Trilinos
version 14.4 can be downloaded
[here](https://github.com/trilinos/Trilinos/archive/refs/heads/trilinos-release-14-4-branch.zip).

The following process will produce a serial Trilinos installation that will
contain only the libraries needed by Xyce. For a distributed-memory parallel
build, be sure to understand the serial build process prior to reading the
[Building Trilinos with MPI Parallelism](#building-trilinos-with-mpi-parallelism) section.

It is recommended to build the necessary packages outside the Trilinos source
code directories. 
1. Create a directory in a convenient location where you will build a serial
   version of Trilinos. 
2. Identify a convenient location for the Trilinos installation. On Unix-like
   systems, the default Trilinos installation location is `/usr/local`.
   Multiple installations of Trilinos can exist on the same system, but they
   must be in different directories. We recommend specifying unique
   sub-directories in `/usr/local`, such as `/usr/local/trilinos_serial`. The
   installation directory can be changed by adding the following flag to the
   CMake invocation.
   ```sh
   -D CMAKE_INSTALL_PREFIX=<path/to/where-you-will-install-Trilinos> \
   ```
3. Verify the location of your compilers. If you have compilers or libraries in
   non-standard locations, see the [Other Trilinos Options](#other-trilinos-options)
   section, below.

A CMake "initial cache" file, called `trilinos-base.cmake`, is included in the
Xyce repository in the `cmake/trilinos` directory. The file contains a typical
set of options for a Xyce-oriented serial Trilinos build.

To configure Trilinos (using the default `/usr/local` install location), enter
the build directory and run:
```sh
cmake -C <path/to/Xyce>/cmake/trilinos/trilinos-base.cmake <path/to/Trilinos>
```
Once the configuration step has completed, run the following in the build
directory to build and install Trilinos:
```sh
cmake --build . -j 2 -t install
```
The "-j 2" designates the number of processors to be used for compiling
Trilinos. Choose an appropriate number for your system.

#### Other Trilinos Options

<details>
<summary>Click here for more details about building Trilinos</summary>

The following are various flags that might be needed when building Trilinos.
These can be added to the command line, or put in a
[configuration script](#configuration-scripts).

CMake will use the first compiler set it finds on your system. You can specify
the compilers by adding the following flags to the CMake invocation:
```sh
-D CMAKE_C_COMPILER=<C-compiler> \
-D CMAKE_CXX_COMPILER=<C++-compiler> \
-D CMAKE_Fortran_COMPILER=<Fortran-compiler> \
```
You may need to use a full path if they cannot be located in your default paths.

If no Fortran compiler is available you can use:
```sh
-DTrilinos_ENABLE_Fortran=OFF \
```
to disable any Fortran-dependent code in Trilinos.

On some systems it may be necessary to add the following to your CMake
invocation if the installed C compiler treats implicitly defined functions as
errors, such as Apple's Clang.
```sh
-DCMAKE_C_FLAGS="-Wno-error=implicit-function-declaration" \
```

Finally, if the third-party libraries (AMD, BLAS and LAPACK) cannot be located
in your default paths, use the following flags to help CMake find the
libraries:
```sh
-D AMD_LIBRARY_DIRS=/path/to/AMD/lib \
-D AMD_INCLUDE_DIRS=/path/to/AMD/include \
-D BLAS_LIBRARY_DIRS=/path/to/BLAS/lib \
-D LAPACK_LIBRARY_DIRS=/path/to/LAPACK/lib \
```
See the [Trilinos Build Reference](https://docs.trilinos.org/files/TrilinosBuildReference.html)
for more details on the Trilinos build options.
</details>

#### Building Trilinos with MPI Parallelism

<details>
<summary>Click here for the MPI parallel Trilinos building instructions</summary>

To enable MPI parallelism in Xyce, Trilinos must be built with MPI enabled. An
"initial cache" file, called `trilinos-MPI-base.cmake`, is provided for MPI
builds. In addition to using the MPI-oriented cache file, the MPI compilers
must be explicitly specified to CMake. The following CMake invocation should
work on most systems:
```sh
cmake \
-C <path/to/Xyce>/cmake/trilinos/trilinos-MPI-base.cmake \
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
</details>

## Building Xyce

> **Note:** If you plan to have multiple builds of Xyce on your system, they
> must be in different directories. We recommend specifying unique
> sub-directories in `/usr/local`, such as `/usr/local/xyce_serial` and
> `/usr/local/xyce_mpi`.

The generalized process for a serial configuration of Xyce, assuming Trilinos
is already [installed](#building-trilinos), can be summarized as:

```sh
cd <your-build-directory>
mkdir xyce-build
cd xyce-build

cmake \
-D CMAKE_INSTALL_PREFIX=<path/to/where-you-will-install-Xyce> \
-D Trilinos_ROOT=</path/to/Trilinos-install-location> \
<path/to/Xyce>
```

The simple CMake invocation above assumes that:
- All the required [prerequisites](#prerequisites) are found in your default paths  
- Trilinos is installed in the `Trilinos_ROOT` directory
- Xyce will be installed in the `CMAKE_INSTALL_PREFIX` directory
- The first compiler found in your default paths is the __same__ one used to
  build Trilinos.

If a different C/C++ compiler was used to build Trilinos, specify that compiler
using these commands:
```sh
-D CMAKE_C_COMPILER=<C-compiler> \
-D CMAKE_CXX_COMPILER=<C++-compiler> \
```
You may need to use a full path if the compilers cannot be located in your
default paths. If additional compiler flags need to be passed to the C/C++
compiler, use these commands:
```sh
-D CMAKE_C_FLAGS=<flags> \
-D CMAKE_CXX_FLAGS=<flags> \
```

If a required prerequiste is not found in your default path, then you must
specify the location of those headers, libraries, or executables. Xyce provides
special CMake options for some of these prerequisites, like flex and bison:
```sh
-D FLEX_EXECUTABLE=<path/to/flex-install-location>/bin/flex \
-D FLEX_INCLUDE_DIR=<path/to/flex-install-location>/include \
-D BISON_EXECUTABLE=<path/to/bison-install-location>/bin/bison \
```

Once the configuration step is done, run the following in the `xyce-build`
directory to compile and install Xyce in the `CMAKE_INSTALL_PREFIX` directory:
```sh
cmake --build . -j 2 -t install
```
The `-j 2` indicates that two processors should be used for compiling Xyce.
Choose an appropriate number for your system.

### Building Xyce with MPI Parallelism

The generalized process for a parallel configuration of Xyce, given that
MPI-enabled Trilinos libraries are already [installed](#building-trilinos), can
be summarized as:
```sh
cd <your-build-directory>
mkdir xyce-mpi-build
cd xyce-mpi-build

cmake \
-D CMAKE_INSTALL_PREFIX=<path/to/where-you-will-install-mpi-Xyce> \
-D Trilinos_ROOT=</path/to/Trilinos-mpi-install-location> \
-D CMAKE_C_COMPILER=<mpi-C-compiler> \
-D CMAKE_CXX_COMPILER=<mpi-C++-compiler> \
<path/to/Xyce>
```

The simple CMake invocation above assumes that:
- All the required [prerequisites](#prerequisites) are found in your default paths  
- Trilinos is installed in the `Trilinos_ROOT` directory
- Xyce will be installed in the `CMAKE_INSTALL_PREFIX` directory
- The MPI compilers provided to `CMAKE_C_COMPILER` and `CMAKE_CXX_COMPILER`
  are the __same__ ones used to build Trilinos.

If a required prerequisite is not found in your default path, like flex or
bison, then you must specify the location of those executables (see the serial
[Building Xyce](#building-xyce) section, above).

As with the serial build, once the configuration step is done, run the
following in the `xyce-build` directory to compile and install Xyce in the
`CMAKE_INSTALL_PREFIX` directory:
```sh
cmake --build . -j 2 -t install
```  

### Adding the Xyce/ADMS Verilog-A Model Compiler

Xyce has a Verilog-A model compiler capability, which uses the "Xyce/ADMS"
compiler tool. See the
[Xyce/ADMS Users Guide](https://xyce.sandia.gov/documentation-tutorials/xyce-adms-users-guide/)
for more information on the capability and for instructions on using Xyce/ADMS.

To enable the feature with CMake, install [ADMS](https://github.com/Qucs/ADMS)
prior to building Xyce. Then, to enable the capability in the Xyce build, add
the following flag to the Xyce CMake invocation:
```sh
-D Xyce_PLUGIN_SUPPORT=ON \
```
The CMake support for the plugin capability is still being developed. As such,
there are some differences from the website:
- The "toys" example is installed in `/path/to/install/share/examples/toys`.
- The `buildxyceplugin.sh` script requires the absolute path to the `.va` files.
- The name of the plugin file can vary by system (e.g., on a Mac, the "toys"
  library will be called, "libtoys.dylib", not "toys.so").

## Uninstalling Xyce

The Xyce CMake does not create an uninstall script. However, on installation it
does produce an "install manifest" file, which lists every installed file with
their full path. To remove a Xyce installation from a Unix-like system, simply
run:
```sh
xargs rm < install_manifest.txt
```
Note that the above method is immediate and permanent. You may also want to
uninstall Trilinos. As with Xyce, an `install_manifest.txt` file is produced
when Trilinos is installed.

If you do not want to keep the build directories, simply copy the
`install_manifest.txt` file(s) to a safe location (the file name can be
changed). Then Xyce and/or Trilinos can be uninstalled at any time using the
`xarg` command, above (be sure to specify the correct filename).

## Running the Test Suite

If you wish to test the Xyce installation, run the
[Xyce Regression Suite](https://github.com/Xyce/Xyce_Regression). See the
[Running the Xyce Regression Suite](https://xyce.sandia.gov/documentation-tutorials/running-the-xyce-regression-suite/)
documentation on the Xyce home page. Note that the test suite is controlled
with Perl and Bash scripts, and some tests require Python with Scipy and Numpy.

## Configuration Scripts

If many flags are being applied as part of a CMake invocation, you might want
to create a configuration script. This is essentially a shell script with the
CMake invocation. The following commands are for Unix-like systems. Analogous
scripts can be created for Windows systems.

As an example, to configure an MPI-enabled Trilinos build, one might create a
file containing:
```sh
#!/bin/sh

cmake \
-C path/to/Xyce/cmake/trilinos/trilinos-MPI-base.cmake \
-D CMAKE_C_COMPILER=mpicc \
-D CMAKE_CXX_COMPILER=mpicxx \
-D CMAKE_Fortran_COMPILER=mpifort \
-D CMAKE_INSTALL_PREFIX=$HOME/install/Trilinos/mpi/ \
path/to/Trilinos
```
More options can be added as needed. The file can be made executable by running
the following on the command line:
```sh
chmod a+x <filename>
```
Running the configuration script in the build directory is the same as running
the CMake configure line directly.

## System-Specific Modifications

The instructions above are applicable to Linux and similar systems. However,
additional instructions are provided here for Windows and macOS. Click on the
following for more information.

<details>
<summary>Windows</summary>

*Build tools:*\
Compiling Xyce on Windows is not a small task at the moment, primarily because
Windows does not have the equivalent of a package manager. The Xyce team uses
the Intel compiler suite with the Intel
[Math Kernel Library](https://software.intel.com/content/www/us/en/develop/tools/oneapi/components/onemkl.html)
(MKL). The MKL provides the BLAS, LAPACK and FFT capabilities (removing the need
for FFTW). At the beginning of 2021, Intel rebranded their tool chains as the
[oneAPI](https://software.intel.com/content/www/us/en/develop/tools/oneapi.html)
Toolkits, and makes them available for free. The
[oneAPI Base Toolkit](https://software.intel.com/content/www/us/en/develop/tools/oneapi/base-toolkit.html)
is sufficient for building Xyce and Trilinos. Running the setvars.bat
(Component Directory Layout) or oneapi-vars.bat (Unified Directory Layout)
script that installs with oneAPI can help ensure that the necessary environment
variables and paths are set for building.

Note that the oneAPI toolkit requires a compatible version of
[Microsoft C++ BuildTools](https://visualstudio.microsoft.com/visual-cpp-build-tools/)
to be installed first to provide required system libraries. Check the 
[Intel compiler compatibility site](https://www.intel.com/content/www/us/en/developer/articles/reference-implementation/intel-compilers-compatibility-with-microsoft-visual-studio-and-xcode.html).
To find a suitable version of BuildTools, you may need to refer to one of the 
Long Term Servicing Channels for Visual Studio
[2022](https://learn.microsoft.com/en-us/visualstudio/releases/2022/release-history#release-dates-and-build-numbers)
or
[2019](https://learn.microsoft.com/en-us/visualstudio/releases/2019/history).
The Microsoft Build Tools will also supply CMake, as well as several supported 
build tool options (msbuild, nmake and ninja).

Trilinos does not test their code on Windows (versions 14.0 and earlier may be
unstable on Windows platforms). The initial configuration file for Trilinos in
`path/to/Xyce/cmake/trilinos/trilinos-base.cmake` may require a few extra
options:
```sh
cmake ^
    -C path\to\Xyce\cmake\trilinos\trilinos-base.cmake ^
    -D HAVE_TEUCHOS_LAPACKLARND=OFF ^
    -D Trilinos_ENABLE_Stokhos=OFF ^
    -D Trilinos_ENABLE_Sacado=OFF ^
    -D Trilinos_ENABLE_Amesos2=OFF ^
    ...
    path\to\Trilinos
```

On Windows, Bison and flex are available via the
[WinFlexBison](https://github.com/lexxmark/winflexbison) package. For
SuiteSparse, see the [Building SuiteSparse](#building-suitesparse) section
above, under [Obtaining the TPLs](#obtaining-the-tpls). Once all the TPLs are
in place, continue with the standard Trilinos and Xyce build processes.

To compile Xyce on Windows, there are additional CMake options that must be 
added in generating the build configuration.
```sh
-D Xyce_USE_FFT=TRUE
-D Xyce_USE_INTEL_FFT=TRUE
```

The Xyce CMake build system also supports a package target. On Windows, this 
target will generate a standalone installer for Xyce. The package target uses
[NSIS](https://nsis.sourceforge.io/Main_Page) (version 3.0 or newer) to generate 
the installer. This is an additional third-party dependency that needs to be 
installed prior to building the package target.

```sh
cmake --build <path to build folder> -t package
```
-- or --
```sh
cpack <path to build folder>
```

*Test environment:*  
The Xyce regression suite is intended to run in a Unix-like environment. The 
engine that runs the regression suite consists primarily of bash, perl, and (to 
a lesser degree) python scripts.

While there is no native support for running the regression suite on Windows, 
there are various Unix compatibility tools that enable running the test suite in
a Windows environment. The Xyce team uses Cygwin to run Xyce regression testing 
on Windows. Other configurations that leverage utilites such as WSL, Minigw, 
or MSYS2 may also be possible, but are untested.

Cygwin installs a minimal set of packages by default, so all of the dependencies 
will have to be explicitly added. To run the full regression suite, the 
recommended set of packages include:
- liblapack-devel
- libopenblas
- gcc-fortran
- gcc-g++
- perl
- perl-libwww-perl
- bash
- python39
- python39-numpy
- python39-devel
- python39-pip
- diff
- cmake

There are few tests in the regression suite that use the python package scipy. 
This package is not included in any of the available cygwin packages. After 
cygwin is in place, you can use pip to install scipy.
</details>

<details>
<summary>macOS</summary>

macOS provides a code development environment through
[Xcode](https://developer.apple.com/xcode/), which includes the Clang C/C++
compiler and BLAS/LAPACK libraries. It does not include a Fortran compiler. For
that, install gcc from one of the package managers (below) or use the
```sh
-D Trilinos_ENABLE_Fortran=OFF \
```
flag when doing the Trilinos CMake configuration.

macOS does not have a native package manager, but there are several third-party
package managers available. The Xyce team has had success with
[MacPorts](https://www.macports.org/) and [Homebrew](https://brew.sh/).

__Xyce/ADMS on macOS__\
The `buildxyceplugin.sh` script relies on the behavior of GNU "readlink", which
is very different than BSD readlink (macOS is based on BSD Unix). You can
obtain GNU readlink by installing the "coreutils" package from MacPorts or
Homebrew.
</details>
