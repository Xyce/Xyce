# Running the Expression Library Unit Tests

## TL;DR

To run the unit tests independently of a full Xyce build, do the following:
```bash
$ mkdir build
$ cd build
$ cmake \
  -D Trilinos_ROOT="/Path/to/Trilinos" \
  -D BUILD_ONLY_UNIT_TESTS=TRUE \
  -D CMAKE_BUILD_TYPE=Debug \
  ..
$ cmake --build . -j 2
$ ctest
```
The tests executables can be run individually as well.

## Explanation

The CMakeLists.txt file basically has two independent ways of building the
Expression Library test suite. The first method, defined in the 
`if(GTest_FOUND)` block, is enabled if 'BUILD_TESTING' is defined in the main
build and CMake is able to locate the GoogleTest framework. This method links
into the `XyceLib` library.

The AST tree uses templates, so it needs a type specified in order to use it;
but, because of the inflexibility of flex/bison, the type needs to be specified
at compile time rather than run time. Currently this is done via the
`USE_TYPE_DOUBLE` variable. It is possible that using flex/bison variants will
allow templates to be used in a more complete way and then the type can be
specified at runtime. At this time, Xyce _requires_ the expression library to
be in its "complex" form. The "complexParserUnitTest" program is written to
test the expression library using `std::complex<double>`. The "parserUnitTest"
program is written to test the expression library using doubles. Therefore, the
tests will not pass when linking to XyceLib, as with the other test programs.
To resolve this, the "parserUnitTest" setup forces a reconstruction/recompile
of the expression parser with `USE_TYPE_DOUBLE` defined, thus overriding the
parser in XyceLib.

The second building method is in the `if(BUILD_ONLY_UNIT_TESTS)` block. It is a
completely separate CMake project, and `BUILD_ONLY_UNIT_TESTS` should never be
enabled within a full Xyce build (the rest of the Xyce CMake knows nothing of
`BUILD_ONLY_UNIT_TESTS`). It compiles the minimum required files to be able to
run the unit test programs. To reduce compilation times, the non-test files are
compiled into a library to which the "complex"-based tests link. The
`parserUnitTest.C` file is compiled separately with the `-DUSE_TYPE_DOUBLE`
flag.

Since CMake assumes an out-of-source build, `cmake...` should be run from a
separate directory. The easiest way to do this is by creating a "build"
directory within `ExpressionTest`. Since the tests are dependent on Trilinos,
the Trilinos install directory must be specified with `Trilinos_ROOT`. Also,
`BUILD_ONLY_UNIT_TESTS=TRUE` must be specified to indicate the non-Xyce build.
Then build and test using the standard commands.

A handful of tests access data files. The files are currently in the following
directory structure:
```
ExpressionTest/SubDir1/test1.dat
ExpressionTest/Sub_Dir/1test_5.dat
```
The CMake script copies the files to the appropriate locations for the tests to
pass (`test1.dat` has to be moved to both the root build directory and in a
`SubDir1` subdirectory of the build directory). The reason for the two
sub-directories containing table file data is to ensure the parser can handle a
variety of pathnames and filenames.

## Miscellanea

The `*.gcov*` makefiles produce coverage information. The files lcovScript and
lcovScriptComplex walk through the entire coverage test process, including
generating html files and opening them in the browser. Nothing has been done to
reproduce those capabilities using CMake. It might make sense to update the
lcovScript and lcovScriptComplex script files so that they use the
CMake-generated buildfiles, which have been instrumented to include enable gcov
instrumentation. Also, see
[here](https://gitlab.kitware.com/cmake/cmake/-/issues/19942) for Kitware's
take. A Google search for "code coverage cmake" also brings up several
examples.

Once the CMake approach has been fully tested out, several files can
be removed:
- all the `Makefile*` files
- Test.C
- lcovScript
- lcovScriptComplex
- `parserCopyUnitTest.C` (Superseded by the other test programs. Can probably
  be removed at any time.)
