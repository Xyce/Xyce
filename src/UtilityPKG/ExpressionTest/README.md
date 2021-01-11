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
Expression Library test suite. The first is in the `if(Xyce_TEST_SUITE)` block.
That links into the `XyceLib` library, and `Xyce_TEST_SUITE` is a known option
in the main build. In this mode, the tests in "parserUnitTest" are not run,
because it needs to be compiled with `-DUSE_TYPE_DOUBLE`. The `USE_TYPE_DOUBLE`
flag sets the `usedType` typedef to `double`, whereas the `usedType` default is
`std::complex<double>`.

The AST tree uses templates, so it needs a type specified in order to use it;
but, because of the inflexibility of flex/bison, the type needs to be specified
at compile time rather than run time. Currently this is done via the
`USE_TYPE_DOUBLE` variable. It is possible that using flex/bison variants will
allow templates to be used in a more complete way and then the type can be
specified at runtime. For the time being, though, Xyce requires the expression
library to be in its "complex" form, so "parserUnitTest" cannot currently be
run with the full Xyce build. The other parser unit test program,
"complexParserUnitTest" is written to test out the expression library with
`std::complex<double>`.

The second building method is in the `if(BUILD_ONLY_UNIT_TESTS)` block. It is a
completely separate project, and `BUILD_ONLY_UNIT_TESTS` should never be
enabled within a full Xyce build (the rest of the Xyce CMake knows nothing of
`BUILD_ONLY_UNIT_TESTS`). It compiles the minimum required files to be able to
run the unit test programs. To reduce compilation times, the non-test files are
compiled into a library to which the "complex"-based tests link. The
`parserUnitTest.C` file is compiled separately with the `-DUSE_TYPE_DOUBLE`
flag.

Since CMake assumes an out-of-source build, cmake should be run from a separate
directory. The easiest way to do this is by creating a "build" directory within
`ExpressionTest`. Since the tests are dependent on Trilinos, the Trilinos
install directory must be specified with `Trilinos_ROOT`. Also,
`BUILD_ONLY_UNIT_TESTS=TRUE` must be specified to indicate the non-Xyce build.
Then build and test using the standard commands.

## Miscellanea

A handful of tests access data files. The files are currently in the following
directory structure:
```
ExpressionTest/SubDir1/test1.dat
ExpressionTest/Sub_Dir/1test_5.dat
```
The CMake script copies the files to the appropriate locations for the tests to
pass. However, the reason for the current structure is not clear. Also,
`test1.dat` has to be moved to both the root build directory and in a `SubDir1`
subdirectory of the build directory.

The `parserCopyUnitTest.C` has been superseded by the other test programs. It
is not built by CMake. Once it has been deemed unnecessary, it should be
removed.

Finally, once the CMake approach has been fully tested out, several files can
be removed:
```
Makefile
Makefile.bsd
Makefile.cmplx
Makefile.cmplx.gcov_gcc8
Makefile.cmplx.valinor
Makefile.gcc8
Makefile.gcc9
Makefile.gcov_gcc8
Makefile.gcov_gcc9
Makefile.reflex
Makefile.valinor
Test.C
lcovScript
lcovScriptComplex
```
