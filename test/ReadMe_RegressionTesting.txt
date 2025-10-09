
Xyce supports running Xyce_Regression tests through cmake's ctest. 

When you configure a Xyce build with cmake and with the optional setting
-DBUILD_TESTING=on, then, Xyce's top level CMakeLists.txt file will look
for a directory called Xyce_Regression at the same directory level as
Xyce.  Optionally, one can specify the location of Xyce_Regression with
the option -DXyce_REGRESSION_DIR set to the location of Xyce_Regression.
If that directory is found then cmake will attempt to add the regression
testing directory at this location in the BUILD directory.  All test
output will be in the BUILD directory.  Finally, if Google Test is
available and the variable, -DXyce_GTEST_UNIT_TESTS=ON is set, the Xyce
unit tests that use the GTest framework will also be included in the
build.

CMake variables that influence testing when cmake is used to configure a Xyce build:

- BUILD_TESTING : If set to "on" then cmake will try to build unit tests and look for 
  the Xyce_Regression directory for all of the regression tests.
  
- Xyce_REGRESSION_DIR : A location on the filesystem where Xyce_Regression can be found

- Xyce_GTEST_UNIT_TESTS : if set to "on" and if the Google Test library is found then
  include Xyce unit tests that utilize GTest in the build.  
  

All of the tags used to annotate Xyce_Regression's tests will be added as labels 
in the ctest system.  Any options set on the testing time limit will be set as
ctest TIMEOUT option for a given test.

Note, some Xyce_Regression tests have tags of the format "required:xxx"
where this test would not be run if the required component was not
available.  Usually this was some external program that was used as part
of the test.  During the CMake configure process these dependencies will
be probed and only tests that have all their required dependencies met
will be run.  For example, some tests require both Python and the Python
package SciPy. If those are not found then the test will be included in
testing run.

After configuring Xyce with cmake and building it with something like 
"cmake --build . --target all" and Xyce has successfully built, you can run tests with
ctest.  Here are some helpful guidlines:

ctest has the following useful command line options (there are many more options. 
see "ctest --help")

  -N,--show-only[=format]      = Disable actual execution of tests.  
                                 Dry run to show which tests will be run.

  -L <regex>, --label-regex <regex>
                               = Run tests with labels matching regular
                                 expression.  With multiple -L, run tests
                                 where each regular expression matches at
                                 least one label.
                                 
  -LE <regex>, --label-exclude <regex>
                               = Exclude tests with labels matching regular
                                 expression.  With multiple -LE, exclude
                                 tests where each regular expression matches
                                 at least one label.
                                 
  -R <regex>, --tests-regex <regex>
                               = Run tests matching regular expression.
                               
  -E <regex>, --exclude-regex <regex>
                               = Exclude tests matching regular expression.
                               
  -j [<level>], --parallel [<level>]
                               = Run tests in parallel, optionally limited to
                                 a given level of parallelism.
                                 
  --print-labels               = Print all available test labels.
  
  -V,--verbose                 = Enable verbose output from tests.
  
  --no-label-summary           = Don't print a summary of all the labels encountered after 
                                 a successful test run.
  
  
So, "ctest -L nightly -L resistor -LE required " would run "nightly" tests that have the
tag "resistor" and exclude tests that have a tag with "required" in the name 
(like required:vpi, required:xdm)

Running with the "-N" option lists the 300+ tests that would be run:

ctest -N -L nightly -L resistor 
Test project /Users/rlschie/src/XyceDevelopment/BUILD/Normal
  Test   #37: BSIMCMG/ac.cir.sh
  Test   #38: ABM_EXPLN/exp_const.cir
  Test   #40: ABM_EXPLN/exp_ln.cir
  Test   #43: ABM_FREQ/RC_simple_hertz.cir.sh
  Test   #44: ABM_FREQ/RC_simple.cir.sh
...
  Test #3570: BSIMCMG_110/ac.cir.sh
  Test #3593: ABM_SQRT/sqrt.cir.sh
Total Tests: 316

Note, the Test number as in #37 is the ctest index number for this test and is generated
at cmake configuration time.  It can be used to run a specific test with the "-I" option

-I [Start,End,Stride,test#,test#|Test file], --tests-information
                               = Run a specific number of tests by number.

as in:
ctest -I 37

The text in ctest's output after the test number is the test name.  One can use the test
name to run a test directly too with the "-R" option 

  -R <regex>, --tests-regex <regex>
                               = Run tests matching regular expression.
as in:
ctest -R ABM_EXPLN/exp_const.cir

or

ctest -R ABM_EXPLN

for all the tests in the ABM_EXPLN directory.

Another useful option is "--print-labels"  Running  "ctest --print-labels" will output
a list of all test labels in alphabetical order.

Note:  Xyce's unit tests were not originally part of the regression testing system 
that used tags which have been converted to labels.  All of Xyce's unittests now
have the label unittest.  So "ctest -N -L unittest " will show the list of available
unittests (over 1600 if Google Test is available at build time.)

When using the options like "-L" and "-LE" keep in mind that ctest applies a logical AND 
to the results.  So "-L resistor" "-L capacitor" will select tests that have BOTH the 
resistor AND the capacitor tag.  If you want to do an OR try -L "resistor|capacitor".  
Experimenting and using the "-N" option will help you select a reasonable number of tests 
to run.

Using "-j" to run tests in parallel is a great way to speed up tests.  But at this time
some tests may fail.  This can occur if two tests use some of the same underlying Xyce
input files or if a test requires the output of a prior test.  We will attempt to fix this
as the problems are identified.  For now try running "ctest --rerun-failed ". If the 
test then passes, it's likely that the test was conflicting with another test run in 
parallel.

