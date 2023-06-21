# this is meant to be invoked via jenkins and assumes that jenkins has
# cloned/updated the source

# hard-code most things for now, later will convert to arguments
# passed in via the ctest command line.

# screen output level
#   0 - almost no screen output (default)
#   5 - all screen output available
if(NOT VERBOSITY)
  set(VERBOSITY 0)
endif()

# WORKSPACE is an environment variable set by jenkins
set(CTEST_SOURCE_DIRECTORY "$ENV{WORKSPACE}/source/Xyce")

# the specified directory must exist or ctest with error out
set(CTEST_BINARY_DIRECTORY "$ENV{WORKSPACE}/build")

# use the "hostname" command for use as the CTEST_SITE variable, which
# is used in the "Site" column on the dashboard
find_program(HNAME NAMES hostname)
execute_process(COMMAND "${HNAME}"
  OUTPUT_VARIABLE CTEST_SITE
  OUTPUT_STRIP_TRAILING_WHITESPACE)

# make sure the custom xyce regression testing script is found
find_program(XYCE_REGR_SCRIPT NAMES run_xyce_regression
  PATHS ${WORKSPACE}/tests/Xyce_Regression/TestScripts
  REQUIRED)
  
# Release or Debug
set(CMAKE_BUILD_TYPE "Release")

# this is used as the "Build Name" column on the dashboard
set(CTEST_BUILD_NAME "glh_Intel64_RHEL7_Serial-cmake-cde-gcc-librarytests")

# used for invocation of parallel make
set(CTEST_BUILD_FLAGS "-j16")

set(CTEST_CMAKE_GENERATOR "Unix Makefiles")

set(MODEL "Experimental")

# the following are likely pretty invariant
set(CTEST_PROJECT_NAME "Xyce")
set(CTEST_DROP_METHOD "http")
set(CTEST_DROP_SITE "joseki-srn.sandia.gov/CDash")
set(CTEST_DROP_LOCATION "/submit.php?project=Xyce")

# begin ctest procedures. MODEL should be one of Nighlty, Weekly, or
# Experimental. this can use custom categories via the GROUP option to
# ctest_start() if desired
ctest_start(${MODEL})

# this runs cmake on xyce
ctest_configure()

# this runs make
ctest_build()

# run the custom xyce regression test script
execute_process(COMMAND ${XYCE_REGR_SCRIPT}
  --output=${WORKSPACE}/build/Xyce_Regression/
  --xyce-test=${WORKSPACE}/tests/Xyce_Regression/
  --xyce_verify=${WORKSPACE}/tests/Xyce_Regression/TestScripts/xyce_verify.pl
  --ignoreparsewarnings 
  --taglist="+serial?klu?weekly?nightly-verbose?noverbose?nonfree?rad?qaspr?athena?fft?stokhos?amesos2basker?amesos2klu2?xdm+library"
  --resultfile=${WORKSPACE}/build/regr_test_results_all
  ${WORKSPACE}/build/src/Xyce
  RESULT_VARAIBLE xyce_reg_result)

# submit results to the dashboard
ctest_submit(BUILD_ID cdash_bld_id
  RETRY_COUNT 10 RETRY_DELAY 30)

message("Build ID for cdash: ${cdash_bld_id}")
