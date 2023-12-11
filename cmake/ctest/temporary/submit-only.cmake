# this is meant to be invoked via jenkins and assumes that jenkins has
# cloned/updated the source

# arguments, specified via "-D"
#   -DVERBOSITY=<0-5>
#   -DDASHSUBMIT=<TRUE|FALSE>    # mostly for debugging to avoid cdash submission
#   -DCMAKE_ARGS_LIST="-DVAR1=VAL1;-DVAR2=VAL2;..."
#   -DCDASHVER=<version of cdash>  # should be either 3.1 or not set

cmake_minimum_required(VERSION 3.23)

# verbosity level
#   0 - no specific screen output (default)
#   5 - all screen output available
if(NOT VERBOSITY)
  set(VERBOSITY 0)
endif()

# default TRUE
if(NOT DEFINED DASHSUBMIT)
  set(DASHSUBMIT TRUE)
endif()

# the version of cdash matters for the custom Test.xml file that is
# generated
if(NOT DEFINED CDASHVER)
  set(CDASHVER 0.0)
endif()

set(CTEST_PROJECT_NAME "Xyce")

set(CTEST_DROP_METHOD "https")
set(CTEST_DROP_SITE "xyce-cdash.sandia.gov")
set(CTEST_DROP_LOCATION "/submit.php?project=Xyce")

set(MODEL "Nightly")

set(CTEST_BINARY_DIRECTORY "/gpfs/glhenni/xyce/pbpr-pbsr/workspace/build")
set(CTEST_SOURCE_DIRECTORY "/gpfs/glhenni/xyce/pbpr-pbsr/workspace/source/Xyce")

ctest_start(${MODEL})

#ctest_submit(FILES "/gpfs/glhenni/xyce/pbpr-pbsr/workspace/build/PbPr_regr_test_results_all.test.xml"
#  "/gpfs/glhenni/xyce/pbpr-pbsr/workspace/build/regr_test_results_all.test.xml")

ctest_submit(FILES "/gpfs/glhenni/xyce/pbpr-pbsr/workspace/build/regr_test_results_all.xml")
