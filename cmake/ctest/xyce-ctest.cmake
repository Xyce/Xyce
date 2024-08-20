# this is meant to be invoked via jenkins and assumes that jenkins has
# cloned/updated the source

# arguments, specified via "-D"
#   -DVERBOSITY=<0-5>
#   -DDASHSUBMIT=<TRUE|FALSE>    # mostly for debugging to avoid cdash submission
#   -DCMAKE_ARGS_LIST="-DVAR1=VAL1;-DVAR2=VAL2;..."  # arguments to pass directly to cmake
#   -DCDASHVER=<version of cdash>  # should be either 3.1 or not set
#   -DRXR_APPEND_TAGS="tags as used by run_xyce_regression script to add"
#   -DMPI_TESTING=<TRUE|FALSE>
#   -DTDEV_BUILD=<TRUE|FALSE>

cmake_minimum_required(VERSION 3.23)

# verbosity level
#   0 - no specific screen output (default)
#   5 - all screen output available
if(NOT VERBOSITY)
  set(VERBOSITY 0)
endif()

# default off
if(NOT ENABLE_TESTING)
  set(ENABLE_TESTING OFF)
endif()

# default TRUE
if(NOT DEFINED DASHSUBMIT)
  set(DASHSUBMIT TRUE)
endif()

# default is to NOT do testing against a Trilinos/develop build
if(NOT DEFINED TDEV_BUILD)
  set(TDEV_BUILD FALSE)
endif()

# default FALSE
if(NOT DEFINED MPI_TESTING)
  set(MPI_TESTING FALSE)
endif()

# the version of cdash matters for the custom Test.xml file that is
# generated
if(NOT DEFINED CDASHVER)
  set(CDASHVER 0.0)
endif()

if(TDEV_BUILD)
  set(xyceSitePostFix "Xyce TrilinosDev")
else()
  # NOT building against Trilinos/develop
  if($ENV{XYCE_MPI})
    if(MPI_TESTING)
      set(xyceSitePostFix "PbPr")
    else()
      set(xyceSitePostFix "PbSr")
    endif()
  else()
    set(xyceSitePostFix "")
  endif()
endif()

# error check
if(NOT DEFINED ENV{MYBUILDNAME})
  message(FATAL_ERROR "ERROR: Required environment varialble \"MYBUILDNAME\" not set")
endif()
if(NOT DEFINED ENV{branch})
  message(FATAL_ERROR "ERROR: Required environment varialble \"branch\" not set")
endif()
if(NOT DEFINED ENV{TESTSET})
  message(FATAL_ERROR "ERROR: Required environment variable \"TESTSET\" not set")
endif()

# WORKSPACE is an environment variable set by jenkins
set(CTEST_SOURCE_DIRECTORY "$ENV{WORKSPACE}/Xyce")

# the specified directory must exist or ctest will error out
set(CTEST_BINARY_DIRECTORY "$ENV{WORKSPACE}/build")

# this should probably be a variable in the environment or passed in,
# but for now it's hard-coded to 1am
set(CTEST_NIGHTLY_START_TIME "01:00:00 MDT")

# use the "hostname" command as the CTEST_SITE variable, which is used
# in the "Site" column on the dashboard
find_program(HNAME NAMES hostname)
execute_process(COMMAND "${HNAME}"
  OUTPUT_VARIABLE HOST_NAME
  OUTPUT_STRIP_TRAILING_WHITESPACE)

if(${HOST_NAME} MATCHES "^ascic[0-9]*")
  set(CTEST_SITE "ascic")
elseif(${HOST_NAME} MATCHES "cee-build[0-9]*")
  set(CTEST_SITE "cee-build")
else()
  set(CTEST_SITE ${HOST_NAME})
endif()

if(VERBOSITY GREATER 4)
  message("[VERB5]: ENV{SNLSYSTEM} = $ENV{SNLSYSTEM}")
endif()
if($ENV{SNLSYSTEM} MATCHES "^cts*")
  set(CTEST_SITE "cts")
endif()

# add any postfix to the site
set(CTEST_SITE "${CTEST_SITE} ${xyceSitePostFix}")

# this is used as the "Build Name" column on the dashboard
set(CTEST_BUILD_NAME "$ENV{MYBUILDNAME}")

# used for invocation of parallel make
if(DEFINED ENV{NUM_JOBS})
  set(CTEST_BUILD_FLAGS "-j$ENV{NUM_JOBS}")
else()
  set(CTEST_BUILD_FLAGS "-j8")
endif()

# note that "Weekly" is just a Nightly category with a different group
# name
if(NOT DEFINED ENV{TESTSET})
  message(FATAL_ERROR "ERROR: You must set the environment variable TESTSET to one of Nighlty, Weekly or Experimental")
endif()

if($ENV{TESTSET} STREQUAL "Nightly")
  set(MODEL "Nightly")
  set(TESTGROUP "Nightly")
elseif($ENV{TESTSET} STREQUAL "Weekly")
  set(MODEL "Nightly")
  set(TESTGROUP "Weekly")
else()
  set(MODEL "Experimental")
  set(TESTGROUP "Experimental")
endif()

set(CTEST_PROJECT_NAME "Xyce")

set(CTEST_DROP_METHOD "https")
set(CTEST_DROP_SITE "xyce-cdash.sandia.gov")
set(CTEST_DROP_LOCATION "/submit.php?project=Xyce")

if(VERBOSITY GREATER 4)
  message("[VERB5]: CMAKE_ARGS_LIST = ${CMAKE_ARGS_LIST}")
endif()
set(XYCE_CMAKE_CONF_ARG "")
foreach(cmakeopt IN LISTS CMAKE_ARGS_LIST)
  set(XYCE_CMAKE_CONF_ARG "${XYCE_CMAKE_CONF_ARG} ${cmakeopt}")
endforeach()
if(VERBOSITY GREATER 4)
  message("[VERB5]: XYCE_CMAKE_CONF_ARG = ${XYCE_CMAKE_CONF_ARG}")
endif()

# generate the cmake configuration command
find_program(CTEST_CMAKE_COMMAND cmake REQUIRED)
set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
set(CTEST_CONFIGURE_COMMAND "${CTEST_CMAKE_COMMAND} --log-level=Debug \"-GUnix Makefiles\"")
set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} ${XYCE_CMAKE_CONF_ARG}")
set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} ${XYCE_CMAKE_CONF_ARG}")
set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} -DENABLE_TESTING=${ENABLE_TESTING}")
set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} ${CTEST_SOURCE_DIRECTORY}")

if(VERBOSITY GREATER 1)
  message("[VERB1]: CTEST_CONFIGURE_COMMAND = ${CTEST_CONFIGURE_COMMAND}")
endif()

# begin ctest procedures. MODEL should be one of Nighlty, Weekly,
# Continuous or Experimental. this can use custom categories via the
# GROUP option to ctest_start() if desired
ctest_start(${MODEL} GROUP ${TESTGROUP})

# this runs cmake on xyce
ctest_configure()

# this runs make
ctest_build(RETURN_VALUE buildReturnVal)

# if the build succeeds, as indicated by a zero return value, proceed,
# otherwise skip to submission
if(buildReturnVal EQUAL 0)
  ctest_test(RETURN_VALUE testReturnVal
    PARALLEL_LEVEL $ENV{NUM_JOBS}
    INCLUDE_LABEL "nightly;serial"
    EXCLUDE_LABEL "^required:.*")
  if(VERBOSITY GREATER 1)
    message("[VERB1]: ctest_test() exited with return value: ${testReturnVal}")
  endif()
endif()

# submit results to the dashboard
if(DASHSUBMIT)
  ctest_submit(RETRY_COUNT 10 RETRY_DELAY 30)
endif()
