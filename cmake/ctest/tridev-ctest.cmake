# this is used to build Trilinos/develop which will later be utilized
# by the xyce build

# arguments:
#   -DVERBOSITY=<0-5>
#   -DTRACK=<Nightly|OnDemand>
#   -DBUILD_WITH_MPI=<ON|OFF>
#   -DBUILD_DIR=<subdirectory in which to actually build trilinos>
#   -DBUILD_NAME=<name of build to be displayed on dashboard>
#   -DCMAKE_ARGS_LIST="-DVAR1=VAL1;-DVAR2=VAL2;..."  # arguments to pass directly to cmake
#   -DSITE_APPEND=<string to append to the site variable, "TrilinosDev" for example>
#   -DOSX_TOOLKIT_VERSION=<toolkit version (see Xyce root CMakeLists.txt for value)>

if(NOT VERBOSITY)
  set(VERBOSIT 0)
endif()

cmake_minimum_required(VERSION 3.23)

set(CTEST_PROJECT_NAME "Xyce")

if(NOT DEFINED TRACK)
  set(TRACK "Nightly")
endif()

# you need this to match the one you'll be using when building xyce so
# if pass it in if it's not the default as set below
if(NOT DEFINED OSX_TOOLKIT_VERSION)
  set(OSX_TOOLKIT_VERSION "12.5")
endif()

set(CTEST_DROP_METHOD "https")
set(CTEST_DROP_SITE "xyce-cdash.sandia.gov")
set(CTEST_DROP_LOCATION "/submit.php?project=Xyce")

if(NOT DEFINED BUILD_DIR)
  message(FATAL_ERROR "ERROR: Must pass in correct \"-DBUILD_DIR=<directory name>\"")
endif()
if(NOT DEFINED BUILD_NAME)
  message(FATAL_ERROR "ERROR: Must pass in \"-DBUILD_NAME=<dashboard build name>\"")
endif()

set(CTEST_BUILD_NAME ${BUILD_NAME})

# default is to build without MPI
if(NOT DEFINED BUILD_WITH_MPI)
  set(BUILD_WITH_MPI "NO")
endif()

# used for invocation of parallel make
if(DEFINED ENV{JOBWEIGHT})
  set(CTEST_BUILD_FLAGS "-j$ENV{JOBWEIGHT}")
else()
  set(CTEST_BUILD_FLAGS "-j8")
endif()

# this should probably be a variable in the environment or passed in,
# but for now it's hard-coded to 1am
set(CTEST_NIGHTLY_START_TIME "01:00:00 MDT")

set(CTEST_SOURCE_DIRECTORY "$ENV{WORKSPACE}/Trilinos")
set(CTEST_BINARY_DIRECTORY "${BUILD_DIR}")

# use the "hostname" command as the CTEST_SITE variable, which is used
# in the "Site" column on the dashboard
find_program(HNAME NAMES hostname)
execute_process(COMMAND "${HNAME}"
  OUTPUT_VARIABLE CTEST_SITE
  OUTPUT_STRIP_TRAILING_WHITESPACE)

# make the hostname generic for all ascic[0-9]* platforms
if(${CTEST_SITE} MATCHES "^ascic[0-9]*")
  set(CTEST_SITE "ascic")
elseif(${CTEST_SITE} MATCHES "cee-build[0-9]*")
  set(CTEST_SITE "cee-build")
endif()

if(VERBOSITY GREATER 4)
  message("[VERB5]: ENV{SNLSYSTEM} = $ENV{SNLSYSTEM}")
endif()
if($ENV{SNLSYSTEM} MATCHES "^cts*")
  set(CTEST_SITE "cts")
endif()

if(NOT DEFINED SITE_APPEND)
  set(CTEST_SITE "${CTEST_SITE} TrilinosDev")
else()
  set(CTEST_SITE "${CTEST_SITE} ${SITE_APPEND}")
endif()
if(VERBOSITY GREATER 0)
  message("[VERB]: CTEST_SITE = ${CTEST_SITE}")
endif()

# The following sets variables as required for the trilinos cmake
# invocation. "ON" vars mean that the cmake invocation for trilinos will
# include them via the command line like "-DVAR_NAME=ON" and similarly
# the "OFF" vars mean that the cmake invocation for trilinos will
# include them via the command line like "-DVAR_NAME=OFF". They are
# set in this manner, one by one, in order to faclitate ease of adding
# or removing as necessary.

# create a list of variables to turn ON for trilinos cmake configuration
set(TRI_CMAKE_VARS_ON "")
list(APPEND TRI_CMAKE_VARS_ON "Trilinos_ENABLE_Amesos")
list(APPEND TRI_CMAKE_VARS_ON "Trilinos_ENABLE_Amesos2")
list(APPEND TRI_CMAKE_VARS_ON "Trilinos_ENABLE_AztecOO")
list(APPEND TRI_CMAKE_VARS_ON "Trilinos_ENABLE_Belos")
list(APPEND TRI_CMAKE_VARS_ON "Trilinos_ENABLE_Epetra")
list(APPEND TRI_CMAKE_VARS_ON "Trilinos_ENABLE_EpetraExt")
list(APPEND TRI_CMAKE_VARS_ON "Trilinos_ENABLE_Ifpack")
list(APPEND TRI_CMAKE_VARS_ON "Trilinos_ENABLE_Isorropia")
list(APPEND TRI_CMAKE_VARS_ON "Trilinos_ENABLE_Kokkos")
list(APPEND TRI_CMAKE_VARS_ON "Trilinos_ENABLE_NOX")
list(APPEND TRI_CMAKE_VARS_ON "Trilinos_ENABLE_Sacado")
list(APPEND TRI_CMAKE_VARS_ON "Trilinos_ENABLE_Stokhos")
list(APPEND TRI_CMAKE_VARS_ON "Trilinos_ENABLE_Teuchos")
list(APPEND TRI_CMAKE_VARS_ON "Trilinos_ENABLE_TrilinosCouplings")
list(APPEND TRI_CMAKE_VARS_ON "Trilinos_ENABLE_Zoltan")

list(APPEND TRI_CMAKE_VARS_ON "Trilinos_ENABLE_CXX11")
list(APPEND TRI_CMAKE_VARS_ON "Trilinos_ENABLE_COMPLEX_DOUBLE")

list(APPEND TRI_CMAKE_VARS_ON "TPL_ENABLE_AMD")
list(APPEND TRI_CMAKE_VARS_ON "TPL_ENABLE_BLAS")
list(APPEND TRI_CMAKE_VARS_ON "TPL_ENABLE_LAPACK")

list(APPEND TRI_CMAKE_VARS_ON "NOX_ENABLE_AztecOO")
list(APPEND TRI_CMAKE_VARS_ON "NOX_ENABLE_Epetra")
list(APPEND TRI_CMAKE_VARS_ON "NOX_ENABLE_Ifpack")

list(APPEND TRI_CMAKE_VARS_ON "Epetra_ENABLE_MPI")

list(APPEND TRI_CMAKE_VARS_ON "EpetraExt_BUILD_BTF")
list(APPEND TRI_CMAKE_VARS_ON "EpetraExt_BUILD_EXPERIMENTAL")
list(APPEND TRI_CMAKE_VARS_ON "EpetraExt_BUILD_GRAPH_REORDERINGS")

list(APPEND TRI_CMAKE_VARS_ON "TPL_ENABLE_AMD")
list(APPEND TRI_CMAKE_VARS_ON "TPL_ENABLE_BLAS")
list(APPEND TRI_CMAKE_VARS_ON "TPL_ENABLE_LAPACK")

list(APPEND TRI_CMAKE_VARS_ON "Amesos_ENABLE_KLU")
list(APPEND TRI_CMAKE_VARS_ON "Amesos2_ENABLE_Basker")
list(APPEND TRI_CMAKE_VARS_ON "Amesos2_ENABLE_KLU2")

# if the NKL library is loaded set things appropriately
set(TRI_CMAKE_VARS_SET "")
if(VERBOSITY GREATER 2)
  message("[VERB]: MKLROOT is $ENV{MKLROOT}")
endif()
if(DEFINED ENV{MKLROOT})
  list(APPEND TRI_CMAKE_VARS_ON "TPL_ENABLE_MKL")

  list(APPEND TRI_CMAKE_VARS_SET "MKL_LIBRARY_DIRS=$ENV{MKLROOT}/lib")
  list(APPEND TRI_CMAKE_VARS_SET "MKL_INCLUDE_DIRS=$ENV{MKLROOT}/include")

  list(APPEND TRI_CMAKE_VARS_SET "BLAS_LIBRARY_NAMES=mkl_rt")
  list(APPEND TRI_CMAKE_VARS_SET "BLAS_LIBRARY_DIRS=$ENV{MKLROOT}/lib")

  list(APPEND TRI_CMAKE_VARS_SET "LAPACK_LIBRARY_NAMES=mkl_rt")
  list(APPEND TRI_CMAKE_VARS_SET "LAPACK_LIBRARY_DIRS=$ENV{MKLROOT}/lib")
endif()

# create a list of variables to turn OFF for trilinos cmake configuration
set(TRI_CMAKE_VARS_OFF  "Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES")
list(APPEND TRI_CMAKE_VARS_OFF "TPL_ENABLE_DLlib")
list(APPEND TRI_CMAKE_VARS_OFF "Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES")

if(VERBOSITY GREATER 3)
  message("[VERB]: TRI_CMAKE_VARS_ON = ${TRI_CMAKE_VARS_ON}")
  message("[VERB]: TRI_CMAKE_VARS_OFF = ${TRI_CMAKE_VARS_OFF}")
  message("[VERB]: TRI_CMAKE_VARS_SET = ${TRI_CMAKE_VARS_SET}")
endif()

# now construct command line options for cmake
set(CMAKE_COMMAND_OPTS "-DTPL_ENABLE_MPI=${BUILD_WITH_MPI}")
foreach(VARNAME ${TRI_CMAKE_VARS_ON})
  set(CMAKE_COMMAND_OPTS "${CMAKE_COMMAND_OPTS} -D${VARNAME}=ON")
endforeach()
foreach(VARNAME ${TRI_CMAKE_VARS_OFF})
  set(CMAKE_COMMAND_OPTS "${CMAKE_COMMAND_OPTS} -D${VARNAME}=OFF")
endforeach()
foreach(VARNAME ${TRI_CMAKE_VARS_SET})
  set(CMAKE_COMMAND_OPTS "${CMAKE_COMMAND_OPTS} -D${VARNAME}")
endforeach()

# add arguments passed in via ctest invocation
foreach(VARARG ${CMAKE_ARGS_LIST})
  set(CMAKE_COMMAND_OPTS "${CMAKE_COMMAND_OPTS} ${VARARG}")
endforeach()

# add the toolkit version
set(CMAKE_COMMAND_OPTS "${CMAKE_COMMAND_OPTS} -DCMAKE_OSX_DEPLOYMENT_TARGET=${OSX_TOOLKIT_VERSION}")

if(VERBOSITY GREATER 2)
  message("[VERB]: CMAKE_COMMAND_OPTS = ${CMAKE_COMMAND_OPTS}")
endif()

# cmake command to configure trilinos
find_program(CTEST_CMAKE_COMMAND cmake REQUIRED)
set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
set(CTEST_CONFIGURE_COMMAND "${CTEST_CMAKE_COMMAND} \"-GUnix Makefiles\"")
set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} ${CMAKE_COMMAND_OPTS} ${CTEST_SOURCE_DIRECTORY}")
if(VERBOSITY GREATER 1)
  message("[VERB]: CTEST_CONFIGURE_COMMAND = ${CTEST_CONFIGURE_COMMAND}")
endif()

# This needs to be set because of CTestConfig.cmake in Trilinos
SET(CMAKE_MODULE_PATH
  "$ENV{WORKSPACE}/Trilinos/cmake/tribits/core/utils")

ctest_start(Nightly GROUP ${TRACK})
ctest_configure()
ctest_build()
ctest_submit(RETRY_COUNT 10 RETRY_DELAY 30)
