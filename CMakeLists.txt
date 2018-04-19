#
# CMakeLists.txt for Xyce
#
# Build Xyce!
#
# Author: Aaron Gibson (asgibso@sandia.gov)
#
# The minimum required version is currently 2.8.8. We need to support RHEL6.
# If you have the chance, it is better to use a newer version of CMake.
#
cmake_minimum_required( VERSION 2.8.8 )

#-------------------------------------------------------
# Define Some Functions and Macros
#-------------------------------------------------------
#
# Define a debugging macro that is useful for debugging CMake.
# This will print all of the variables in CMake at the time of the call.
function( debug_print )
	get_cmake_property(_variableNames VARIABLES)
	list (SORT _variableNames)
	foreach (_variableName ${_variableNames})
		message(STATUS "${_variableName}=${${_variableName}}")
	endforeach()
endfunction()

#
# option_message( option_var_name option_name )
#
# This function prints a standardized message that displays the contents of the
# given option, as stored in the parameter "option_var_name". The option displayed
# should be a simple enable/disable option.
#
# The inputs are:
# - _option_var_name -- The name of the CMake variable storing the option.
# - _option_name -- The human-readable name of the option.
#
function( option_message _option_var_name _option_name )
	set( VAL_TO_PRINT "${${_option_var_name}}" )
	if( NOT ${_option_var_name} )
		set( VAL_TO_PRINT OFF )
	endif()
	message("-- Enable ${_option_name} (${_option_var_name}): ${VAL_TO_PRINT}" )
endfunction()

#-------------------------------------------------------
# Find Trilinos and Compilers
#-------------------------------------------------------
#
# Be nice... Users are likely to set Trilinos_DIR to the root install directory
# and not the directory that explicitly contains the "TrilinosConfig.cmake" file,
# which is what is expected by find_package(). We try and process that here.
if( Trilinos_DIR )
	if( NOT EXISTS "${Trilinos_DIR}/TrilinosConfig.cmake" )
		if( EXISTS "${Trilinos_DIR}/lib/cmake/Trilinos" )
			# We only set this field if the TrilinosConfig.cmake doesn't exist in the 
			# first place, so we don't need to worry too much about overwriting what the
			# user may have explicitly set prior in the cache.
			set( Trilinos_DIR "${Trilinos_DIR}/lib/cmake/Trilinos" CACHE PATH
				"Path to the directory containing TrilinosConfig.cmake." FORCE
			)
		endif()
	endif()
endif()

find_package( Trilinos CONFIG REQUIRED )

option( USE_TRILINOS_COMPILERS
	"Enable to use all compiler settings from the Trilinos Install. NOTES: This only makes sense if you are building Trilinos on the same machine as Xyce!" FALSE )

if( USE_TRILINOS_COMPILERS )
	message( "We will determine compiler information from the existing Trilinos install." )

	set( CMAKE_C_COMPILER ${Trilinos_C_COMPILER} )
	set( CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER} )
	set( CMAKE_Fortran_COMPILER ${Trilinos_Fortran_COMPILER} )

	set( MPI_C_COMPILER ${Trilinos_MPI_C_COMPILER} )
	set( MPI_CXX_COMPILER ${Trilinos_MPI_CXX_COMPILER} )
	set( MPI_Fortran_COMPILER ${Trilinos_MPI_Fortran_COMPILER} )
	
	# Also inherit Trilinos compiler flags
	set( CMAKE_C_FLAGS ${Trilinos_C_COMPILER_FLAGS} )
	set( CMAKE_CXX_FLAGS ${Trilinos_CXX_COMPILER_FLAGS} )
	set( CMAKE_Fortran_COMPILER_FLAGS ${Trilinos_Fortran_COMPILER_FLAGS} )
endif()


#-------------------------------------------------------
# CMake Options and Includes
#-------------------------------------------------------
#
# First, define the Xyce project, which will initialize the compiler information.
project( Xyce )

# Define some macros here that are used throughout the Xyce build.
# Note that we use CMakeParseArguments, available in CMake 2.8.3 and newer.
include( CMakeParseArguments )
include( CheckIncludeFileCXX )
include( CheckCXXSymbolExists )
include( CheckCXXSourceCompiles )
include( CheckCXXSourceRuns )
include( CheckCCompilerFlag )
include( CheckCXXCompilerFlag )

#
# check_cxx_source_file_compiles
#
# This function implements a convenience wrapper around the trycompile stuff in CMake.
# This adds the ability to set any include directories and external libraries to link
# against in the build process, since some of the test programs may require additional
# stuff. This function requires the following arguments:
# * FILE         -- The file to try to compile and link.
# * VAR_TO_STORE -- The CMake variable name to store the results to.
#
# This function also accepts the following optional arguments:
# * EXTERNAL_INCLUDE_DIRS -- A list of include directories to add when compiling.
# * EXTERNAL_LIBS -- A list of libraries to try and link against when linking.
#
# A sample usage of this function is shown below:
# check_cxx_source_file_compiles( FILE "${CMAKE_CURRENT_SOURCE_DIR}/CheckRandom.cxx"
#     VAR_TO_STORE HAS_RANDOM
#     INCLUDE_DIRS
#         ${Boost_INCLUDE_DIRS}
#         ${Trilinos_DIRS}
#     EXTERNAL_LIBS
#         ${Boost_LIBRARIES}
#         ${Trilinos_LIBRARIES}
# )
#
# The sample will try and compile the file "CheckRandom.cxx", including and linking 
# against Boost and Trilinos. The result will be stored in the variable HAS_RANDOM.
#
function( check_cxx_source_file_compiles )
	set( _OPTIONS_ARGS )
	set( _ONE_VALUE_ARGS VAR_TO_STORE FILE )
	set( _MULTI_VALUE_ARGS EXTERNAL_LIBS INCLUDE_DIRS )

	cmake_parse_arguments( _COMPILETEST "${_OPTIONS_ARGS}" "${_ONE_VALUE_ARGS}" "${_MULTI_VALUE_ARGS}" "${ARGN}")

	# Check the file to compile.
	if( NOT _COMPILETEST_FILE )
		message( FATAL_ERROR "Function requires a FILE argument to compile!" )
	endif()
	
	# Check the variable to store the results to.
	if( NOT _COMPILETEST_VAR_TO_STORE )
		message( FATAL_ERROR "Function requires a VAR_TO_STORE argument to compile!" )
	endif()

	if( _COMPILETEST_INCLUDE_DIRS )
		# Add our include directories.
		set( CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES}
			${_COMPILETEST_INCLUDE_DIRS}
		)
	endif()
	if( _COMPILETEST_EXTERNAL_LIBS )
		# Add our external libraries
		set( CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES}
			${_COMPILETEST_EXTERNAL_LIBS}
		)
	endif()
	
	# Read in the file to our variable.
	file( READ "${_COMPILETEST_FILE}" TMPSRC )
	
	# Try and compile (and link!) it.
	check_cxx_source_compiles( "${TMPSRC}" ${_COMPILETEST_VAR_TO_STORE} )
endfunction()

# Set the main directory for the source of our custom platform check programs.
set( CUSTOM_CMAKE_SCRIPT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules" )
set( PLATFORM_CHECK_SRC_DIR  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/platformCheck" )

# Add our custom CMake scripts to the module search path.
list( APPEND CMAKE_MODULE_PATH "${CUSTOM_CMAKE_SCRIPT_DIR}" )

# Store a variable with our linking libraries for Xyce.
list( APPEND Xyce_LINK_LIBS )


# Get some relevant system information.
include( ProcessorCount )
ProcessorCount( BUILD_SYSTEM_PROCESSOR_COUNT )

# Check whether we are building for 64-bit.
set( Xyce_BUILD_64BIT TRUE )
if( NOT CMAKE_SIZEOF_VOID_P EQUAL 8 )
	set( Xyce_BUILD_64BIT FALSE )
endif()

#-------------------------------------------------------
# Handle Trilinos Details
#-------------------------------------------------------
include_directories( ${Trilinos_INCLUDE_DIRS} )
link_directories( ${Trilinos_LIBRARY_DIRS} )

# Check if Trilinos is compiled to have shared libraries, or if the static
# libraries were compiled in a manner to be used in shared libraries (i.e.
# with the -fPIC flag or equivalent). 
set( SHOULD_BUILD_SHARED_LIBS TRUE )
if( NOT Trilinos_BUILD_SHARED_LIBS )
	if( NOT "${Trilinos_CXX_COMPILER_FLAGS}" MATCHES "-f(PIC)|(pic)" )
		set( SHOULD_BUILD_SHARED_LIBS FALSE )
	endif()
endif()
set( Xyce_BUILD_SHARED ${SHOULD_BUILD_SHARED_LIBS} CACHE BOOL "Build the Xyce Shared libraries." )
set( Xyce_BUILD_STATIC TRUE CACHE BOOL "Build Xyce Static libraries." )

# If we aren't building the shared libraries, then we can't link against it either. So
# recalculate the Xyce_LINK_EXE_WITH_SHARED option, even if it is already in the cache.
if( NOT Xyce_BUILD_SHARED )
	unset( Xyce_LINK_EXE_WITH_SHARED CACHE )
endif()
set( Xyce_LINK_EXE_WITH_SHARED ${Xyce_BUILD_SHARED} CACHE BOOL "Link Xyce executable with shared library." )

set( Trilinos_USES_INTEL FALSE )
if( "${Trilinos_CXX_COMPILER_ID}" STREQUAL "Intel" )
	set( Trilinos_USES_INTEL TRUE )
endif()

set( Xyce_USES_INTEL FALSE )
if( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel" )
	set( Xyce_USES_INTEL TRUE )
endif()

# Print a warning if the compiler IDs for Xyce and Trilinos do not match.
# This warns the user in the event that Trilinos is using the Intel compiler
# (for example), while the Xyce build is using GCC or Clang.
if( NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "${Trilinos_CXX_COMPILER_ID}" )
	message( WARNING
		"WARNING!!! Current Compiler ID (${CMAKE_CXX_COMPILER_ID}) does NOT match "
		"Trilinos compiler ID (${Trilinos_CXX_COMPILER_ID})! This may result in "
		"subtle build errors later on! Do this only at your own risk."
	)
endif()

# Check if Trilinos is compiled to use OpenMPI. If so, we try and "reuse" those
# settings here for consistency.
set( SHOULD_USE_MPI FALSE )
if( Trilinos_MPI_EXEC )
	set( SHOULD_USE_MPI TRUE )
endif()
set( Xyce_PARALLEL_MPI ${SHOULD_USE_MPI} CACHE BOOL "Build Xyce to use MPI." )

# Handle the additional details for a parallel build.
#
# This may require setting custom compiler flags and so on as well as adding additional
# include paths and libraries depending on compiler. We adjust all of that here, and 
# also sanity-check that the compiler paths make sense.
if( Xyce_PARALLEL_MPI )
	# Note that the caller may have explicitly set the MPI_<lang>_COMPILER
	# flags as a hint already. We only default to use what Trilinos is using
	# if these parameters weren't explicitly set.
	if( NOT MPI_C_COMPILER OR NOT MPI_CXX_COMPILER )
		# See if our current compiler can be made to work with OpenMPI.
		set( MPI_C_COMPILER "${CMAKE_C_COMPILER}" )
		set( MPI_CXX_COMPILER "${CMAKE_CXX_COMPILER}" )
	endif()
	#if( NOT MPI_Fortran_COMPILER )
	#	set( MPI_Fortran_COMPILER "${Trilinos_Fortran_COMPILER}" )
	#endif()
	
	# Try to find MPI with our current settings. 
	find_package( MPI )
		
	# If we found our libraries, set that we should build with MPI.
	if( NOT MPI_C_FOUND )
		message( FATAL_ERROR "Xyce is expected to build with MPI, but no MPI C compiler found! "
		"Please specify this directly via: CMAKE_C_COMPILER." )
	endif()
	if( NOT MPI_CXX_FOUND )
		message( FATAL_ERROR "Xyce is expected to build with MPI, but no MPI C++ compiler found! "
		"Please specify this directly via: CMAKE_CXX_COMPILER." )
	endif()
	
	include_directories( SYSTEM ${MPI_CXX_INCLUDE_DIRS} )
	list( APPEND Xyce_LINK_LIBS ${MPI_CXX_LIBRARIES} )
endif()

# Check if Trilinos is compiled with C++11 enabled. We need to build Xyce with
# the same standard for C++, so check that before doing anything else.
#
# Note that we assume C++11 is enabled by default, because this is what is supported
# by most compilers. We will only disable it if Trilinos is explicitly configured with
# C++11 set false. 
#
# Note that there are a few corner-cases that this may not cover and which cannot be
# easily detected, such as gcc 4.4, where C++11 is disabled, but no flags are passed.
# Since Trilinos doesn't store the compiler version, this is hard to detect, so the
# users may need to explicitly set C++11 false themselves via a cmake flag:
# .. -DXyce_ENABLE_CXX11=FALSE
#
# Also, analyzing compiler flags is not sufficient because GCC 6.0 enables C++14 by
# default with no flags and again, Trilinos doesn't store the compiler version.
#
# However, this flow should handle the most common cases, and frankly, we should
# just build the default with C++11 anyway.
#
set( SHOULD_USE_CXX11 TRUE )
if( DEFINED Trilinos_ENABLE_CXX11 )
	if( NOT Trilinos_ENABLE_CXX11 )
		set( SHOULD_USE_CXX11 FALSE )
		message( "-- Trilinos installation has explicitly disabled C++11." )
	endif()
endif()
set( Xyce_ENABLE_CXX11 ${SHOULD_USE_CXX11} CACHE BOOL "Enable C++11 for Xyce" )

if( Xyce_ENABLE_CXX11 )
	set( CMAKE_CXX_STANDARD 11 )
else()
	# Note that we explicitly set this because newer compilers may need to
	# explicitly disable C++11 if they compile to C++11 by default.
	set( CMAKE_CXX_STANDARD 98 )
endif()

# Check for our threading libraries.
set( THREADS_PREFER_PTHREAD_FLAG ON )
find_package( Threads )
list( APPEND Xyce_LINK_LIBS Threads::Threads )

message( "System Details:" )
message( "-- OS Info: ${CMAKE_SYSTEM}" )
message( "-- Procesor Type: ${CMAKE_SYSTEM_PROCESSOR}" )
message( "-- Processor Count: ${BUILD_SYSTEM_PROCESSOR_COUNT}" )
message( "-- 64 bit Build: ${Xyce_BUILD_64BIT}" )
message( " " )
message( "Build Details:" )
option_message( Xyce_ENABLE_CXX11 "C++11" )
option_message( Xyce_BUILD_SHARED "Shared libraries" )
option_message( Xyce_BUILD_STATIC "Static Libraries" )
option_message( Xyce_LINK_EXE_WITH_SHARED "Link Xyce with Shared" )
option_message( Xyce_PARALLEL_MPI "MPI" )
message( " " )
message( "Compiler Details:" )
message( "-- C Compiler  : ${CMAKE_C_COMPILER}" )
message( "-- CXX Compiler: ${CMAKE_CXX_COMPILER}" )
message( " " )


#-------------------------------------------------------
# Find Other Packages (Non-Trilinos)
#-------------------------------------------------------
#
# OpenMP
#
find_package( OpenMP )

#-------------------------------------------------------
# Flex and Bison
find_package( FLEX )
find_package( BISON )

#-------------------------------------------------------
# Git
#
# We use git to calculate the XYCEBUILDTIMESTAMP compiler define.
find_package( Git )

#-------------------------------------------------------
# Dakota
#
# Be nice... Users are likely to set Dakota_DIR to the root install directory
# and not the directory that explicitly contains the "DakotaConfig.cmake" file,
# which is what is expected by find_package(). We try and process that here.
if( Dakota_DIR )
        if( NOT EXISTS "${Dakota_DIR}/DakotaConfig.cmake" )
                if( EXISTS "${Dakota_DIR}/lib/cmake/Dakota" )
                        # We only set this field if the DakotaConfig.cmake doesn't exist in the 
                        # first place, so we don't need to worry too much about overwriting what the
                        # user may have explicitly set prior in the cache.
                        set( Dakota_DIR "${Dakota_DIR}/lib/cmake/Dakota" CACHE PATH
                                "Path to the directory containing DakotaConfig.cmake." FORCE
                        )
                endif()
        endif()
endif()

# Search for Dakota
find_package( Dakota CONFIG QUIET )

#-------------------------------------------------------
# Boost
#
# Some models might want to use boost, so we search for that as well.
#
find_package( Boost COMPONENTS system filesystem serialization signals regex QUIET )

set( SHOULD_USE_BOOST FALSE )
if( Boost_FOUND )
	set( SHOULD_USE_BOOST TRUE )
endif()
set( Xyce_USE_BOOST ${SHOULD_USE_BOOST} CACHE BOOL "Use Boost Libraries, potentially needed by some models." )

if( Xyce_USE_BOOST )
	if( NOT Boost_FOUND )
		message( SEND_ERROR "Xyce_USE_BOOST enabled, but Boost was not found!" )
	endif()
	include_directories( ${Boost_INCLUDE_DIRS} )

	list( APPEND Xyce_LINK_LIBS
		${Boost_LIBRARIES}
#		Boost::system
#		Boost::filesystem
#		Boost::serialization
#		Boost::signals
#		Boost::regex
	)
endif()

#-------------------------------------------------------
# Intel Math Kernel Libraries (MKL)
#
# Search for the Intel MKL Libraries. Note that it is technically possible
# to use some versions of the MKL libraries with a GCC compiler, so we allow
# it, but we will only search for the MKL stuff automatically for the Intel
# compilers.
#
# The Intel MKL stuff should be searched for using the "MKLROOT" variable.
# We check if it is set in our current environment, or otherwise search our
# existing (Intel) compiler path. This transparently handles cases where the
# Intel compilers we are using happen to be in non-standard locations.
#
# https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor
#

# Set the default search path for MKLROOT
if( DEFINED $ENV{MKLROOT} )
	set( MKLROOT_PATH_VAR "$ENV{MKLROOT}" )
elseif( Xyce_USES_INTEL )
	get_filename_component( INTEL_BASE_PATH ${CMAKE_CXX_COMPILER} DIRECTORY )
	set( MKLROOT_PATH_VAR "${INTEL_BASE_PATH}/../mkl" )
else()
	set( MKLROOT_PATH_VAR "MKL-NOTFOUND" )
endif()
set( MKLROOT "${MKLROOT_PATH_VAR}" CACHE PATH "Path to the Intel Math Kernel Libraries (MKL)" )

# Search for Intel Math Kernel Library stuff.
find_package( MKL )

# Find the FFTW library
find_package( FFTW )

#-------------------------------------------------------
# Check for CURL Stuff 
find_package( CURL )


#-------------------------------------------------------
# Print Third-Party Library Results
#-------------------------------------------------------
message( "Third Party Libraries:" )
message( "-- Trilinos Include Dir: ${Trilinos_INCLUDE_DIRS}" )
if( OPENMP_FOUND OR OpenMP_FOUND )
	message( "-- Found OpenMP (${OpenMP_VERSION})" )
	message( "-- * C   Flags: ${OpenMP_C_FLAGS}" )
	message( "-- * CXX Flags: ${OpenMP_CXX_FLAGS}" )
else()
	message( "-- Could not find OpenMP for this compiler." )
endif()

if( FLEX_FOUND )
	message( "-- Found Flex (${FLEX_VERSION}): ${FLEX_EXECUTABLE}" )
else()
	message( "-- Could not find Flex." )
endif()
if( BISON_FOUND )
	message( "-- Found Bison (${BISON_VERSION}): ${BISON_EXECUTABLE}" )
else()
	message( "-- Could not find Bison." ) 
endif()
if( Git_FOUND OR GIT_FOUND )
	message( "-- Found Git (${GIT_VERSION_STRING})" )
	message( "-- * Path: ${GIT_EXECUTABLE}" )
else()
	message( "-- Could not find Git. Set the path to git via: -DGIT_EXECUTABLE=.." )
endif()
if( Boost_FOUND )
	message( "-- Found Boost (${Boost_VERSION})" )
	message( "-- * Include Dir: ${Boost_INCLUDE_DIR}" )
else()
	message( "-- Could not find Boost." )
endif()
if( MKL_FOUND )
	message( "-- Found Intel Math Kernel Libraries (MKL)" )
	message( "-- * MKLROOT: ${MKLROOT}" )
	message( "-- * Include Dir: ${MKL_INCLUDE_DIRS}" )
else()
	message( "-- Could not find Intel MKL. Set MKLROOT appropriately to find them." )
endif()
if( FFTW_FOUND )
	message( "-- Found FFTW Libraries" )
	message( "-- * Include Dir: ${FFTW_INCLUDE_DIRS}" )
	message( "-- * Libraries: ${FFTW_LIBRARIES}" )
else()
	message( "-- Could not find FFTW Libraries" )
endif()

if( Dakota_FOUND OR DAKOTA_FOUND )
	message( "-- Found Dakota Libraries" )
	message( "-- * Include Dir: ${Dakota_INCLUDE_DIRS}" )
	message( "-- * Executable : ${Dakota_EXECUTABLE}" )
else()
	message( "-- Could not find Dakota. Set Dakota_DIR appropriately to find it." )
endif()

if( CURL_FOUND )
	message( "-- Found cURL (${CURL_VERSION_STRING})" )
	message( "-- * Include Dir: ${CURL_INCLUDE_DIRS}" )
	message( "-- * Libraries: ${CURL_LIBRARIES}" )
else()
	message( "-- Could not find cURL libraries." )
endif()

# For formatting
message( " " )

#-------------------------------------------------------
# CMake checks
#-------------------------------------------------------
#
# Before running all of these checks, we should set the compiler flags
# to be what they would be when we perform our full compile. Unfortunately,
# CMake doesn't do this by default, because it is hypothetically possible
# for some projects that the flags used are related to the results of the
# check.
# In our case, we just need to know if we are invoking the compiler in a
# C++11 setting or not and add that flag. We know the version we are using
# based on what Trilinos used, since that is the massive library that we
# are linking against.
set( CMAKE_REQUIRED_FLAGS "${CMAKE_CXX_FLAGS}" )
if( Xyce_ENABLE_CXX11 )
	set( CMAKE_REQUIRED_DEFINITIONS ${CMAKE_CXX11_STANDARD_COMPILE_OPTION} )
endif()

# Check for C++11 features like: <unordered_map>, <unordered_set>, std::iota
if( Xyce_ENABLE_CXX11 )
	# Only check for these headers if we are building with C++11. This is
	# necessary to avoid legacy compilers that might have C++11 headers
	# in the standard include path, even though the current compiler can't
	# handle them itself. We treat as if we never found them in that case.
	check_include_file_cxx( "unordered_map" HAVE_UNORDERED_MAP )
	check_include_file_cxx( "unordered_set" HAVE_UNORDERED_SET )

	# std::iota is (usually) a C++11 feature, but can be offered on other
	# compilers as an extension. We'll only check for it if C++11 is enabled.
	check_cxx_source_file_compiles( FILE "${PLATFORM_CHECK_SRC_DIR}/CheckForIota.cxx"
		VAR_TO_STORE HAVE_IOTA )
else()
	# We need to unset these variables if C++11 isn't enabled, because
	# these variables are in the cache and won't otherwise be re-evaluated
	# if "Xyce_ENABLE_CXX11" is changed from it was originally.
	#
	# Hopefully, this prevents some frustrating issues that can happen if
	# Xyce_ENABLE_CXX11 is set true, but later set false, where unordered_map may not
	# be supported.
	unset( HAVE_UNORDERED_MAP CACHE )
	unset( HAVE_UNORDERED_SET CACHE )
	unset( HAVE_IOTA CACHE )
endif()
check_include_file_cxx( "tr1/unordered_map" HAVE_TR1_UNORDERED_MAP )
check_include_file_cxx( "tr1/unordered_set" HAVE_TR1_UNORDERED_SET )

# Check for C++ common headers.
check_include_file_cxx( "numeric"    HAVE_NUMERIC )
check_include_file_cxx( "algorithm"  HAVE_ALGORITHM )
check_include_file_cxx( "functional" HAVE_FUNCTIONAL )
check_include_file_cxx( "cstdlib"    HAVE_CSTDLIB )
check_include_file_cxx( "iostream"   HAVE_IOSTREAM )
check_include_file_cxx( "fstream"    HAVE_FSTREAM )
check_include_file_cxx( "unistd.h"   HAVE_UNISTD_H )
check_include_file_cxx( "dlfcn.h"    HAVE_DLFCN_H )
check_include_file_cxx( "float.h"    HAVE_FLOAT_H )
check_include_file_cxx( "strings.h"  HAVE_STRINGS_H )
check_include_file_cxx( "string.h"   HAVE_STRING_H )
check_include_file_cxx( "inttypes.h" HAVE_INTTYPES_H ) 
check_include_file_cxx( "limits.h"   HAVE_LIMITS_H )
check_include_file_cxx( "malloc.h"   HAVE_MALLOC_H )
check_include_file_cxx( "math.h"     HAVE_MATH_H )
check_include_file_cxx( "mathimf.h"  HAVE_MATHIMF_H )

# Check for Windoze header.
check_include_file_cxx( "Windows.h" HAVE_WINDOWS_H )

# Check for sys/* headers.
check_include_file_cxx( "sys/resource.h" HAVE_SYS_RESOURCE_H )
check_include_file_cxx( "sys/stat.h"     HAVE_SYS_STAT_H )
check_include_file_cxx( "sys/types.h"    HAVE_SYS_TYPES_H )

# Check for linux exceptions
check_cxx_source_file_compiles( VAR_TO_STORE FOUND_FP_EXCEPTIONS
	FILE "${PLATFORM_CHECK_SRC_DIR}/CheckForLinuxExceptions.cxx" FOUND_FP_EXCEPTIONS
)


#-------------------------------------------------------
# Determine TPL Options (excluding Trilinos)
#-------------------------------------------------------
#
# Determine whether to build the Reaction parser
set( SHOULD_USE_REACTION_PARSER FALSE )
if( FLEX_FOUND AND BISON_FOUND )
	# TODO -- Check if we are using the correct version of Flex and Bison ?
	set( SHOULD_USE_REACTION_PARSER TRUE )
endif()
set( Xyce_REACTION_PARSER ${SHOULD_USE_REACTION_PARSER} CACHE BOOL "Enable the Reaction Parser." )

set( HAVE_LINUX_EXCEPTIONS FALSE CACHE BOOL "Set to true to enable floating point exceptions." )

if( HAVE_LINUX_EXCEPTIONS AND (NOT FOUND_FP_EXCEPTIONS ) )
	# Warn if we explicitly enable HAVE_LINUX_EXCEPTIONS, but don't actually have the
	# appropriate headers.
	# Also, forcibly set HAVE_LINUX_EXCEPTIONS back to false.
	message( WARNING "Floating point exceptions (HAVE_LINUX_EXCEPTIONS) were enabled, "
		"but the appropriate headers were not found. Disabling floating point exceptions."
	)
	set( HAVE_LINUX_EXCEPTIONS FALSE CACHE BOOL
		"Set to true to enable floating point exceptions." FORCE )
		
endif()

#-------------------------------------------------------
# FFT Stuff
#
# We can configure the build to use FFTW or Intel's Math Kernel Libraries. We 
# have various options here that let us configure this stuff appropriately.
#
# Due to the potentially complicated Intel MKL stuff, we add an option to explicitly
# add the Intel MKL Compiler and Linker flags manually. This will override whatever
# is found by the MKL module above, if set. Note that these manual options should NOT
# set the optimization or debug levels (i.e. -O3 -g) because those flags are added by
# CMake directly, based on build type. The only flags that should be added are the
# flags from Intel's MKL link line advisor:
# https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor
#
# This is a workaround when the "FindMKL" module isn't working as desired.
#
# Declare a place to store the FFT libraries that we are going to use.
set( FFT_LIBRARIES_TO_USE )

# Declare some options for setting Intel MKL options manually.
set( Xyce_INTEL_MKL_CXX_FLAGS "" CACHE STRING "Compile flags for Intel MKL." )
set( Xyce_INTEL_MKL_LINK_FLAGS "" CACHE STRING "Linker flags for Intel MKL." )

# Determine whether the options were explicitly set or not.
set( Xyce_MANUAL_MKL FALSE )
if( Xyce_INTEL_MKL_CXX_FLAGS OR Xyce_INTEL_MKL_LINK_FLAGS )
	set( Xyce_MANUAL_MKL TRUE )
endif()

# If we are manually setting the link line for Intel's MKL, we create an interface
# library that stores the resulting options. Then, all we need to do to enable these
# options is link with the created Interface library.
#
# This is the appropriate "CMake" way to add optional compiler and linker flags more
# easily.
#
if( Xyce_MANUAL_MKL )
	# If we are manually setting the MKL flags, then we create an imported target
	# with the appropriate settings, then set them as our FFT_LIBRARIES.
	add_library( MKL_Manual INTERFACE IMPORTED )
	set_property( TARGET MKL_Manual
		PROPERTY INTERFACE_COMPILE_OPTIONS "${Xyce_INTEL_MKL_CXX_FLAGS}"
	)
	set_property( TARGET MKL_Manual
		PROPERTY INTERFACE_LINK_LIBRARIES "${Xyce_INTEL_MKL_LINK_FLAGS}"
	)
endif()

# Check if we can use Intel's FFT.
if( Xyce_MANUAL_MKL )
	set( TEMP_MKL_LIBS MKL_Manual )
else()
	set( TEMP_MKL_LIBS ${MKL_LIBRARIES} )
endif()

check_cxx_source_file_compiles( VAR_TO_STORE FOUND_INTEL_FFT
	FILE "${PLATFORM_CHECK_SRC_DIR}/CheckForIntelFFT.cxx"
	INCLUDE_DIRS
		${MKL_INCLUDE_DIRS}
	EXTERNAL_LIBS
		${TEMP_MKL_LIBS}
)

# Print an error if one FFT type is explicitly requested, but not found.
if( Xyce_USE_FFTW AND NOT FFTW_FOUND )
	message( SEND_ERROR "Xyce_USE_FFTW enabled, but could not find FFTW libraries. "
		"Set path to FFTW libraries."
	)
endif()

if( Xyce_USE_INTEL_FFT AND NOT FOUND_INTEL_FFT )
	message( SEND_ERROR "Xyce_USE_INTEL_FFT set, but we could not find or compile with Intel FFT. "
		"Make sure MKLROOT is set, and/or set the Intel MKL link line manually via: "
		"Xyce_INTEL_MKL_LINK_FLAGS and Xyce_INTEL_MKL_CXX_FLAGS"
	)
endif()

# We should enable Intel FFT if we found it, or if we have manually set the 
# MKL compile and linker flags.
set( HAS_INTEL_FFT FALSE )
if( FOUND_INTEL_FFT OR Xyce_MANUAL_MKL )
	set( HAS_INTEL_FFT TRUE )
endif()

# Determine what FFT library we should use in the build, if any.
set( SHOULD_ENABLE_FFT FALSE )
if( FOUND_INTEL_FFT OR FFTW_FOUND )
	set( SHOULD_ENABLE_FFT TRUE )
endif()

set( Xyce_USE_FFT "${SHOULD_ENABLE_FFT}" CACHE BOOL "Enable FFT for Xyce." )

# Now, choose the library, if it is enabled.
if( Xyce_USE_FFT )
	
	# IF FFT is enabled, we have options to choose between Intel FFT or FFTW.
	set( Xyce_USE_INTEL_FFT "${HAS_INTEL_FFT}" CACHE BOOL "Enable the Intel FFT." )
	set( Xyce_USE_FFTW "${FFTW_FOUND}" CACHE BOOL "Enable the FFTW library for FFT." )

	# More logic to sanity-check that only one of Xyce_USE_INTEL_FFT or Xyce_USE_FFTW
	# are used. This handles cases where the user manually sets these flags, bypassing
	# our hint checks above.
	if( Xyce_USE_INTEL_FFT AND Xyce_USE_FFTW )
		message( WARNING "Only one of Xyce_USE_INTEL_FFT or Xyce_USE_FFTW should be enabled at a time! "
			"Disable one of them, or disable FFT entirely with: -DXyce_USE_FFT=FALSE. "
			"For now, we will assume the Intel FFT libraries and disable FFTW explicitly."
		)
		set( Xyce_USE_FFTW FALSE CACHE BOOL "Enable the FFTW library for FFT." FORCE )
	endif()
	
	# Set which FFT libraries we are using to FFT_LIBRARIES_TO_USE
	if( Xyce_USE_INTEL_FFT )
		# If we are using Intel's FFT, include the MKL libraries.
		include_directories( "${MKL_INCLUDE_DIRS}" )
		
		# If we have manually set the compile and link lines, use our imported target
		if( Xyce_MANUAL_MKL )
			set( FFT_LIBRARIES_TO_USE MKL_Manual )
		else()
			set( FFT_LIBRARIES_TO_USE "${MKL_LIBRARIES}" )
		endif()
	elseif( Xyce_USE_FFTW )
		# If we are using FFTW, make sure to add them to the include path and our list
		# of libraries to link against in the final targets.
		include_directories( ${FFTW_INCLUDE_DIRS} )
		
		set( FFT_LIBRARIES_TO_USE ${FFTW_LIBRARIES} )
	endif()
else()
	# If we disable FFT entirely, we should unset the other options too. That way, 
	# the options get reset if re-enabled.
	# This also prevents these from getting set in "Xyce_config.h" accidentally.
	#
	unset( Xyce_USE_INTEL_FFT CACHE )
	unset( Xyce_USE_FFTW CACHE )
	
	# Also unset our option from "trycompile".
	unset( FOUND_INTEL_FFT CACHE )
endif()

# Finally, add our FFT libraries to the final link line.
if( FFT_LIBRARIES_TO_USE )
	list( APPEND Xyce_LINK_LIBS ${FFT_LIBRARIES_TO_USE} )
endif()

#-------------------------------------------------------
# Dakota Stuff
#

set( SHOULD_USE_DAKOTA FALSE )
if( Dakota_FOUND OR DAKOTA_FOUND )
	if( Xyce_USE_BOOST )
		# Compile a sample program to see if we can actually link against Dakota with our
		# current libraries. Remember that Dakota might require Boost, and we need to make
		# sure that we link properly against it.
		check_cxx_source_file_compiles( VAR_TO_STORE DAKOTA_COMPILES
			FILE "${PLATFORM_CHECK_SRC_DIR}/CheckDakota.cxx"
			INCLUDE_DIRS
				${Dakota_INCLUDE_DIRS}
			EXTERNAL_LIBS
				${Xyce_LINK_LIBS}
				${Dakota_LIBRARIES}
				${Boost_LIBRARIES}
		)

		if( DAKOTA_COMPILES )			
			set( SHOULD_USE_DAKOTA TRUE )
		else()
			message( WARNING "Dakota and Boost were both found, but the test program "
				"did not compile. It is likely that the found Boost installation is not "
				"compatible with the found Dakota installation. Try setting BOOST_ROOT "
				"or Dakota_DIR to versions that work."
			)
		endif()
	endif()
endif()

# TODO -- Should we try and compile a sample using dakota?
# By default, we already check for some of the dependencies, but it might be better
# to have a complete compile to be absolutely sure.
set( Xyce_Dakota "${SHOULD_USE_DAKOTA}" CACHE BOOL "Build with Dakota support." )

# In the meantime, we will warn if it does not immediately appear that Dakota's
# dependencies are appropriately satisfied.
if( Xyce_Dakota )
	if( NOT Xyce_USE_BOOST )
		message( SEND_ERROR "Dakota was enabled, but Boost (a required dependency of Dakota) was not. "
				"Disable Dakota: -DXyce_Dakota=OFF, or enable boost: -DXyce_USE_BOOST=ON"
		)
	endif()
	
	# Include Dakota directories
	include_directories( ${Dakota_INCLUDE_DIRS} )
	
	# Include Dakota libraries.
	list( APPEND Xyce_LINK_LIBS ${Dakota_LIBRARIES} )
	list( APPEND Xyce_LINK_LIBS Boost::filesystem )
endif()

#-------------------------------------------------------
# Tracking URL Stuff
#
# Handle the tracking URL here.
#
set( ENABLE_CURL_TRACKING FALSE )
if( CURL_FOUND )
	if( NOT DEFINED ENV{XYCE_NO_TRACKING} )
		if( DEFINED TRACKING_URL )
			set( ENABLE_CURL_TRACKING TRUE )
		elseif( DEFINED ENV{TRACKING_URL} )
			set( TRACKING_URL "$ENV{TRACKING_URL}" )
			set( ENABLE_CURL_TRACKING TRUE )
		endif()
	endif()
endif()
if( NOT ENABLE_CURL_TRACKING )
	set( TRACKING_URL "" )
endif()

set( Xyce_TRACKING_URL "${TRACKING_URL}" CACHE STRING "Xyce Usage Tracking URL" )
set( Xyce_USE_CURL     ${ENABLE_CURL_TRACKING} CACHE BOOL "Enable Xyce usage tracking." )

#---------------------------------
# Add Generic (Non-Trilinos) TPL Libraries
#
# Include the dynamic library loader libraries, if applicable.
if( HAVE_DLFCN_H )
	list( APPEND Xyce_LINK_LIBS ${CMAKE_DL_LIBS} )
endif()

#---------------------------------
# Print our found options
#
message( "Non-Trilinos Options:" )
option_message( HAVE_LINUX_EXCEPTIONS "Enable floating point exceptions." )
option_message( Xyce_REACTION_PARSER "Reaction Parser" )
option_message( Xyce_USE_FFT "FFT" )
if( Xyce_USE_FFT )
	message( "-- * Using Intel FFT (Xyce_USE_INTEL_FFT): ${Xyce_USE_INTEL_FFT}" )
	message( "-- * > Manual CXX  Flags (Xyce_INTEL_MKL_CXX_FLAGS) : ${Xyce_INTEL_MKL_CXX_FLAGS}" )
	message( "-- * > Manual Link Flags (Xyce_INTEL_MKL_LINK_FLAGS): ${Xyce_INTEL_MKL_LINK_FLAGS}" )
	message( "-- * Using FFTW (Xyce_USE_FFTW): ${Xyce_USE_FFTW}" )
else()
	message( "-- * Set Xyce_USE_FFTW=TRUE or Xyce_USE_INTEL_FFT=TRUE to enable." )
	message( "-- * Note that the Intel MKL flags can be set explicitly via:" )
	message( "-- * > Xyce_INTEL_MKL_CXX_FLAGS" )
	message( "-- * > Xyce_INTEL_MKL_LINK_FLAGS" )
endif()

option_message( Xyce_USE_BOOST "Boost libraries." )
option_message( Xyce_Dakota "Dakota" )
message( "-- Enable Usage Tracking: ${Xyce_USE_CURL}" )
message( "-- * Tracking URL: ${Xyce_TRACKING_URL}" )

# For formatting
message( " " )

#-------------------------------------------------------
# Trilinos Checks
#-------------------------------------------------------
#
# Define a special macro to compile source with Trilinos.
macro( check_cxx_source_compiles_with_trilinos COMPILATION_FILE CMAKE_VAR_TO_STORE )
	check_cxx_source_file_compiles( VAR_TO_STORE "${CMAKE_VAR_TO_STORE}"
		FILE "${COMPILATION_FILE}"
		INCLUDE_DIRS
			${Trilinos_INCLUDE_DIRS}
		EXTERNAL_LIBS
			${Xyce_LINK_LIBS}
			${Trilinos_LIBRARIES}
	)
endmacro()

# Define a special macro to search for Trilinos headers.
macro( check_trilinos_header HeaderFileToCheck CMAKE_VAR_TO_STORE )
	# Note that we cannot simply use 'check_include_file_cxx' because
	# CMake tries to link the result, and this macro doesn't let us 
	# specify the libraries to link against.
	# 
	# Trilinos may implicitly define some symbols that we need to link
	# against for it to properly compile, even if we only include the
	# header. This is due to the way that it initializes things. So, we
	# instead invoke our macro:
	# check_cxx_source_compiles_with_trilinos() to go ahead and do the
	# full compile, but in a way that let's us include the Trilinos 
	# libraries during the linking stage.
	#
	# The code below borrows heavily from what the CheckIncludeFileCXX
	# defines anyway:
	set( FILE_TO_COMPILE
		"${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/CheckIncludeFile.cxx" )
	# Setup our configuration file with the actual include to test.
	set(CHECK_INCLUDE_FILE_VAR ${HeaderFileToCheck})
	configure_file( ${PLATFORM_CHECK_SRC_DIR}/CheckIncludeFile.cxx.in
		"${FILE_TO_COMPILE}" )

	check_cxx_source_compiles_with_trilinos( "${FILE_TO_COMPILE}" ${CMAKE_VAR_TO_STORE} )
endmacro()

if( NOT Trilinos_Teuchos_FOUND )
	message( FATAL_ERROR "Xyce requires the Teuchos package of Trilinos.\n"
		"Try setting: Trilinos_ENABLE_Teuchos=ON in your Trilinos build."
	)
endif()

check_cxx_source_compiles_with_trilinos(
	"${PLATFORM_CHECK_SRC_DIR}/CheckTeuchosComplex.cxx" HAS_TEUCHOS_COMPLEX
)
if( NOT HAS_TEUCHOS_COMPLEX )
	# Remember to unset this variable, so we run the check again next time.
	# We do this here explicitly, because this is a required Xyce parameter.
	unset( HAS_TEUCHOS_COMPLEX CACHE )
	message( FATAL_ERROR
		"Error with Trilinos install: Teuchos was not built with COMPLEX enabled.\n"
		"To enable, make sure to set 'Teuchos_ENABLE_COMPLEX=ON' in your Trilinos build."
	)
endif()

# BELOS -- Add this library if we found it.
set( SHOULD_USE_BELOS OFF )
if ( Trilinos_Belos_FOUND )
	set( SHOULD_USE_BELOS ON )
endif()
set( Xyce_BELOS ${SHOULD_USE_BELOS} CACHE BOOL "Use Belos" )

# NOX FOUND -- If NOX stuff was found, enable features.
set( SHOULD_USE_NOX OFF )
if( Trilinos_NOX_FOUND )
	set( SHOULD_USE_NOX ON )
endif()
set( Xyce_NOX_USE_VECTOR_COPY ${SHOULD_USE_NOX} CACHE BOOL "Enable NOX Vector Copy" )

# Check for Amesos libraries
if( Trilinos_Amesos_FOUND )
	check_trilinos_header( "Amesos2.hpp" HAS_AMESOS2_HEADER )
endif()
set( Xyce_AMESOS2 ${HAS_AMESOS2_HEADER} CACHE BOOL "Enable Amesos2 Library." )

# Check for libraries that depend on EPetraExt
if( NOT Trilinos_EpetraExt_FOUND )
	message( FATAL_ERROR "Trilinos requires the EpetraExt libraries to be enabled.\n"
		"Try rebuilding Trilinos with: Trilinos_ENABLE_EpetraExt=ON\n"
		"You will also likely need to set: EpetraExt_BUILD_BTF=ON\n"
		"Epetra_BUILD_EXPERIMENTAL=ON, EpetraExt_BUILD_GRAPH_REORDERINGS=ON"
	)
endif()

# Check for the EpetraExt_BTF_CrsGraph.h header.
check_trilinos_header( "EpetraExt_BTF_CrsGraph.h" HAS_EPETRA_BTF )

# Check for the EpetraExt_AMD_CrsGraph.h header.
check_trilinos_header( "EpetraExt_AMD_CrsGraph.h" HAS_EPETRA_AMD_HEADER )

# Check for the EpetraExt_AmesosBTFGlobal_LinearProblem.h header.
check_trilinos_header( "EpetraExt_AmesosBTFGlobal_LinearProblem.h"
	HAS_EPETRA_AMESOS_BTF_GLOBAL_LP_HEADER )

# Check for the EpetraExt_Isorropia_CrsGraph.h header.
check_trilinos_header( "EpetraExt_Isorropia_CrsGraph.h" HAS_EPETRAEXT_ISOROPPIA )

# If Xyce_PARALLEL_MPI is set, we also need a few additional headers in
# order for Xyce to build:
set( SHOULD_USE_ISORROPIA OFF )
if( HAS_EPETRAEXT_ISOROPPIA AND HAS_EPETRA_AMESOS_BTF_GLOBAL_LP_HEADER )
	set( SHOULD_USE_ISORROPIA ON )
endif()
if( Xyce_PARALLEL_MPI AND NOT SHOULD_USE_ISOROPPIA )
	if( NOT HAS_EPETRAEXT_ISOROPPIA )
		message( SEND_ERROR
		"Parallel builds of Xyce require Trilinos to have the following header\n"
		"EpetraExt_Isorropia_CrsGraph.h, but it was not found!"
		)
	endif()
	if( NOT HAS_EPETRA_AMESOS_BTF_GLOBAL_LP_HEADER )
		message( SEND_ERROR
		"Parallel builds of Xyce require Trilinos to have the following header\n"
		"EpetraExt_AmesosBTFGlobal_LinearProblem.h, but it was not found!"
		)
	endif()
endif()
set( Xyce_USE_ISORROPIA ${SHOULD_USE_ISORROPIA} CACHE BOOL "Enable Isorropia in Xyce. (Needed for Parallel builds!)" )

#	# Check for the ML header.
#	check_trilinos_header( "ml_MultiLevelPreconditioner.h" HAS_ML_HEADER )

# Check for KSparse header
check_trilinos_header( "Epetra_CrsKundertSparse.h" HAS_KSPARSE_HEADER )

# Check for Trilinos Zoltan. 
# NOTE: Zoltan is required for Parallel Xyce builds!
set( SHOULD_ENABLE_ZOLTAN FALSE )
if( Trilinos_Zoltan_FOUND )
	set( SHOULD_ENABLE_ZOLTAN TRUE )
endif()
set( Xyce_USE_ZOLTAN ${SHOULD_ENABLE_ZOLTAN} CACHE BOOL "Use Zoltan library from Trilinos." )

if( Xyce_PARALLEL_MPI AND NOT Xyce_USE_ZOLTAN )
	message( SEND_ERROR "Parallel builds of Xyce require the Zoltan library in Trilinos.\n"
			"Try enabling this when building Trilinos: Trilinos_ENABLE_ZOLTAN=ON"
	)
endif()


if( Xyce_ENABLE_CXX11 )
	# We only check for the ShyLU header if C++11 is enabled.
	# Check for the ShyLU header
	check_trilinos_header( "Ifpack_ShyLU.h" HAS_SHYLU_HEADER )
endif()

# Note that Isorropia is required if Parallel MPI is enabled for Xyce!
if( Xyce_PARALLEL_MPI AND NOT Trilinos_Isorropia_FOUND )
	message( SEND_ERROR "Parallel Builds of Xyce require the Isorropia library \n"
				"in Trilinos! Try and enable it in Trilinos using:\n"
				"Trilinos_ENABLE_Isorropia=ON"
	)
endif()

# Xyce ML is deprecated and should likely be removed.
#set( Xyce_ML      ${HAS_ML_HEADER}      CACHE BOOL "Use ML Solver")

set( Xyce_KSPARSE FALSE CACHE BOOL "Include Kundert's Sparse direct solver." )
if( Xyce_KSPARSE AND NOT HAS_KSPARSE_HEADER )
	message( SEND_ERROR "KParse was explicitly enabled, but missing required header! To disable, use -DXyce_KSPARSE=OFF" )
endif()

# Set the ShyLU solver off by default. We can update this later.
#
# NOTE: I think this currently expects to use OpenMPI, because when I enabled it, I got
# errors involving redefinitions because "#include <mpi.h>" is (indirectly) included by
# the ShyLU stuff. Bear this in mind when enabling?
#
set( Xyce_SHYLU FALSE CACHE BOOL "Use ShyLU Solver" )
if( Xyce_SHYLU AND NOT HAS_SHYLU_HEADER )
	message( SEND_ERROR "ShyLU Solver is explicitly enabled, but the header does not exist! To disable, use -DXyce_SHYLU=NO" )
endif()

message( "Trilinos Details:" )
option_message( Xyce_AMESOS2 "Amesos2" )
option_message( Xyce_BELOS "Belos" )
option_message( Xyce_KSPARSE "KSparse" )
option_message( Xyce_SHYLU "ShyLU Solver" )
option_message( Xyce_USE_ISORROPIA "Isorropia" )
option_message( Xyce_USE_ZOLTAN "Zoltan" )
option_message( Xyce_NOX_USE_VECTOR_COPY "NOX Vector Copy" )
message( " " )

#-------------------------------------------------------
# Xyce Version Information
#-------------------------------------------------------
#
# NOTE: To minimize dependencies and the need to rerun CMake, we create custom CMake
# script that will call "configure_file" for the header file: "Xyce_version.h" in the
# current directory.
#
# This will prevent a rebuild of the entire project when updating the version or
# timestamp information, because CMake can identify header dependencies, and only the
# files that include the "Xyce_version.h" header will need to be rebuilt. If we set this
# via a preprocessor/compiler definition, every source file will need to be rebuilt,
# which seems like overkill of what we are doing.
#
set( Xyce_VERSION "D::7.0" )
set( Xyce_VERSION_MAJOR "7" )
set( Xyce_VERSION_MINOR "0" )
set( Xyce_VERSION_PATCH "0" )

add_definitions( "-DVERSION=\"${Xyce_VERSION}\"" )


#
# Create Xyce_version.h
#
# This file is created via a custom script, since we need to update the time of the
# build, as well as the Git hash (if applicable)
#
# The code below will create a file "GIT_HASH.txt" with the current git hash info
# from the last commit. If the hash changes, for whatever reason, then this will re-run
# cmake to make sure that the latest hash is actually available in GIT_HASH.txt.
#
set( GIT_HASH_FILE "${CMAKE_CURRENT_BINARY_DIR}/GIT_HASH.txt" )
if( Git_FOUND OR GIT_FOUND )
	include( GetGitRevisionDescription )
	get_git_short_commit_hash( GIT_HASH )
else()
	set( GIT_HASH "" )
endif()

# Set the current build time as well.
string( TIMESTAMP BUILD_TIME_STR "%Y%m%d%H%M" )
if( GIT_HASH )
	set( XYCEBUILDTIMESTAMP "${BUILD_TIME_STR}-g-${GIT_HASH}" )
else()
	set( XYCEBUILDTIMESTAMP "${BUILD_TIME_STR}" )
endif()

# Setup the Xyce_version.h file and install it.
set( XYCE_VERSION_FILEPATH "${CMAKE_CURRENT_BINARY_DIR}/src/Xyce_version.h" )
configure_file( "${CMAKE_CURRENT_SOURCE_DIR}/src/Xyce_version.h.cmake"
	"${XYCE_VERSION_FILEPATH}"
)
install( FILES ${XYCE_VERSION_FILEPATH}
	DESTINATION include
)

#-------------------------------------------------------
# Xyce Debug Options
#-------------------------------------------------------
#
option( Xyce_DEBUG "Enable all Xyce debugging." OFF )

# Define a macro to make defining debug and verbose options easier.
macro(xyce_debug_option argVar argName helpStr defaultVal)
	if( Xyce_DEBUG )
		set( ${defaultVal} ON )
	endif()
	option( ${argVar} "${helpStr}" ${defaultVal} )
	message("-- Enable ${argName} debugging (${argVar}): ${${argVar}}" )
endmacro()

# Since the macros will print the option status as well, let's start printing the
# header here.
message( "Xyce Debug Options:" )

# Note that option is shorthand for:
# set( <var> <default val> CACHE BOOL <doc string> )
xyce_debug_option( Xyce_DEBUG_DEVICE      "Device"      "Enable Device Debug."    OFF )
xyce_debug_option( Xyce_DEBUG_ANALYSIS    "Analysis"    "Enable Analysis Debug."  OFF )
xyce_debug_option( Xyce_DEBUG_HB          "HB"          "Enable HB Debug."        OFF )
xyce_debug_option( Xyce_DEBUG_MPDE        "MPDE"        "Enable MPDE Debug."      OFF )
xyce_debug_option( Xyce_DEBUG_MOR         "MOR"         "Enable MOR Debug."       OFF )
xyce_debug_option( Xyce_DEBUG_IO          "IO"          "Enable Xyce IO Debug."   OFF )
xyce_debug_option( Xyce_DEBUG_EXPRESSION  "Expression"  "Enable Expression Debug." OFF )
xyce_debug_option( Xyce_DEBUG_CONDUCTANCE "Conductance" "Enable Conductance Debug." OFF )
xyce_debug_option( Xyce_DEBUG_RESTART     "Restart"     "Enable Debug Restart." OFF )
xyce_debug_option( Xyce_DEBUG_TIME        "Time"        "Enable Time Debug." OFF )
xyce_debug_option( Xyce_DEBUG_CIRCUIT     "Circuit"     "Enable Circuit Debug."   OFF )
xyce_debug_option( Xyce_DEBUG_NONLINEAR   "Nonlinear"   "Enable Nonlinear Debug." OFF )
xyce_debug_option( Xyce_DEBUG_LINEAR      "Linear"      "Enable Linear Debug."    OFF )
xyce_debug_option( Xyce_DEBUG_VOLTLIM     "Voltlim"     "Enable Voltlim Debug."   OFF )
xyce_debug_option( Xyce_DEBUG_TOPOLOGY    "Topology"    "Enable Topology Debug."  OFF )
xyce_debug_option( Xyce_DEBUG_PARALLEL    "Parallel"    "Enable Parallel processing debug." OFF )
xyce_debug_option( Xyce_DEBUG_DISTRIBUTION "Distribution" "Enable Distribution debug." OFF )
#option( Xyce_DEBUG_ALL_PROCS_SAME_WD "Enable something?" OFF )
#option( Xyce_Dakota_Debug "Enable Dakota Debug." OFF )

# For formatting.
message( " " )

macro( xyce_verbose_option argVar argName helpStr defaultVal)
	if( Xyce_DEBUG )
		set( ${defaultVal} ON )
	endif()
	option( ${argVar} "${helpStr}" ${defaultVal} )
	message("-- Enable ${argName} verbose output (${argVar}): ${${argVar}}" )
endmacro()

message( "Xyce Verbosity Options:")
xyce_verbose_option( Xyce_VERBOSE_CONDUCTANCE "Conductance" "Enable verbose conductance output." OFF )
xyce_verbose_option( Xyce_VERBOSE_TIME "Time" "Enable verbose time output." OFF )
xyce_verbose_option( Xyce_VERBOSE_NONLINEAR "Nonlinear" "Enable verbose nonlinear output." OFF )
xyce_verbose_option( Xyce_VERBOSE_LINEAR "Linear" "Enable verbose linear output." OFF )

# For formatting.
message( " " )

#-------------------------------------------------------
# Platform-Specific Details
#-------------------------------------------------------
#
# Setting up the RPATH stuff.
#
set( CMAKE_MACOSX_RPATH 1 )

# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 
set( CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib" )

# Don't add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set( CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE )

list( FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if( "${isSystemDir}" STREQUAL "-1" )
	set( CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib" )
endif()

#-------------------------------------------------------
# Xyce Libraries
#-------------------------------------------------------
#
# The code for the Xyce libraries are listed out in the subdirectories below.
# These subdirectories should use the "add_to_xyce_library" macro to add their files 
# into the final library, so that these objects are registered in the
# "XYCE_LIBRARY_IMPL_LIST" CMake variable that is used to actually build Xyce.
#
add_subdirectory( src )

#-------------------------------------------------------
# Xyce Installer Details.
#-------------------------------------------------------
#
# First, we should calculate the license that we are using.
#
set( Xyce_LICENSE ECI )
set( Xyce_RPM_LICENSE "Export Controlled" )
if( NOT Xyce_RAD_MODELS )
	if( NOT Xyce_NONFREE_MODELS )
		# If we aren't adding radiation or "nonfree" models, the license is opensource.
		set( Xyce_LICENSE OS )
		set( Xyce_RPM_LICENSE "GPLv3" )
	else()
		set( Xyce_LICENSE NonFree )
		set( Xyce_RPM_LICENSE "NDA Required" )
	endif()
endif()

# Set the path to the license file, and also install it for reference purposes.
set( Xyce_LICENSE_FILE "${CMAKE_CURRENT_SOURCE_DIR}/distribution/CPack.${Xyce_LICENSE}License.txt" )
install( FILES ${Xyce_LICENSE_FILE}
	DESTINATION share

	# Rename the file to "LICENSE" so it is more easily recognizable as to what it is.
	RENAME "LICENSE.txt"
)

# Print out our packaging details.
message( " " )
message( "Packaging Details:" )
message( "-- License type: ${Xyce_LICENSE} (${Xyce_RPM_LICENSE})" )
message( "-- * Path: ${Xyce_LICENSE_FILE}" )
message( "-- Install Prefix: ${CMAKE_INSTALL_PREFIX}" )

# Include any required system libraries in our installation.
include (InstallRequiredSystemLibraries)

#------------------------------
# Determine the CPack Generators to use.
#
set( CPACK_GENERATOR )

# Define a macro to make adding installer options more systematic.
# Arguments are:
# - _name_var -- The human-readable name of the installer.
# - _install_var -- The variable to store whether its enabled.
# - _cpack_gen_var -- The name of the installer generator for CPack (i.e. ZIP, TGZ, etc.)
# - _desc_var -- The help/description text for the installer.
# - _def_var -- The default value to assign.
#
macro( add_installer_option _name_var _install_var _cpack_gen_var _desc_var _def_var )
	option( ${_install_var} "${_desc_var}" ${_def_var} )
	message( "-- Enable ${_name_var} installer (${_install_var}): ${${_install_var}}" )

	# Add the option to the generator list, if it is enabled.
	if( ${_install_var} )
		list( APPEND CPACK_GENERATOR ${_cpack_gen_var} )
	endif()
endmacro()

# Zip and TGZ should be enabled on all systems.
add_installer_option( Zip Xyce_INSTALL_ZIP ZIP "Enable to create a Zip installer." ON )
add_installer_option( TGZ Xyce_INSTALL_TGZ TGZ "Enable to create a TGZ installer." ON )

# Self-contained script installers should be enabled on UNIX OS's.
if( UNIX )
	add_installer_option( "Script Installer" Xyce_INSTALL_STGZ STGZ
		"Enable for a self-contained script installer." ON )

	# Check for the RPM installer.
	add_installer_option( "RPM Installer" Xyce_INSTALL_RPM RPM
		"Enable for an RPM installer." OFF )
endif()


# We should try to add the NSIS installer on Windows.
if( WIN32 )
	add_installer_option( "Script Installer" Xyce_INSTALL_NSIS NSIS
		"Enable for a self-contained script installer." OFF )
endif()


# Set the CPack Package information for Xyce.
set( CPACK_PACKAGE_NAME "Xyce" )
set( CPACK_PACKAGE_VENDOR "Sandia National Laboratories" )
set( CPACK_PACKAGE_VERSION_MAJOR "${Xyce_VERSION_MAJOR}" )
set( CPACK_PACKAGE_VERSION_MINOR "${Xyce_VERSION_MINOR}" )
set( CPACK_PACKAGE_VERSION_PATCH "${Xyce_VERSION_PATCH}" )
set( CPACK_PACKAGE_INSTALL_DIRECTORY "Xyce" )
set( CPACK_PACKAGE_VERSION "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}" )

# Set some common resource files
set( CPACK_RESOURCE_FILE_LICENSE "${Xyce_LICENSE_FILE}" )
set( CPACK_RESOURCE_FILE_WELCOME "${CMAKE_CURRENT_SOURCE_DIR}/distribution/CPack.Welcome.txt" )
set( CPACK_RESOURCE_FILE_README  "${CMAKE_CURRENT_SOURCE_DIR}/distribution/CPack.Description.txt" )

# Set CPack RPM Information
set( CPACK_RPM_PACKAGE_LICENSE "${Xyce_RPM_LICENSE}" )
set( CPACK_RPM_PACKAGE_DESCRIPTION " For more information, visit http://xyce.sandia.gov" )
set( CPACK_RPM_PACKAGE_GROUP "Applications/Engineering" )
set( CPACK_RPM_PACKAGE_NAME "${CPACK_PACKAGE_NAME}" )
set( CPACK_RPM_EXCLUDE_FROM_AUTO_FILELIST_ADDITION "/usr/local" )

# Set the CPackage resources.
set( CPACK_PACKAGE_RELOCATABLE TRUE )

# Set parameters specific to NSIS
set( CPACK_NSIS_DISPLAY_NAME "Xyce" )
set( CPACK_NSIS_PACKAGE_NAME "Xyce" )

message( " " )

if( PRINT_DEBUG )
	debug_print()
	# Unset the PRINT_DEBUG variable, since we don't want the variables to keep printing
	# all the time.
	unset( PRINT_DEBUG CACHE )
endif()

# NOTE: This MUST be the last line in our CMake build!
include( CPack )

