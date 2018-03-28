#
# CMakeLists.txt for Xyce
#
# Build Xyce!
#

# The minimum required version is currently 2.8.8. We need to support RHEL6.
# If you have the chance, it is better to use a newer version of CMake.
cmake_minimum_required( VERSION 2.8.8 )

project( Xyce )

# Add our custom CMake scripts to the module search path.
list( APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" )

# Store a variable with our linking libraries for Xyce.
list( APPEND Xyce_LINK_LIBS )

#-------------------------------------------------------
# CMake Options and Includes
#-------------------------------------------------------
#
# Define some macros here that are used throughout the Xyce build.
# Note that we use CMakeParseArguments, available in CMake 2.8.3 and newer.
include( CMakeParseArguments )
include( CheckIncludeFileCXX )
include( CheckCXXSymbolExists )
include( CheckCXXSourceCompiles )

# Set the main directory for the source of our custom platform check programs.
set( PLATFORM_CHECK_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/cmake/platformCheck" )

#-------------------------------------------------------
# Find Trilinos
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
			# user may have explicitly set prior.
			set( Trilinos_DIR "${Trilinos_DIR}/lib/cmake/Trilinos" CACHE PATH "Path to the directory containing TrilinosConfig.cmake." FORCE )
		endif()
	endif()
endif()

find_package( Trilinos CONFIG REQUIRED )
include_directories( ${Trilinos_INCLUDE_DIRS} )

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
set( Xyce_LINK_EXE_WITH_SHARED ${Xyce_BUILD_SHARED} CACHE BOOL "Link Xyce executable with shared library." )

# Check if Trilinos is compiled to use OpenMPI. If so, we try and "reuse" those
# settings here for consistency.
set( SHOULD_USE_MPI FALSE )
if( EXISTS "${Trilinos_MPI_EXEC}" )
	set( SHOULD_USE_MPI TRUE )
endif()
set( Xyce_PARALLEL_MPI ${SHOULD_USE_MPI} CACHE BOOL "Build Xyce to use MPI." )

if( Xyce_PARALLEL_MPI )
	# Note that the caller may have explicitly set the MPI_<lang>_COMPILER
	# flags as a hint already. We only default to use what Trilinos is using
	# if these parameters weren't explicitly set.
	if( NOT MPI_C_COMPILER )
		set( MPI_C_COMPILER "${Trilinos_C_COMPILER}" )
	endif()
	if( NOT MPI_CXX_COMPILER )
		set( MPI_CXX_COMPILER "${Trilinos_CXX_COMPILER}" )
	endif()
	#if( NOT MPI_Fortran_COMPILER )
	#	set( MPI_Fortran_COMPILER "${Trilinos_Fortran_COMPILER}" )
	#endif()
	
	# Try to find MPI with our current settings. 
	find_package( MPI )
		
	# If we found our libraries, set that we should build with MPI.
	if( NOT MPI_C_FOUND )
		message( FATAL_ERROR "Xyce is expected to build with MPI, but no MPI C compiler found! Please specify this directly via: MPI_C_COMPILER." )
	endif()
	if( NOT MPI_CXX_FOUND )
		message( FATAL_ERROR "Xyce is expected to build with MPI, but no MPI C++ compiler found! Please specify this directly via: MPI_CXX_COMPILER." )
	endif()

	set( CMAKE_C_COMPILER "${MPI_C_COMPILER}" )
	set( CMAKE_CXX_COMPILER "${MPI_CXX_COMPILER}" )
	include_directories( SYSTEM ${MPI_CXX_INCLUDE_DIRS} )
	list( APPEND Xyce_LINK_LIBS ${MPI_CXX_LIBRARIES} )
endif()

# Check if Trilinos is compiled with C++11 enabled. We need to build Xyce with
# the same standard for C++, so check that before doing anything else.
set( SHOULD_USE_CXX11 TRUE )
if( Trilinos_ENABLE_CXX11 )
	message( "-- Enabling CXX11 for Xyce, since it is enabled in Trilinos." )
	set( SHOULD_USE_CXX11 TRUE )
endif()
set( Xyce_ENABLE_CXX11 ${SHOULD_USE_CXX11} CACHE BOOL "Enable C++11 for Xyce" )

if( Xyce_ENABLE_CXX11 )
	set( CMAKE_CXX_STANDARD 11 )
endif()

#-------------------------------------------------------
# Find Other Packages (Non-Trilinos)
#-------------------------------------------------------
#

find_package( CURL )

# Handle the tracking URL here.
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

# Check for any external FFT libraries.
# Check Intel FFT.
#if ( Xyce_USE_INTEL_FFT )
#	set( FFT_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/N_UTL_IntelFFT_Interface.C )
#endif()
#if ( Xyce_USE_FFTW OR NOT Xyce_USE_INTEL_FFT )
#	set( FFTW_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/N_UTL_FFTW_Interface.C )
#endif()

# Look for BLAS

# Look for LAPACK

# Look for MPI

# Check for Intel headers.
check_include_file_cxx( "mkl_dfti.h" Xyce_USE_INTEL_FFT )
if( Xyce_USE_INTEL_FFT )
	message("-- Will enable Intel FFT")
endif()

find_package( FFTW )
if( FFTW-NOTFOUND )
	message( "-- Could not find FFTW libraries." )
	set( FOUND_FFTW FALSE )
else()
	include_directories( ${FFTW_INCLUDE_DIRS} )

	set( FOUND_FFTW TRUE )
endif()
set( Xyce_USE_FFTW ${FOUND_FFTW} CACHE BOOL "Enable the FFTW library for FFT." )

find_package( FLEX )
find_package( BISON )

find_package( Dakota CONFIG QUIET)
if( Dakota-NOTFOUND )
	message( "-- Did NOT find Dakota installation." )
	# Don't use Dakota.
endif()

#-------------------------------------------------------
# Print Initial Config and Third-Party Library Results
#-------------------------------------------------------
message( "Build Details:" )
message( "-- Enable C++11: ${Xyce_ENABLE_CXX11}" )
message( "-- Build Shared libraries: ${Xyce_BUILD_SHARED}" )
message( "-- Build Static libraries: ${Xyce_BUILD_STATIC}" )
message( "-- Link Xyce against Shared: ${Xyce_LINK_EXE_WITH_SHARED}" )
message( "-- Enable MPI: ${Xyce_PARALLEL_MPI}" )
message( " " )
message( "Compiler Details:")
message( "-- C Compiler  : ${CMAKE_C_COMPILER}" )
message( "-- CXX Compiler: ${CMAKE_CXX_COMPILER}" )
message( " " )
message( "Third Party Libraries:" )
message( "-- Trilinos Include Dir: ${Trilinos_INCLUDE_DIRS}" )
if( FLEX_FOUND )
	message( "-- Found Flex (${FLEX_VERSION}): ${FLEX_EXECUTABLE}" )
endif()
if( BISON_FOUND )
	message( "-- Found Bison (${BISON_VERSION}): ${BISON_EXECUTABLE}" )
endif()
message( "-- Found FFTW: ${Xyce_USE_FFTW}" )
if( Xyce_USE_FFTW )
	message( "-- * FFTW Include Path: ${FFTW_INCLUDE_DIRS}" )
	message( "-- * FFTW Libraries: ${FFTW_LIBRARIES}" )
endif()
message( "-- Enable Usage Tracking: ${Xyce_USE_CURL}" )
if( Xyce_USE_CURL )
	message( "-- * Curl Version: ${CURL_VERSION_STRING}" )
	message( "-- * Curl Include Path: ${CURL_INCLUDE_DIRS}" )
	message( "-- * Curl Library Path: ${CURL_LIBRARIES}" )
	message( "-- * Tracking URL: ${Xyce_TRACKING_URL}" )
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

# This macro invokes check_cxx_source_compiles on the given file.
# 
# Unfortunately, check_cxx_source_compiles() expects the actual source code that
# is to be compiled to be stored in the variable, not just the filepath. This 
# quietly reads the file automatically before invoking check_cxx_source_compiles()
macro( check_cxx_source_compiles_file COMPILATION_FILE CMAKE_VAR_TO_STORE )
	file( READ "${COMPILATION_FILE}" TMPSRC )
	check_cxx_source_compiles( "${TMPSRC}" ${CMAKE_VAR_TO_STORE} )
endmacro()

# Check for <unordered_map>
check_include_file_cxx( "unordered_map" HAVE_UNORDERED_MAP )
check_include_file_cxx( "tr1/unordered_map" HAVE_TR1_UNORDERED_MAP )

# Check for <unordered_set>
check_include_file_cxx( "unordered_set" HAVE_UNORDERED_SET )
check_include_file_cxx( "tr1/unordered_set" HAVE_TR1_UNORDERED_SET )

# Check for C++ common headers.
check_include_file_cxx( "numeric" HAVE_NUMERIC )
check_include_file_cxx( "algorithm" HAVE_ALGORITHM )
check_include_file_cxx( "functional" HAVE_FUNCTIONAL )
check_include_file_cxx( "cstdlib" HAVE_CSTDLIB )
check_include_file_cxx( "iostream" HAVE_IOSTREAM )
check_include_file_cxx( "fstream" HAVE_FSTREAM )

check_include_file_cxx( "unistd.h" HAVE_UNISTD_H )
check_include_file_cxx( "dlfcn.h" HAVE_DLFCN_H )
check_include_file_cxx( "float.h" HAVE_FLOAT_H )
check_include_file_cxx( "strings.h" HAVE_STRINGS_H )
check_include_file_cxx( "string.h" HAVE_STRING_H )
check_include_file_cxx( "inttypes.h" HAVE_INTTYPES_H ) 
check_include_file_cxx( "limits.h" HAVE_LIMITS_H )
check_include_file_cxx( "malloc.h" HAVE_MALLOC_H )
check_include_file_cxx( "math.h" HAVE_MATH_H )
check_include_file_cxx( "mathimf.h" HAVE_MATHIMF_H )

# Check for Windoze header.
check_include_file_cxx( "Windows.h" HAVE_WINDOWS_H )

# Check for sys/* headers.
check_include_file_cxx( "sys/resource.h" HAVE_SYS_RESOURCE_H )
check_include_file_cxx( "sys/stat.h" HAVE_SYS_STAT_H )
check_include_file_cxx( "sys/types.h" HAVE_SYS_TYPES_H )

# Check for std::iota function.
check_cxx_source_compiles_file( "${PLATFORM_CHECK_SRC_DIR}/CheckForIota.cxx" HAVE_IOTA )

# Check for linux exceptions
check_cxx_source_compiles_file( "${PLATFORM_CHECK_SRC_DIR}/CheckForLinuxExceptions.cxx" HAVE_LINUX_EXCEPTIONS )

#-------------------------------------------------------
# Trilinos Checks
#-------------------------------------------------------
#
# Define a special macro to compile source with Trilinos on the include path.
macro( check_cxx_source_compiles_with_trilinos COMPILATION_FILE CMAKE_VAR_TO_STORE )
	set( _BACKUP_INC ${CMAKE_REQUIRED_INCLUDES} )
	set( _BACKUP_LIB ${CMAKE_REQUIRED_LIBRARIES} )

	set( CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES} ${Trilinos_INCLUDE_DIRS} )
	set( CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${Trilinos_LIBRARIES} )

	file( READ "${COMPILATION_FILE}" TMPSRC )
	check_cxx_source_compiles( "${TMPSRC}" ${CMAKE_VAR_TO_STORE} )

	set( CMAKE_REQUIRED_INCLUDES ${_BACKUP_INC} )
	set( CMAKE_REQUIRED_LIBRARIES ${_BACKUP_LIB} )
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
	message( SEND_ERROR "Xyce requires the Teuchos package of Trilinos.\n"
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
if ( Trilinos_Belos_FOUND )
	set( Xyce_BELOS ON CACHE BOOL "Use Belos" )
endif()

# NOX FOUND -- If NOX stuff was found, enable features.
if( Trilinos_NOX_FOUND )
	set( Xyce_NOX_USE_VECTOR_COPY ON )
endif()

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
message( "-- Enable Amesos2: ${Xyce_AMESOS2}" )
message( "-- Enable Belos: ${Xyce_BELOS}" )
message( "-- Enable KSparse: ${Xyce_KSPARSE}" )
message( "-- Enable ShyLU Solver: ${Xyce_SHYLU}" )
message( "-- Enable Isorropia: ${Xyce_USE_ISORROPIA}" )
message( "-- Enable Zoltan: ${Xyce_USE_ZOLTAN}" )
message( "-- Using NOX Vector Copy: ${Xyce_NOX_USE_VECTOR_COPY}" )
message( " " )

#-------------------------------------------------------
# Xyce Version Information
#-------------------------------------------------------
set( Xyce_VERSION "D::7.0" )
set( Xyce_VERSION_MAJOR "7" )
set( Xyce_VERSION_MINOR "0" )
set( Xyce_VERSION_PATCH "0" )
add_definitions( "-DVERSION=\"${Xyce_VERSION}\"" )

add_definitions( -DXYCEBUILDTIMESTAMP )

#-------------------------------------------------------
# Platform-Specific Details
#-------------------------------------------------------
if( APPLE )
	set( CMAKE_MACOSX_RPATH 1 )
endif()

# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 
set( CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib" )


# Don't add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set( CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE )

#list( FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
#if( "${isSystemDir}" STREQUAL "-1" )
#	set( CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib" )
#endif()


# Now, create a configuration file that we will use.
set( XYCE_CONFIG_FILEPATH "${CMAKE_CURRENT_BINARY_DIR}/Xyce_config.h" )
configure_file( "${CMAKE_CURRENT_SOURCE_DIR}/src/Xyce_config.h.cmake" "${XYCE_CONFIG_FILEPATH}" )
install( FILES "${XYCE_CONFIG_FILEPATH}"
	DESTINATION include
)
include_directories( "${CMAKE_CURRENT_BINARY_DIR}" )

# Useful for debugging!
if( PRINT_DEBUG )
get_cmake_property(_variableNames VARIABLES)
list (SORT _variableNames)
foreach (_variableName ${_variableNames})
    message(STATUS "${_variableName}=${${_variableName}}")
endforeach()
endif()


#Xyce_NOX_USE_VECTOR_COPY

#-------------------------------------------------------
# Macro Definitions
#-------------------------------------------------------

#
# add_to_xyce_library macro
#
# Usage:
# add_to_xyce_library( NAME XycePackageA
#     FILES
#         packageAFile.cpp
#         packageBFile.cpp
# )
#
# This macro is designed to add source files in some package into the final Xyce library
# that is compiled. This macro will automatically handle both the shared and static
# library linking case.
#
macro( add_to_xyce_library )
	set( _OPTIONS_ARGS )
	set( _ONE_VALUE_ARGS )
	set( _MULTI_VALUE_ARGS NAME FILES )

	cmake_parse_arguments( _ADDXYCELIB "${_OPTIONS_ARGS}" "${_ONE_VALUE_ARGS}" "${_MULTI_VALUE_ARGS}" "${ARGN}")

	if( NOT _ADDXYCELIB_NAME )
		message( FATAL_ERROR "add_to_xyce_library: 'NAME' argument required.")
	endif()
	
	list( GET _ADDXYCELIB_NAME 0 LIBNAME )

	# WHEW, CMake object libraries are supported in CMake 2.8.8, which is just barely
	# our required version. The Xyce build system relies on this feature heavily.
	#
	# I spoke too soon -- Note that we could use the same object library for both the
	# static and shared libraries IF AND ONLY IF we permit the static library to be
	# compiled with -fPIC or similar (allow position-independent code). This flag is
	# necessary for shared libraries, and only minorly affects some optimization for
	# static libraries, so we'll add it here.
	#
	# The CMake way to do this is via the "POSITION_INDEPENDENT_CODE" property, but this
	# property only exists in CMake 2.8.9 onwards, so we need to handle a stupid case
	# just to be safe. If the CMake version is incremented, we can get rid of the check.
	#
	add_library( "${LIBNAME}" OBJECT "${_ADDXYCELIB_FILES}")
	if( ${CMAKE_VERSION} EQUAL "2.8.8" )
		# The stupid case here.
		set_property( TARGET "${LIBNAME}" PROPERTY COMPILE_FLAGS "-fPIC" )
	else()
		# The usual/desired case here. Using the CMake property directly is safer,
		# since it should work across multiple platforms with no additional checks.
		set_property( TARGET "${LIBNAME}" PROPERTY POSITION_INDEPENDENT_CODE ON )
	endif()
	
	# Add the created object library into our global list of libraries, so that we
	# include it when we create the final Xyce libraries.
	set_property( GLOBAL APPEND PROPERTY xyce_object_library_list "${LIBNAME}" )
endmacro()

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
set( CPACK_PACKAGE_NAME "Xyce" )
set( CPACK_PACKAGE_VENDOR "Sandia National Laboratories" )
set( CPACK_PACKAGE_VERSION_MAJOR "${Xyce_VERSION_MAJOR}" )
set( CPACK_PACKAGE_VERSION_MINOR "${Xyce_VERSION_MINOR}" )
set( CPACK_PACKAGE_VERSION_PATCH "${Xyce_VERSION_PATCH}" )
set( CPACK_PACKAGE_INSTALL_DIRECTORY "Xyce" )
set( CPACK_PACKAGE_VERSION "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}" )

# Unset the PRINT_DEBUG variable, since we don't want the variables to keep printing
# all the time.
unset( PRINT_DEBUG CACHE )

# NOTE: This MUST be the last line in our CMake build!
include( CPack )

