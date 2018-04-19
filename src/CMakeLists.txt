#
# CMakeLists.txt for src
#
# Define all of the packages to include in the current build of Xyce
#
# Author: Aaron Gibson (asgibso@sandia.gov)
#

# Include the current binary directory, because we may generate some configuration
# files that are needed everywhere else.
include_directories( "${CMAKE_CURRENT_BINARY_DIR}" )

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
	set( _ONE_VALUE_ARGS NAME )
	set( _MULTI_VALUE_ARGS FILES )

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

#
# add_xyce_subdirectory
#
# This calls "add_subdirectory", but adds that directory as an include path
# before calling it.
#
# This is necessary due to the way that the header files are handled in Xyce.
#
macro( add_xyce_subdirectory DIR_TO_ADD )
	set_property( GLOBAL APPEND PROPERTY xyce_subpackages "${DIR_TO_ADD}" )
endmacro()

#-------------------------------------------------------
# Parse Xyce Object Libraries
#-------------------------------------------------------
#
add_xyce_subdirectory( AnalysisPKG )
add_xyce_subdirectory( CircuitPKG )
add_xyce_subdirectory( DakotaLinkPKG )
add_xyce_subdirectory( ErrorHandlingPKG )
add_xyce_subdirectory( LinearAlgebraServicesPKG )
add_xyce_subdirectory( LoaderServicesPKG )
add_xyce_subdirectory( MultiTimePDEPKG )
add_xyce_subdirectory( NonlinearSolverPKG )
add_xyce_subdirectory( ParallelDistPKG )
add_xyce_subdirectory( TimeIntegrationPKG )
add_xyce_subdirectory( TopoManagerPKG )
add_xyce_subdirectory( UtilityPKG )

# Unfortunately, we need to treat this package differently.
# I'm not entirely sure why we have an "Output" directory that also has include files
# that we need elsewhere, but we do :/
#
# We should just consider moving this into a different PKG directory altogether for
# simplicity (maybe IOInterfaceOutputPKG?)
# As a kludgy workaround, we'll include the nested 'Output' directory here directly.
#
include_directories( IOInterfacePKG/Output )
add_xyce_subdirectory( IOInterfacePKG )

# Handle some other projects that may recursively call add_xyce_subdirectory.
# We don't want to add "DeviceModelPKG" to the include path directly because there
# aren't any headers to include there.
add_subdirectory( DeviceModelPKG )

#-------------------------------------------------------
# Build Xyce from Object Libraries
#-------------------------------------------------------
#
# Iterate through all of our directories. We need to iterate to setup the include paths first. 
# Then, we iterate over everything again to actually build the object libraries.
get_property( xyce_subpackage_dirs GLOBAL PROPERTY xyce_subpackages )
foreach( DIR ${xyce_subpackage_dirs} )
	include_directories( "${DIR}" )
endforeach()
foreach( DIR ${xyce_subpackage_dirs} )
	# We need to explicitly determine a binary directory, since some of these
	# directories may be included from outside the source fo the original
	# build tree.
	get_filename_component( XYCE_BINARY_PKG_NAME "${DIR}" NAME_WE )
	add_subdirectory( "${DIR}" "${CMAKE_CURRENT_BINARY_DIR}/${XYCE_BINARY_PKG_NAME}" )
endforeach()

#-------------------------------------------------------
# Create Xyce configuration files
#-------------------------------------------------------
#
# Create Xyce_config.h
#
# NOTE: We CANNOT do this any earlier, because some of our "add_subdirectory" commands
# above may add configuration options themselves that should be visible in this file!
# The file will be generated at configure time, so we shouldn't need to worry about
# any problems that arise.
set( XYCE_CONFIG_FILEPATH "${CMAKE_CURRENT_BINARY_DIR}/Xyce_config.h" )
configure_file( "${CMAKE_CURRENT_SOURCE_DIR}/Xyce_config.h.cmake" "${XYCE_CONFIG_FILEPATH}" )
install( FILES "${XYCE_CONFIG_FILEPATH}"
        DESTINATION include
)

#
# Build the final Xyce library artifacts
#
# We do this by parsing through "xyce_object_library_list" to determine all of the 
# objects that we should include in the final Xyce library artifacts. The proper syntax
# to use "Object Libraries" in CMake requires using generator expressions, which can 
# appear a bit cryptic.
#
# The for-loop below will store the resulting generator expressions in the variable:
# "FINAL_XYCE_LIBRARY_OBS", so that we can use CMake's 'add_library' command as expected.
#
# More details about object libraries in CMake can be found here:
# https://cmake.org/Wiki/CMake/Tutorials/Object_Library
#
get_property(xyce_obj_list GLOBAL PROPERTY xyce_object_library_list )
foreach( LIB ${xyce_obj_list} )
	set( FINAL_XYCE_LIBRARY_OBS "${FINAL_XYCE_LIBRARY_OBS};\$<TARGET_OBJECTS:${LIB}>" )
endforeach()

# Build the Shared libraries, if applicable.
if( Xyce_BUILD_SHARED )
	add_library( XyceShared SHARED ${FINAL_XYCE_LIBRARY_OBS} )
	target_link_libraries( XyceShared
		${Trilinos_LIBRARIES}
		${FFTW_LIBRARIES}
		${Xyce_LINK_LIBS}
	)
	set_target_properties( XyceShared PROPERTIES OUTPUT_NAME Xyce )
	install( TARGETS XyceShared
		DESTINATION lib
	)
endif()

# Build the Static libraries, if applicable.
if( Xyce_BUILD_STATIC )
	add_library( XyceStatic STATIC ${FINAL_XYCE_LIBRARY_OBS} )
	target_link_libraries( XyceStatic
		${Trilinos_LIBRARIES}
		${FFTW_LIBRARIES}
		${Xyce_LINK_LIBS}
	)
	
	if( WIN32 )
		# Change the static library to use the "libXyce.lib" convention on Windows.
		# This is necessary, because Windows uses the ".lib" file extension to denote
		# DLL import libraries and static libraries both. Windoze SUX!
		set_target_properties( XyceStatic PROPERTIES OUTPUT_NAME libXyce )
	else()
		set_target_properties( XyceStatic PROPERTIES OUTPUT_NAME Xyce )
	endif()
	
	install( TARGETS XyceStatic
		DESTINATION lib
	)
endif()

# Create the main executable, and make sure it is named "Xyce" as we'd expect.
add_executable( XyceExe
	Xyce.C
)
set_target_properties( XyceExe PROPERTIES OUTPUT_NAME Xyce )

# If configured to do so, link against the shared library. Otherwise,
# default to linking against the static library.
if( Xyce_LINK_EXE_WITH_SHARED )
	list( APPEND Xyce_LINK_LIBS XyceShared )
else()
	list( APPEND Xyce_LINK_LIBS XyceStatic )
endif()

# Link the executable.
target_link_libraries( XyceExe
	${Xyce_LINK_LIBS}
)

# Install it.
install( TARGETS XyceExe
	BUNDLE DESTINATION .
	RUNTIME DESTINATION bin
)

