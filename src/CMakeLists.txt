#
# CMakeLists.txt for src
#
# Define all of the packages to include in the current build of Xyce
#

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

# These packages are always required.
#
# NOTE: The order of these packages matters greatly!
#
# Only change the order if you know what you are doing...
#
add_xyce_subdirectory( AnalysisPKG )
add_xyce_subdirectory( CircuitPKG )
add_xyce_subdirectory( DakotaLinkPKG )

add_xyce_subdirectory( DeviceModelPKG/Core )
add_xyce_subdirectory( DeviceModelPKG/OpenModels )
add_xyce_subdirectory( DeviceModelPKG/TCADModels )
add_xyce_subdirectory( DeviceModelPKG/EXTSC )

add_xyce_subdirectory( ErrorHandlingPKG )

#add_xyce_subdirectory( IOInterfacePKG/Output )
include_directories( IOInterfacePKG/Output)
add_xyce_subdirectory( IOInterfacePKG )

add_xyce_subdirectory( LinearAlgebraServicesPKG )
add_xyce_subdirectory( LoaderServicesPKG )
add_xyce_subdirectory( MultiTimePDEPKG )
add_xyce_subdirectory( NonlinearSolverPKG )
add_xyce_subdirectory( ParallelDistPKG )
add_xyce_subdirectory( TimeIntegrationPKG )
add_xyce_subdirectory( TopoManagerPKG )
add_xyce_subdirectory( UtilityPKG )

# Iterate through all of our directories. We need to iterate to setup the include paths first. 
# Then, we iterate over everything again to actually build the object libraries.
get_property( xyce_subpackage_dirs GLOBAL PROPERTY xyce_subpackages )
foreach( DIR ${xyce_subpackage_dirs} )
	include_directories( "${DIR}" )
endforeach()
foreach( DIR ${xyce_subpackage_dirs} )
	add_subdirectory( "${DIR}" )
endforeach()


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

# Create the main executable.
add_executable( XyceExe
	Xyce.C
)
# If configured to do so, link against the shared library. Otherwise,
# default to linking against the static library.
if( Xyce_LINK_EXE_WITH_SHARED )
	list( APPEND Xyce_LINK_LIBS XyceShared )
else()
	list( APPEND Xyce_LINK_LIBS XyceStatic )
endif()

if( HAVE_DLFCN_H )
	list( APPEND Xyce_LINK_LIBS ${CMAKE_DL_LIBS} )
endif()

target_link_libraries( XyceExe
	${Xyce_LINK_LIBS}
)
set_target_properties( XyceExe PROPERTIES OUTPUT_NAME Xyce )

install( TARGETS XyceExe
	DESTINATION bin
)

