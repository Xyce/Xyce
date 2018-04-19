#
# FindMKL.cmake
#
# Author: Aaron Gibson (asgibso@sandia.gov)
#
# I borrowed much of this file from:
# https://gist.github.com/scivision/5108cf6ab1515f581a84cd9ad1ef72aa
#
# but modified it to be more CMake-friendly, by using imported targets and so 
# on instead of directly touching the c++ flags and link lines.
#
#
# Find the Math Kernel Library from Intel
#
#  MKL_FOUND - System has MKL
#  MKL_INCLUDE_DIRS - MKL include files directories
#  MKL_LIBRARIES - The MKL libraries
#  MKL_INTERFACE_LIBRARY - MKL interface library
#  MKL_CORE_LIBRARY - MKL core library
#  MKL_OPENMP_LAYER_LIBRARY - MKL OpenMP Threading layer library.
#
#  Don't use?
#  MKL_SEQUENTIAL_LAYER_LIBRARY - MKL sequential layer library
#
#  The environment variables MKLROOT and INTEL are used to find the library.
#
#  Example usage:
#
#  find_package(MKL)
#  if(MKL_FOUND)
#    target_link_libraries(TARGET ${MKL_LIBRARIES})
#  endif()

# If already in cache, be silent
if( MKL_INCLUDE_DIRS
	AND MKL_LIBRARIES
)
	set (MKL_FIND_QUIETLY TRUE)
endif()

# MKLROOT is an important variable, so make sure it is set.
if( NOT DEFINED ${MKLROOT} )
	set( MKLROOT $ENV{MKLROOT} )
endif()

# Set the names of the static libraries.
set( INT_LIB_STATIC "libmkl_intel_ilp64.a" )
set( SEQ_LIB_STATIC "libmkl_sequential.a" )
set( THR_LIB_STATIC "libmkl_intel_thread.a" )
set( COR_LIB_STATIC "libmkl_core.a" )

# Set the names of the shared libraries.
set( INT_LIB "mkl_intel_ilp64" )
set( SEQ_LIB "mkl_sequential" )
set( THR_LIB "mkl_intel_thread" )
set( COR_LIB "mkl_core" )

# Test whether we need to add the ABI flag of -m64.
if (NOT DEFINED ENV{CRAY_PRGENVPGI} AND
	NOT DEFINED ENV{CRAY_PRGENVGNU} AND
	NOT DEFINED ENV{CRAY_PRGENVCRAY} AND
	NOT DEFINED ENV{CRAY_PRGENVINTEL})
	set(ABI "-m64")
endif()

# Set the variables we want to set here.
set( MKL_LIBRARIES )

find_path(MKL_INCLUDE_DIR 
	NAMES mkl.h
	HINTS ${MKLROOT}/include
)

#
# Check for OpenMP
#
find_package( OpenMP )

#if( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel" )
	set( ADD_MKL_PARALLEL_OPT TRUE )
#endif()

# 
# find_mkl_library
#
# A helper macro that implements the common code when searching for the 
# different MKL libraries.
#
# This macro will create an imported target with a name of: ${_target_name_var}
# by searching for a library with name: ${_lib_name_var}.
# The path of the target (if found), is stored in ${_target_path_var}.
#
macro( find_mkl_library _target_name_var _lib_name_var _target_path_var )
	find_library( ${_target_path_var}
		NAMES ${_lib_name_var}
		PATHS ${MKLROOT}/lib/intel64
			${INTEL}/mkl/lib/intel64
		NO_DEFAULT_PATH
	)

	# If we found the library, add it as an imported target.
	if( ${_target_path_var} )
		# Set the name of the target here for reference.
		set( _target_name "MKL::${_target_name_var}" )

		add_library( ${_target_name} SHARED IMPORTED )
		set_target_properties( ${_target_name} PROPERTIES
			IMPORTED_LOCATION ${${_target_path_var}}
		)

		set( custom_compile_opts "-DMKL_ILP64" )
		if( ABI )
			set( custom_compile_opts "${custom_compile_opts} -m64" )
		endif()
		
		# Check if OpenMP is available.
		if( OpenMP_FOUND OR OPENMP_FOUND )
			set( custom_compile_opts "${custom_compile_ops} ${OpenMP_CXX_FLAGS}" )
		endif()
		
		# If we are using the Intel compiler, we might need to add "-mkl=parallel"
		if( ADD_MKL_PARALLEL_OPT )
			set( custom_compile_opts "${custom_compile_ops} ${OpenMP_CXX_FLAGS} -mkl=parallel" )
		endif()

		set_property( TARGET ${target_name} PROPERTY INTERFACE_COMPILE_OPTIONS "${custom_compile_opts}" )
		
		list( APPEND MKL_LIBRARIES ${_target_name} )
	endif()
endmacro()

find_mkl_library( mkl_interface    "${INT_LIB}"  MKL_INTERFACE_LIBRARY )
find_mkl_library( mkl_sequential   "${SEQ_LIB}"  MKL_SEQUENTIAL_LAYER_LIBRARY )
find_mkl_library( mkl_core         "${COR_LIB}"  MKL_CORE_LIBRARY )
find_mkl_library( mkl_intel_thread "${THR_LIB}"  MKL_OPENMP_THREAD_LIBRARY )

set( MKL_INCLUDE_DIRS ${MKL_INCLUDE_DIR} )

# Handle the QUIETLY and REQUIRED arguments and set MKL_FOUND to TRUE if
# all listed variables are TRUE.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args( MKL
	DEFAULT_MSG
	MKL_LIBRARIES
	MKL_INCLUDE_DIRS
	MKL_INTERFACE_LIBRARY
	MKL_CORE_LIBRARY
	MKL_OPENMP_THREAD_LIBRARY

#	MKL_SEQUENTIAL_LAYER_LIBRARY
)

mark_as_advanced( MKL_INCLUDE_DIRS
	MKL_LIBRARIES
	MKL_INTERFACE_LIBRARY
	MKL_CORE_LIBRARY
	MKL_OPENMP_THREAD_LIBRARY

#	MKL_SEQUENTIAL_LAYER_LIBRARY
)

