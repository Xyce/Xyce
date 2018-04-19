#.rst:
# FindFFTW
# -------
#
# Finds the FFTW library
#
# This will define the following variables::
#
#    FFTW_FOUND - True if the FFTW library was found.
#    FFTW_VERSION - The version of FFTW that was found.
#    FFTW_INCLUDE_DIR - The include path to fftw (if found).
#    FFTW_LIBRARIES - The fftw libraries (if found).
#
include( LibFindMacros )

libfind_pkg_check_modules( FFTW_PKGCONF fftw )

# Include dir.
find_path( FFTW_INCLUDE_DIR
	NAMES fftw3.h
	PATHS ${FFTW_PKGCONF_INCLUDE_DIRS}
)

# Library itself.
find_library( FFTW_LIBRARY
	NAMES fftw3
	PATHS ${FFTW_PKGCONF_LIBRARY_DIRS}
)

# Set include dir and libraries.
set( FFTW_PROCESS_INCLUDES FFTW_INCLUDE_DIR )
set( FFTW_PROCESS_LIBS FFTW_LIBRARY )

libfind_process( FFTW )

set( FFTW_INCLUDE_DIRS ${FFTW_INCLUDE_DIR} )
set( FFTW_LIBRARIES ${FFTW_LIBRARY} )


include( FindPackageHandleStandardArgs )
find_package_handle_standard_args( FFTW
	REQUIRED_VARS
		FFTW_INCLUDE_DIRS
		FFTW_LIBRARIES
)
