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
#

find_path( FFTW_INCLUDE_DIRS fftw
	NAMES fftw3.h
)
find_library( FFTW_LIBRARIES
	NAMES fftw3
)

# Search for PkgConfig to make this easier.
#find_package(PkgConfig)
#pkg_check_modules(PC_FFTW QUIET FFTW )
#find_path( FFTW_INCLUDE_DIR
#	NAMES fftw.h
#	PATHS ${PC_FFTW_INCLUDE_DIRS}
#	PATH_SUFFIXES FFTW
#)
#find_library( FFTW_LIBRARY
#	NAMES fftw
#	PATHS ${PC_FFTW_LIBRARY_DIRS}
#)
#include(FindPackageHandleStandardArgs)
#find_package_handle_standard_args(FFTW
#	FOUND_VAR FFTW_FOUND
#	REQUIRED_VARS
#		FFTW_LIBRARY
#		FFTW_INCLUDE_DIR
#	VERSION_VAR FFTW_VERSION
#)


