### Checks for POSIX-defined features ###

# Used with device plugins:
check_include_file_cxx( "dlfcn.h"    HAVE_DLFCN_H )
# Provides high-quality random numbers:
check_include_file_cxx( "drand48"    HAVE_DRAND48 )
# Miscellaneous constants and functions associated with a POSIX OS:
check_include_file_cxx( "unistd.h"   HAVE_UNISTD_H )

# Check for sys/* headers
# This is for timing of runs:
check_include_file_cxx( "sys/resource.h" HAVE_SYS_RESOURCE_H )
# Gives informations about files:
check_include_file_cxx( "sys/stat.h"     HAVE_SYS_STAT_H )


### Check for Windows features ###
check_include_file_cxx( "Windows.h"  HAVE_WINDOWS_H )

# !!!! This was in the old CMakeLists.txt file; use it? !!!
#    if ( HAVE_UNISTD_H )
#         CHECK_FUNCTION_EXISTS ( _chdir HAVE_WIN_CHDIR )
#         CHECK_FUNCTION_EXISTS ( _getcwd HAVE_WIN_GETCWD )
#         if ( HAVE_WIN_CHDIR AND HAVE_WIN_GETCWD )
#              add_definitions ( -DHAVE_WIN_DIRCOMMANDS )
#         endif ( HAVE_WIN_CHDIR AND HAVE_WIN_GETCWD )
#    endif ( HAVE_UNISTD_H )


### Checks for C++11 features, and their alternatives ###

# These are C++11 standards
check_include_file_cxx( "unordered_map" HAVE_UNORDERED_MAP )
check_include_file_cxx( "unordered_set" HAVE_UNORDERED_SET )
# These were what were available prior to the above
check_include_file_cxx( "tr1/unordered_map" HAVE_TR1_UNORDERED_MAP )
check_include_file_cxx( "tr1/unordered_set" HAVE_TR1_UNORDERED_SET )
# Do the above need to be in some sort of logic?
# They are temporary, so as long as everything works, probably not.

# iota became part of "numeric" in the C++11 standard
# This doesn't seem to work  WHY?????????????!!!!!!!!!!!!!!!!  (JCV)
#CHECK_FUNCTION_EXISTS ( iota HAVE_IOTA )
# Commenting out for right now, and setting it to be true in "simple_features.cmake"
# I (JCV) also tried the following, which also does not work:
# check_cxx_symbol_exists ( iota numeric HAVE_IOTA )
# May have to use check_cxx_source_compiles, if we can't figure it out.

# erf became part of "cmath" in the C++11 standard
CHECK_FUNCTION_EXISTS ( erf HAVE_ERF )

# erfc became part of "cmath" in the C++11 standard
CHECK_FUNCTION_EXISTS ( erfc HAVE_ERFC )

# isnan and isinf functions
# these became part of "cmath" in the C++11 standard
# This one is a little messed up from the ideal; but it works, so we can keep
# doing it until we require C++11
#
# For all standard Sandia systems, they should be available in the std
# namespace; i.e., test to see if `std::isnan` and `std::isinf` are available
# Maybe use these?
#    CHECK_FUNCTION_EXISTS ( isnan HAVE_ISNAN )
#    CHECK_FUNCTION_EXISTS ( isinf HAVE_ISINF )
# If they _are not_ available, then the old CMake system effectively assumes it
# is a Windows build. In this case, the following CMake snippet is used:
#    if( HAVE_FLOAT_H )
#      CHECK_FUNCTION_EXISTS( _isnan FUNCTION__ISNAN )
#      CHECK_FUNCTION_EXISTS( _finite FUNCTION__FINITE )
#      if( FUNCTION__ISNAN AND FUNCTION__FINITE )
#        set( HAVE__ISNAN_AND__FINITE_SUPPORT on )
#      endif( FUNCTION__ISNAN AND FUNCTION__FINITE )
#    endif( HAVE_FLOAT_H )
#
# If both of those fail, then the CMake configuration should fail !!!!


