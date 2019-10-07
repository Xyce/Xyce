### Checks for POSIX-defined features ###

# Used with device plugins:
check_include_file_cxx( "dlfcn.h" HAVE_DLFCN_H )

# Miscellaneous constants and functions associated with a POSIX OS:
check_include_file_cxx( "unistd.h" HAVE_UNISTD_H )

# Provides high-quality random numbers:
check_cxx_symbol_exists( drand48 "cstdlib" HAVE_DRAND48 )

# Check for sys/* headers
# This is for timing of runs:
check_include_file_cxx( "sys/resource.h" HAVE_SYS_RESOURCE_H )
# Gives informations about files:
check_include_file_cxx( "sys/stat.h" HAVE_SYS_STAT_H )

# see `src/UtilityPKG/N_UTL_CheckIfValidFile.C` for more stuff about
# HAVE_SYS_STAT_H that should be here.

### Check for Windows features ###

check_include_file_cxx( "Windows.h" HAVE_WINDOWS_H )

# The Windows directory commands are supposedly inside the "direct.h" header.
# See, e.g.,
# <https://docs.microsoft.com/en-us/cpp/c-runtime-library/reference/chdir-wchdir>
#
# Look in the `initializeEarly` function in `src/CircuitPKG/N_CIR_Xyce.C`. The
# logic indicates that, on a Windows system, the Windows directory commands
# (_getcwd and _chdir) are available in unistd.h.
#
# It almost has to be the case that direct.h gets included by unistd.h, somehow.
#
# As an aside, note that, according to this link,
# <https://docs.microsoft.com/en-us/cpp/c-runtime-library/reference/getcwd?view=vs-2019>
# the normal getcwd that would be accessed via unistd.h is deprecated.
#
# At any rate, this follows the proper logic for what is in the code right now,
# and it works.  The commented out part (for "direct.h") might make more sense
# in the future, if things go awry.
#
#check_include_file_cxx( "direct.h" HAVE_DIRECT_H )
#if (HAVE_DIRECT_H)
#     check_cxx_symbol_exists ( _chdir "direct.h" HAVE_WIN_CHDIR )
#     check_cxx_symbol_exists ( _getcwd "direct.h" HAVE_WIN_GETCWD )
if (HAVE_UNISTD_H)
     check_cxx_symbol_exists ( _chdir "unistd.h" HAVE_WIN_CHDIR )
     check_cxx_symbol_exists ( _getcwd "unistd.h" HAVE_WIN_GETCWD )
     if ( HAVE_WIN_CHDIR AND HAVE_WIN_GETCWD )
          set (HAVE_WIN_DIRCOMMANDS TRUE)
     endif ()
endif ()
if(NOT HAVE_UNISTD_H)
     set(YY_NO_UNISTD_H TRUE)
endif()


# Finally we must have POSIX or Windows
# This logic is also in `src/UtilityPKG/N_UTL_CPUTime.C`.  It should probably
# be removed from there.  The error message should probably be made more
# meaningful, too.
if (NOT (HAVE_WINDOWS_H OR (HAVE_UNISTD_H AND HAVE_SYS_RESOURCE_H)))
     message(FATAL_ERROR "Neither Windows nor POSIX features are available."
          "Unable to define cpu_time() for an unknown OS.")
endif ()
if ( WIN32 )
   set (_CRT_SECURE_NO_WARNINGS TRUE)
endif()
