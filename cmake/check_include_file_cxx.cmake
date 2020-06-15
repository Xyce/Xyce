### Checks for POSIX-defined features ###

# Used with device plugins:
check_include_file_cxx( "dlfcn.h" HAVE_DLFCN_H )

# Miscellaneous constants and functions associated with a POSIX OS:
check_include_file_cxx( "unistd.h" HAVE_UNISTD_H )

# This is to keep flex from inserting "#include <unistd.h>" in its .h files
if(NOT HAVE_UNISTD_H)
     set(YY_NO_UNISTD_H TRUE)
endif()

# Provides high-quality random numbers:
check_cxx_symbol_exists( drand48 "cstdlib" HAVE_DRAND48 )
# Should this be changed to use the C++11 <random> library?

# Check for sys/* headers
# This is for timing of runs:
check_include_file_cxx( "sys/resource.h" HAVE_SYS_RESOURCE_H )
# Gives informations about files:
check_include_file_cxx( "sys/stat.h" HAVE_SYS_STAT_H )

# Check for headers and functions needed in stat reporting
check_include_file_cxx( "malloc.h" HAVE_MALLOC_H )
check_include_file_cxx( "pwd.h" HAVE_PWD_H )
check_include_file_cxx( "sys/utsname.h" HAVE_SYS_UTSNAME_H)
#if( HAVE_MALLOC_H )
#  set( HAVE_MALLINFO TRUE )
#else( HAVE_MALLOC_H )

check_cxx_symbol_exists(mallinfo "malloc.h" HAVE_MALLINFO)
check_cxx_symbol_exists(getpwuid "pwd.h" HAVE_GETPWUID)
check_cxx_symbol_exists(gethostname "unistd.h" HAVE_GETHOSTNAME)
check_cxx_symbol_exists(getdomainname "unistd.h" HAVE_GETDOMAINNAME)
check_cxx_symbol_exists(uname "sys/utsname.h" HAVE_UNAME)
check_cxx_symbol_exists(sysconf "unistd.h" HAVE_SYSCONF)
if(EXISTS /proc/self/stat)
  set( HAVE__PROC_SELF_STAT TRUE )
endif(EXISTS /proc/self/stat)

# see `src/UtilityPKG/N_UTL_CheckIfValidFile.C` for more stuff about
# HAVE_SYS_STAT_H that should be here.

### Check for Windows features ###

# This should be able to replace dlfcn.h, but we have not enabled that
# capability
check_include_file_cxx( "Windows.h" HAVE_WINDOWS_H )

# If compiling with a compiler for the native Windows C library, the following
# should be set to TRUE.  Otherwise, the compiler will have warnings for every
# occurence of sprintf.
# This check could be made more specific, but it's harmless to have it TRUE for
# Windows, in general.
if ( WIN32 )
   set (_CRT_SECURE_NO_WARNINGS TRUE)
endif()

#
# packed (binary) restart files are not supported under Windows.
# automatically disable under Windows
if ( WIN32 )
   set (Xyce_RESTART_NOPACK TRUE)
endif()

# Look for the Windows directory commands
# This is to support the `initializeEarly` function in
# `src/CircuitPKG/N_CIR_Xyce.C`. The actual code is called only if
# DEBUG_ALL_PROCS_SAME_WD is declared to be TRUE.
#
check_include_file_cxx( "direct.h" HAVE_DIRECT_H )
if (HAVE_DIRECT_H)
     check_cxx_symbol_exists ( _chdir "direct.h" HAVE_WIN_CHDIR )
     check_cxx_symbol_exists ( _getcwd "direct.h" HAVE_WIN_GETCWD )
     if ( HAVE_WIN_CHDIR AND HAVE_WIN_GETCWD )
          set (HAVE_WIN_DIRCOMMANDS TRUE)
     endif ()
endif ()

# We must have POSIX or Windows for time reporting
# This logic is for reporting time in `src/UtilityPKG/N_UTL_CPUTime.C`. The
# error message should probably be made more meaningful.
if (NOT (HAVE_WINDOWS_H OR (HAVE_UNISTD_H AND HAVE_SYS_RESOURCE_H)))
     message(FATAL_ERROR "Neither Windows nor POSIX features are available."
          "Unable to define cpu_time() for an unknown OS.")
endif ()
