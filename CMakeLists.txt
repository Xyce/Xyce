
# ----------------------------------------------------------------------------
# CMake Setup          
# ----------------------------------------------------------------------------

# -- configure the CMake environment -----------------------------------------

# set CMake minimum version
cmake_minimum_required ( VERSION 2.6 )


# set makefile verbosity level 
set ( CMAKE_VERBOSE_MAKEFILE ON )

# accept new policy while it is in transition
# this policy "now prefer to link the check executable 
# to the libraries listed in the CMAKE_REQUIRED_LIBRARIES variable."
# rather than just check that the library exists.'
# cmake_policy(SET CMP0075 NEW)
# can't do this yet as CMake is old on some of our platforms.

# load CMake macros
include ( CheckIncludeFileCXX )
include ( CheckFortranFunctionExists )
include ( CheckFunctionExists )
include ( CheckLibraryExists )
include ( CheckSymbolExists )
include ( CheckCXXSourceCompiles )



# Xyce library macro 
macro ( XYCE_FIND_LIBRARY LIBNAME )
  message ( STATUS "Looking for library ${LIBNAME}" )
  find_library ( LIBRARY_${LIBNAME} 
                 NAMES ${LIBNAME} lib${LIBNAME}
                 PATHS ${CMAKE_LIBRARY_PATH} )

  if ( LIBRARY_${LIBNAME} )
    message ( STATUS "Looking for library ${LIBNAME} - found" )

    # DEBUG:  /link flag broken on win32; use full paths
    if ( CMAKE_HOST_WIN32 )
      set ( ARCH_LIBRARIES ${LIBRARY_${LIBNAME}} ${ARCH_LIBRARIES} )
    else ( CMAKE_HOST_WIN32 )
      set ( ARCH_LIBRARIES ${LIBNAME} ${ARCH_LIBRARIES} )
    endif ( CMAKE_HOST_WIN32 )
  else ( LIBRARY_${LIBNAME} )
    message ( FATAL_ERROR "Looking for library ${LIBNAME} - not found" )
  endif ( LIBRARY_${LIBNAME} )
endmacro ( XYCE_FIND_LIBRARY LIBNAME )

# Xyce find library from list macro.
# Any arguments passed to this macro are assumed to be 
# library names for which to search.  They are processed 
# with the variable ARGV which is all the unnamed values
# passed in.  When a working library is found, the search
# is exited rather than continue.  
macro ( XYCE_FIND_LIBRARY_FROM_LIST  )
  set (LIB_FOUND 0)
  foreach( LIBNAME ${ARGV})
    if( ${LIB_FOUND} EQUAL 0)
      # only search if we haven't found a library yet.
      message ( STATUS "Looking for library ${LIBNAME}" )
      find_library ( LIBRARY_${LIBNAME} 
                   NAMES ${LIBNAME} lib${LIBNAME}
                   PATHS ${CMAKE_LIBRARY_PATH} )
                   
      if ( LIBRARY_${LIBNAME} )
        message ( STATUS "Looking for library ${LIBNAME} - found" )
        set (LIB_FOUND 1)

        # DEBUG:  /link flag broken on win32; use full paths
        if ( CMAKE_HOST_WIN32 )
          set ( ARCH_LIBRARIES ${LIBRARY_${LIBNAME}} ${ARCH_LIBRARIES} )
        else ( CMAKE_HOST_WIN32 )
          set ( ARCH_LIBRARIES ${LIBNAME} ${ARCH_LIBRARIES} )
        endif ( CMAKE_HOST_WIN32 )
      endif ( LIBRARY_${LIBNAME} )
    endif( ${LIB_FOUND} EQUAL 0)
  endforeach( LIBNAME ${ARGV})
  
  if ( $(LIB_FOUND) EQUAL 0 )
    message ( FATAL_ERROR "Looking for library ${ARGV} - not found" )
  endif ( $(LIB_FOUND) EQUAL 0 )
endmacro ( XYCE_FIND_LIBRARY_FROM_LIST )

macro ( XYCE_FIND_OPTIONAL_LIBRARY LIBNAME )
  message ( STATUS "Looking for library ${LIBNAME}" )
  find_library ( LIBRARY_${LIBNAME} 
                 NAMES ${LIBNAME} lib${LIBNAME}
                 PATHS ${CMAKE_LIBRARY_PATH} )

  if ( LIBRARY_${LIBNAME} )
    message ( STATUS "Looking for library ${LIBNAME} - found" )

    # DEBUG:  /link flag broken on win32; use full paths
    if ( CMAKE_HOST_WIN32 )
      set ( ARCH_LIBRARIES ${LIBRARY_${LIBNAME}} ${ARCH_LIBRARIES} )
    else ( CMAKE_HOST_WIN32 )
      set ( ARCH_LIBRARIES ${LIBNAME} ${ARCH_LIBRARIES} )
    endif ( CMAKE_HOST_WIN32 )
  endif ( LIBRARY_${LIBNAME} )

endmacro ( XYCE_FIND_OPTIONAL_LIBRARY LIBNAME )

set(Xyce_CMAKE_OPTIONS "" CACHE FILEPATH "List of option files to load")

# load configuration files;  cmd line flags override;  no trap for missing file
foreach (CONFIG_OPTS_FILE ${Xyce_CMAKE_OPTIONS})
  message ( STATUS "Reading in configuration options from ${CONFIG_OPTS_FILE} ..." )
  include ( ${CONFIG_OPTS_FILE} )
endforeach (CONFIG_OPTS_FILE ${Xyce_CMAKE_OPTIONS})

# All of the default options are gathered in one file.
message ( STATUS "Loading ${CMAKE_HOME_DIRECTORY}/defaultConfig.cmake" )
include ( "${CMAKE_HOME_DIRECTORY}/defaultConfig.cmake" )



# -- process build option dependencies; notify user of configuration changes--

# check for consistent library and binary flags                         

if ( Xyce_STATIC_BINARY AND CMAKE_HOST_APPLE )
  message ( STATUS "STATIC_BINARY is not supported on OSX.  Disabling static linkage." )
  set ( Xyce_STATIC_BINARY OFF )
endif ( Xyce_STATIC_BINARY AND CMAKE_HOST_APPLE )

if ( Xyce_STATIC_BINARY AND NOT Xyce_ENABLE_STATIC ) 
  message ( STATUS "STATIC_BINARY requires ENABLE_STATIC.  Enabling static libaries." )
  set ( Xyce_ENABLE_STATIC ON ) 
endif ( Xyce_STATIC_BINARY AND NOT Xyce_ENABLE_STATIC ) 

if ( Xyce_STATIC_BINARY AND Xyce_DISABLE_BINARY )
  message ( STATUS "STATIC_BINARY conflicts with DISABLE_BINARY.  Enabling binary." )
  set ( Xyce_DISABLE_BINARY OFF )
endif ( Xyce_STATIC_BINARY AND Xyce_DISABLE_BINARY )

if ( Xyce_DISABLE_BINARY AND NOT Xyce_ENABLE_SHARED AND NOT Xyce_ENABLE_STATIC )
  message ( STATUS "No targets are enabled.  Resetting defaults." )
  set ( Xyce_DISABLE_BINARY OFF )
  set ( Xyce_ENABLE_STATIC ON )
endif ( Xyce_DISABLE_BINARY AND NOT Xyce_ENABLE_SHARED AND NOT Xyce_ENABLE_STATIC )

set ( Xyce_BELOS ON )
set ( Xyce_USE_ISORROPIA ON )
set ( Xyce_USE_ZOLTAN ON  )

# disable options that conflict with MPI 
if ( Xyce_SUPERLUDIST AND NOT Xyce_PARALLEL_MPI )
  message ( STATUS "SUPERLUDIST requires PARALLEL_MPI.  Disabling SUPERLUDIST" )      
  set ( Xyce_SUPERLUDIST OFF )
endif ( Xyce_SUPERLUDIST AND NOT Xyce_PARALLEL_MPI )


# disable Sandia models if the source directory is missing
if ( Xyce_RAD_MODELS AND NOT EXISTS ${Xyce_SOURCE_DIR}/src/DeviceModelPKG/SandiaModels ) 
  message ( STATUS "SandiaModels directory cannot be found.  Disabling RAD_MODELS" )      
  set ( Xyce_RAD_MODELS OFF )  
endif ( Xyce_RAD_MODELS AND NOT EXISTS ${Xyce_SOURCE_DIR}/src/DeviceModelPKG/SandiaModels ) 

# disable non-free models if the source directory is missing
if ( Xyce_NONFREE_MODELS AND NOT EXISTS ${Xyce_SOURCE_DIR}/src/DeviceModelPKG/Xyce_NonFree ) 
  message ( STATUS "Xyce_NonFree directory cannot be found.  Disabling NONFREE_MODELS" )      
  set ( Xyce_NONFREE_MODELS OFF )  
endif ( Xyce_NONFREE_MODELS AND NOT EXISTS ${Xyce_SOURCE_DIR}/src/DeviceModelPKG/Xyce_NonFree )

if ( NOT EXISTS ${Xyce_SOURCE_DIR}/.git )
	set ( Xyce_REPO_ACCESS OFF )
	message ( STATUS ".git not found.  Disabling git SHA-1 timestamp." )
else ( NOT EXISTS ${Xyce_SOURCE_DIR}/.git )
	set ( Xyce_REPO_ACCESS ON )
	message ( STATUS ".git found.  Enabling git SHA-1 timestamp." )
endif ( NOT EXISTS ${Xyce_SOURCE_DIR}/.git )


# Set up RPATH for shared libraries based on link path
# and location of "lib" directory in install path
SET( CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# -- prepare CMake environment for system inspection -------------------------

# enable Fortran test   
if ( Xyce_FORTRAN_TEST )
  enable_language ( Fortran )
endif ( Xyce_FORTRAN_TEST )

# Only if Xyce_ARCHDIR is defined...
if ( DEFINED Xyce_ARCHDIR )
  # inspect ARCHDIR
  if (    NOT EXISTS ${Xyce_ARCHDIR}
       OR NOT EXISTS ${Xyce_ARCHDIR}/include
       OR NOT EXISTS ${Xyce_ARCHDIR}/lib )
    message ( FATAL_ERROR "${Xyce_ARCHDIR} failed inspection." )
  endif ()

  # Add it to CMake's prefix paths.
  list ( APPEND CMAKE_PREFIX_PATH ${Xyce_ARCHDIR} )

  # add common include file directories
  set ( CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH}
      ${Xyce_ARCHDIR}/include
      ${Xyce_ARCHDIR}/include/trilinos )


  # add common library directories
  set ( CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH}
      ${Xyce_ARCHDIR}/lib
      ${Xyce_BINARY_DIR}/src )
endif ( DEFINED Xyce_ARCHDIR )


# set Dakota support header and library paths
if ( Xyce_Dakota AND DEFINED Xyce_ARCHDIR )

  if ( EXISTS ${Xyce_ARCHDIR}/dakota/include )
    set ( CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} ${Xyce_ARCHDIR}/dakota/include )
  endif ( EXISTS ${Xyce_ARCHDIR}/dakota/include )

  if ( EXISTS ${Xyce_ARCHDIR}/reference/lib ) 
    set ( CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} ${Xyce_ARCHDIR}/reference/lib )
  endif ( EXISTS ${Xyce_ARCHDIR}/reference/lib ) 

  if ( EXISTS ${Xyce_ARCHDIR}/dakota/lib )
    set ( CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} ${Xyce_ARCHDIR}/dakota/lib )
  endif ( EXISTS ${Xyce_ARCHDIR}/dakota/lib )

endif ( Xyce_Dakota AND DEFINED Xyce_ARCHDIR )

# set header and linker flags
set (CMAKE_REQUIRED_INCLUDES ${CMAKE_INCLUDE_PATH} )
include_directories( ${CMAKE_INCLUDE_PATH} )
link_directories( ${CMAKE_LIBRARY_PATH} )


# get compiler vendor
get_property ( COMPILER_VENDOR VARIABLE PROPERTY CMAKE_CXX_COMPILER_ID ) 


# set Windows (native & cygwin) specific flags
if ( CMAKE_HOST_WIN32 )

  set ( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -nologo" )
  set ( WINDOWS_BUILD ON )

  # disable "safe" features of MSVS
  add_definitions ( -D_SECURE_SCL=0 )

  # Make "cmath" and math.h include constants like M_PI
  add_definitions ( -D_USE_MATH_DEFINES=1 )


  ###### Eric R does this in his local version with a comment "DEBUG ELR: 
  ######  broken static builds!

  add_definitions ( -MD )

  # replace stl components with their MSVS "safe" counterparts
  ## add_definitions ( -D_CRT_SECURE_CPP_OVERLOAD_STANDARD_NAMES=1 -D_CRT_SECURE_CPP_OVERLOAD_STANDARD_NAMES_COUNT=1 )

  set ( Xyce_RESTART_NOPACK ON )

  # remove linker conflicts
  ## set ( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /NODEFAULTLIB:libcmt.lib /NODEFAULTLIB:libcmtd.lib" )
  #set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /NODEFAULTLIB:libmmd.lib /DEFAULTLIB:libmmt.lib")

endif ( CMAKE_HOST_WIN32 )




# -- system inspection --------------------------------------------------------

# -- search for required programs

if ( Xyce_REACTION_PARSER )
  find_package ( BISON 2.3 )
  find_package ( FLEX )
  if ( NOT FLEX_FOUND OR NOT BISON_FOUND )
    message ( STATUS "Could not find bison/flex.  Disabling REACTION_PARSER" )
    set ( Xyce_REACTION_PARSER OFF )  
  endif ( NOT FLEX_FOUND OR NOT BISON_FOUND )
endif ( Xyce_REACTION_PARSER )


# -- search for libraries

# base math library (no usually available under windows)
find_library( LIB_M m )

if ( Xyce_USE_FFTW )
  XYCE_FIND_LIBRARY ( fftw3 )
endif ( Xyce_USE_FFTW )

if ( Xyce_USE_HDF5 )
  XYCE_FIND_LIBRARY ( hdf5 sz z )
  # need to update this so that sz and z are only used if needed 
  # as we don't use those components
  # and if hdf5 is not compiled with parallel support, then don't
  # activate it for mpi parallel builds
  XYCE_FIND_LIBRARY ( sz )
  XYCE_FIND_LIBRARY ( z )
endif ( Xyce_USE_HDF5 )


if ( Xyce_USE_CURL )
  XYCE_FIND_LIBRARY ( curl )
  if ( DEFINED  TRACKING_URL)
     add_definitions( -DXyce_TRACKING_URL="${TRACKING_URL}" )
  else (DEFINED TRACKING_URL)
     add_definitions( -DXyce_TRACKING_URL="http://localhost/" )
  endif (DEFINED TRACKING_URL)
endif ( Xyce_USE_CURL )


# add lapack / blas and associated dependencies
if ( COMPILER_VENDOR STREQUAL "Intel" )

  # use the pure linking mesthod from Intel 10+

  find_library ( LIB_MKL_INTERFACE   
                 NAMES mkl_intel mkl_intel_lp64 mkl_intel_c
                 PATHS ENV "LIBRARY_PATH" ENV "LD_LIBRARY_PATH" ENV "DYLD_LIBRARY_PATH" )

  find_library ( LIB_MKL_THREADING         
                 NAMES mkl_sequential
                 PATHS ENV "LIBRARY_PATH" ENV "LD_LIBRARY_PATH" ENV "DYLD_LIBRARY_PATH" )

  find_library ( LIB_MKL_COMPUTATION       
                 NAMES mkl_core                        
                 PATHS ENV "LIBRARY_PATH" ENV "LD_LIBRARY_PATH" ENV "DYLD_LIBRARY_PATH" )

  if ( Xyce_Dakota )
    find_library ( LIB_FORTRAN               
                   NAMES ifcore                           
                   PATHS ENV "LIBRARY_PATH" ENV "LD_LIBRARY_PATH" ENV "DYLD_LIBRARY_PATH" )
  endif ( Xyce_Dakota )

  set ( MKL_LIBRARIES 
        ${LIB_MKL_INTERFACE} ${LIB_MKL_THREADING} ${LIB_MKL_COMPUTATION}
        ${LIB_MKL_INTERFACE} ${LIB_MKL_THREADING} ${LIB_MKL_COMPUTATION}
        ${LIB_FORTRAN} )


#if ( CMAKE_HOST_APPLE AND Xyce_PARALLEL_MPI )
#  find_library ( LIBIRC
#                 NAMES irc
#                 PATHS $ENV{INTEL_PATH}/$ENV{INTEL_COM_VER}/lib )
#endif ( CMAKE_HOST_APPLE AND Xyce_PARALLEL_MPI )


  # Windows (native & cygwin) and Open MPI do not require additional threading support
  if ( NOT Xyce_PARALLEL_MPI AND NOT CMAKE_HOST_WIN32 )
    find_library ( LIB_PTHREAD               
                   NAMES pthread                         
                   PATHS ENV "LIBRARY_PATH" ENV "LD_LIBRARY_PATH" ENV "DYLD_LIBRARY_PATH" )
  endif ( NOT Xyce_PARALLEL_MPI AND NOT CMAKE_HOST_WIN32 )


  # collect all libraries
  set ( EXTERNAL_LIBRARIES 
        ${MKL_LIBRARIES} ${LIB_PTHREAD} ${LIBIRC} )



else ( COMPILER_VENDOR STREQUAL "Intel" )

  find_library ( LIB_LAPACK lapack )
  find_library ( LIB_BLAS blas )

  # This stuff is unreliable and should be replaced by proper probing of
  # the fortran link requirements.  
  # Autoconf will run the Fortran compiler with options to display what 
  # libraries it needs, and add those to the link line.  CMake won't do
  # that automatically unless your project contains Fortran files.
  # This mess tries to add -lgfortran for specific system combinations, but
  # it gets it wrong for some (e.g. FreeBSD with clang).  Needed to add
  # a similar wrong block to get mingw builds to work, too.
  # DEBUG:  fix for peculiar CMake behavior for static linked binary and FBSD
  if ( Xyce_STATIC_BINARY OR CMAKE_SYSTEM_NAME STREQUAL "FreeBSD" )
    set ( LIB_FORTRAN "gfortran" )
  endif ( Xyce_STATIC_BINARY OR CMAKE_SYSTEM_NAME STREQUAL "FreeBSD" )

  if ( CMAKE_SYSTEM_NAME STREQUAL "Windows" AND COMPILER_VENDOR STREQUAL "GNU" )
    set ( LIB_FORTRAN "gfortran" )

  # as long as we've got this mingw-specific conditional, make sure to set
  # this, too:
    set ( Xyce_RESTART_NOPACK ON )

  endif ( CMAKE_SYSTEM_NAME STREQUAL "Windows" AND COMPILER_VENDOR STREQUAL "GNU" )

  # collect all libraries
  set ( EXTERNAL_LIBRARIES 
        ${LIB_LAPACK} ${LIB_BLAS} ${LIB_FORTRAN} )

endif ( COMPILER_VENDOR STREQUAL "Intel" )


if ( Xyce_PARALLEL_MPI )

  # expecting openmpi to be present in /usr/local/* or /usr/local/mpi/*
  message ( STATUS "Looking for library mpi" )                 
  file ( GLOB OPENMPI_SEARCH_PATH "/usr/local/openmpi*" "/usr/local/mpi/openmpi*" )
  find_library ( LIB_MPI 
                 NAMES mpi 
                 PATHS ${OPENMPI_SEARCH_PATH}
                 PATH_SUFFIXES lib lib/openmpi 
                 NO_DEFAULT_PATH )
  find_path ( INC_MPI 
                 NAMES mpi mpi.h
                 PATHS ${OPENMPI_SEARCH_PATH}
                 PATH_SUFFIXES include         
                 NO_DEFAULT_PATH )

  if ( NOT LIB_MPI )
    message ( STATUS "Looking for library mpi - not found" )                 
  else ( NOT LIB_MPI )
    message ( STATUS "Looking for library mpi - found" )                 
    set ( ARCH_LIBRARIES ${LIB_MPI} ${ARCH_LIBRARIES} )
    include_directories ( ${INC_MPI} )
  endif ( NOT LIB_MPI )

  if ( Xyce_USE_PARMETIS )
    XYCE_FIND_LIBRARY ( metis ) 
    XYCE_FIND_LIBRARY ( parmetis )
  endif ( Xyce_USE_PARMETIS )
 
  if ( Xyce_USE_ZOLTAN OR Xyce_USE_ISORROPIA )
    if ( DEFINED Xyce_ARCHDIR AND EXISTS ${Xyce_ARCHDIR}/include/zoltan )
      include_directories ( ${Xyce_ARCHDIR}/include/zoltan )
    endif ( DEFINED Xyce_ARCHDIR AND EXISTS ${Xyce_ARCHDIR}/include/zoltan )
    XYCE_FIND_LIBRARY ( zoltan )
  endif ( Xyce_USE_ZOLTAN OR Xyce_USE_ISORROPIA )

  add_definitions ( -DXyce_PARALLEL_MPI )
endif ( Xyce_PARALLEL_MPI )

# This check for the Amesos_Umfpack.h header needs to be performed BEFORE
# the check for umfpack and amd libraries.

find_path(HAVE_AMESOS_UMFPACK_H Amesos_Umfpack.h)
if (NOT HAVE_AMESOS_UMFPACK_H)
   message ( STATUS "Xyce_UMFPACK was set on, but Amesos_Umfpack.h not found.  Disabling Xyce_UMFPACK.")
   set (Xyce_UMFPACK OFF)
else (NOT HAVE_AMESOS_UMFPACK_H)
   message ( STATUS "Amesos_Umfpack.h found.  Enabling Xyce_UMFPACK.")
   set (Xyce_UMFPACK ON)
endif (NOT HAVE_AMESOS_UMFPACK_H) 

if ( Xyce_AMD OR Xyce_UMFPACK )
  XYCE_FIND_LIBRARY ( amd )
endif ( Xyce_AMD OR Xyce_UMFPACK )

# We no longer support Trilinos 11.0 in CMake. Upgrade to 11.2.X or later.
# XYCE_FIND_LIBRARY ( teuchos )
# The old single teuchos library of <=11.0 has been replaced by these 5:
XYCE_FIND_LIBRARY ( teuchoscore )
XYCE_FIND_LIBRARY ( teuchosremainder )
XYCE_FIND_LIBRARY ( teuchosnumerics )
XYCE_FIND_LIBRARY ( teuchoscomm )
XYCE_FIND_LIBRARY ( teuchosparameterlist )

XYCE_FIND_LIBRARY ( triutils )
XYCE_FIND_LIBRARY ( epetra )
XYCE_FIND_LIBRARY ( epetraext )

if ( Xyce_PARALLEL_MPI AND Xyce_USE_ISORROPIA )
  XYCE_FIND_LIBRARY ( isorropia ) 
endif ( Xyce_PARALLEL_MPI AND Xyce_USE_ISORROPIA )

if ( Xyce_SUPERLUDIST )
  XYCE_FIND_LIBRARY ( superludist )
endif ( Xyce_SUPERLUDIST )

if ( Xyce_UMFPACK )
  XYCE_FIND_OPTIONAL_LIBRARY (suitesparseconfig)
  XYCE_FIND_OPTIONAL_LIBRARY ( colamd )
  XYCE_FIND_OPTIONAL_LIBRARY ( cholmod )
  XYCE_FIND_LIBRARY ( umfpack )
endif ( Xyce_UMFPACK )

if ( Xyce_SHYLU )
  XYCE_FIND_LIBRARY ( shylu )
endif ( Xyce_SHYLU )

if ( Xyce_AMESOS2 )
  XYCE_FIND_LIBRARY ( amesos2 )
endif ( Xyce_AMESOS2 )

if ( Xyce_SUPERLU )
  XYCE_FIND_LIBRARY ( superlu )
endif ( Xyce_SUPERLU )

# This is absolutely NOT correct, but it'll do for now.  Trilinos
# 12.12 and later have moved some of the functions that amesos needs
# into trilinosss.  If this library exists, assume we need it, and if
# it doesn't assume we're using an older version.  This assumption
# will be WRONG if someone has a busted Trilinos install and will
# result in a link failure at the end, but it is close enough if the 
# install isn't busted.
#
# This would be better if we were just picking up amesos' dependencies
# properly. That would require doing this CMakeLists correctly (with
# FIND_PACKAGE instead of FIND_LIBRARY).

XYCE_FIND_OPTIONAL_LIBRARY ( trilinosss )
XYCE_FIND_LIBRARY ( amesos )
XYCE_FIND_LIBRARY ( trilinoscouplings )

if ( Xyce_ML )
  XYCE_FIND_LIBRARY ( ml )
endif ( Xyce_ML )

XYCE_FIND_LIBRARY ( ifpack )
XYCE_FIND_LIBRARY ( aztecoo )

if ( Xyce_BELOS )
  XYCE_FIND_LIBRARY ( belos )
  XYCE_FIND_LIBRARY ( belosepetra )
endif ( Xyce_BELOS )

XYCE_FIND_LIBRARY ( nox )
XYCE_FIND_LIBRARY ( loca )

if (CMAKE_HOST_WIN32)
  XYCE_FIND_LIBRARY(advapi32)
  XYCE_FIND_LIBRARY(user32)
  XYCE_FIND_LIBRARY(ws2_32)
  XYCE_FIND_LIBRARY(wldap32)
endif (CMAKE_HOST_WIN32)

# check for the dakota library if the user wants to link to it
# this library list is for Dakota 5.2.  We need a more flexible 
# way to get Dakota's dependent libs automatically.  The Dakota
# team is working on this for Dakota 5.2+
if ( Xyce_Dakota  ) 
  XYCE_FIND_LIBRARY ( dakota_src )
  XYCE_FIND_LIBRARY ( nidr )
  XYCE_FIND_LIBRARY ( pecos ) 
  XYCE_FIND_LIBRARY ( pecos_src )
  XYCE_FIND_LIBRARY ( lhs )
  XYCE_FIND_LIBRARY ( mods )
  XYCE_FIND_LIBRARY ( mod )
  XYCE_FIND_LIBRARY ( dfftpack )
  XYCE_FIND_LIBRARY ( sparsegrid )
  XYCE_FIND_LIBRARY ( surfpack )
  XYCE_FIND_LIBRARY ( surfpack )
  XYCE_FIND_LIBRARY ( utilib )
  XYCE_FIND_LIBRARY ( colin )
  XYCE_FIND_LIBRARY ( interfaces )
  XYCE_FIND_LIBRARY ( scolib )
  XYCE_FIND_LIBRARY ( 3po )
  XYCE_FIND_LIBRARY ( pebbl )
  XYCE_FIND_LIBRARY ( tinyxml )
  XYCE_FIND_LIBRARY ( conmin )
  XYCE_FIND_LIBRARY ( dace )
  XYCE_FIND_LIBRARY ( analyzer )
  XYCE_FIND_LIBRARY ( random )
  XYCE_FIND_LIBRARY ( sampling )
  XYCE_FIND_LIBRARY ( bose )
  XYCE_FIND_LIBRARY ( fsudace )
  XYCE_FIND_LIBRARY ( hopspack )
  XYCE_FIND_LIBRARY ( jega )
  XYCE_FIND_LIBRARY ( jega_fe )
  XYCE_FIND_LIBRARY ( moga )
  XYCE_FIND_LIBRARY ( soga )
  XYCE_FIND_LIBRARY ( eutils )
  XYCE_FIND_LIBRARY ( utilities )
  XYCE_FIND_LIBRARY ( ncsuopt )
  XYCE_FIND_LIBRARY ( nlpql )
  XYCE_FIND_LIBRARY ( cport ) 
  XYCE_FIND_LIBRARY ( npsol )
  XYCE_FIND_LIBRARY ( optpp )
  XYCE_FIND_LIBRARY ( psuade )
  XYCE_FIND_LIBRARY ( amplsolver )
  
  XYCE_FIND_LIBRARY ( boost_regex ) 
  XYCE_FIND_LIBRARY ( boost_signals ) 
  XYCE_FIND_LIBRARY ( boost_filesystem )
  XYCE_FIND_LIBRARY ( boost_system )
              
endif ( Xyce_Dakota  )

# -- search for header files

# system headers
CHECK_INCLUDE_FILE_CXX ( Windows.h HAVE_WINDOWS_H )
CHECK_INCLUDE_FILE_CXX ( fftw3.h HAVE_FFTW3_H )
CHECK_INCLUDE_FILE_CXX ( float.h HAVE_FLOAT_H )
CHECK_INCLUDE_FILE_CXX ( hdf5.h HAVE_HDF5_H )
CHECK_INCLUDE_FILE_CXX ( limits.h HAVE_LIMITS_H )
CHECK_INCLUDE_FILE_CXX ( mathimf.h HAVE_MATHIMF_H )
CHECK_INCLUDE_FILE_CXX ( unistd.h HAVE_UNISTD_H )
CHECK_INCLUDE_FILE_CXX ( values.h HAVE_VALUES_H )
CHECK_INCLUDE_FILE_CXX ( unordered_map HAVE_UNORDERED_MAP)
CHECK_INCLUDE_FILE_CXX ( tr1/unordered_map HAVE_TR1_UNORDERED_MAP)
CHECK_INCLUDE_FILE_CXX ( unordered_set HAVE_UNORDERED_SET)
CHECK_INCLUDE_FILE_CXX ( tr1/unordered_set HAVE_TR1_UNORDERED_SET)
CHECK_INCLUDE_FILE_CXX ( functional HAVE_FUNCTIONAL)
CHECK_INCLUDE_FILE_CXX ( tr1/functional HAVE_TR1_FUNCTIONAL)
CHECK_INCLUDE_FILE_CXX ( sys/resource.h HAVE_SYS_RESOURCE_H )
CHECK_INCLUDE_FILE_CXX ( sys/stat.h HAVE_SYS_STAT_H )
CHECK_INCLUDE_FILE_CXX ( curl/curl.h HAVE_CURL_CURL_H )

# state information
CHECK_INCLUDE_FILE_CXX ( unistd.h HAVE_UNISTD_H )
CHECK_INCLUDE_FILE_CXX ( malloc.h HAVE_MALLOC_H )
CHECK_INCLUDE_FILE_CXX ( pwd.h HAVE_PWD_H )
CHECK_INCLUDE_FILE_CXX ( sys/utsname.h HAVE_SYS_UTSNAME_H )
CHECK_FUNCTION_EXISTS( mallinfo HAVE_MALLINFO )
CHECK_FUNCTION_EXISTS( getpwuid HAVE_GETPWUID )
CHECK_FUNCTION_EXISTS( gethostname HAVE_GETHOSTNAME )
CHECK_FUNCTION_EXISTS( getdomainname HAVE_GETDOMAINNAME )
CHECK_FUNCTION_EXISTS( uname HAVE_UNAME )
CHECK_FUNCTION_EXISTS( erf HAVE_ERF )
CHECK_FUNCTION_EXISTS( erfc HAVE_ERFC )


# Trilinos headers
CHECK_INCLUDE_FILE_CXX ( Epetra_MultiVector.h HAVE_EPETRA_MULTIVECTOR_H )
#Don't do :
#CHECK_INCLUDE_FILE_CXX ( Ifpack_CrsRiluk.h HAVE_IFPACK_CRSRILUK_H )
# as it tries to link a test program, which fails due to dependence on teuchos
# libraries.  check_include_file_cxx can't be forced to add libraries to its
# link line!
find_path ( HAVE_IFPACK_CRSRILUK_H Ifpack_CrsRiluk.h )
find_path ( HAVE_AMESOS_KLU_H Amesos_Klu.h )
find_path ( HAVE_LOCA_PARAMETER_VECTOR_H LOCA_Parameter_Vector.H  )
find_path ( HAVE_NOX_ABSTRACT_VECTOR_H NOX_Abstract_Vector.H )

#if ( Xyce_PARALLEL_MPI AND Xyce_USE_ZOLTAN AND NOT Xyce_USE_HMETIS )
  find_path ( HAVE_ZOLTAN_QUERYOBJECT_H Zoltan_QueryObject.h )
#endif ( Xyce_PARALLEL_MPI AND Xyce_USE_ZOLTAN AND NOT Xyce_USE_HMETIS )

#if ( Xyce_SUPERLU )
   find_path(HAVE_AMESOS_SUPERLU_H Amesos_Superlu.h)
#endif ( Xyce_SUPERLU )

#if ( Xyce_Dakota )
  find_path ( HAVE_DAKOTA_STRATEGY_H DakotaStrategy.H ) 
#endif ( Xyce_Dakota )


# -- search for symbols

CHECK_FUNCTION_EXISTS ( strcasecmp HAVE_STRCASECMP )

#    set ( CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} )
if( HAVE_FLOAT_H )
  CHECK_FUNCTION_EXISTS( _isnan FUNCTION__ISNAN )
  CHECK_FUNCTION_EXISTS( _finite FUNCTION__FINITE )
  if( FUNCTION__ISNAN AND FUNCTION__FINITE ) 
    set( HAVE__ISNAN_AND__FINITE_SUPPORT on )
  endif( FUNCTION__ISNAN AND FUNCTION__FINITE ) 
endif( HAVE_FLOAT_H )


if (NOT (FUNCTION__ISNAN AND FUNCTION__FINITE))
    set ( CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${LIB_M} )
    CHECK_FUNCTION_EXISTS ( isnan FUNCTION_ISNAN )
    CHECK_FUNCTION_EXISTS ( isinf FUNCTION_ISINF )
    if ( FUNCTION_ISNAN AND FUNCTION_ISINF )
      set ( HAVE_NAN_INF_SUPPORT ON )
    endif ( FUNCTION_ISNAN AND FUNCTION_ISINF )
endif(NOT (FUNCTION__ISNAN AND FUNCTION__FINITE))


# message ( STATUS "Warning: Your compiler doesn't support isnan() and isinf().  We will supply a default checker but it is *NOT* guaranteed to work on your platform!" )

if ( HAVE_FFTW3_H )
  CHECK_SYMBOL_EXISTS ( fftw_plan_dft_1d fftw3.h Xyce_USE_FFTW )
endif( HAVE_FFTW3_H )

if ( HAVE_HDF5_H )
  CHECK_SYMBOL_EXISTS ( H5Fclose hdf5.h Xyce_USE_HDF5 )
endif( HAVE_HDF5_H )

if ( HAVE_CURL_CURL_H )
  CHECK_SYMBOL_EXISTS ( CURL curl/curl.h Xyce_USE_CURL )
endif( HAVE_CURL_CURL_H )

# search for Intel  MKL FFT header and library compatibility
CHECK_INCLUDE_FILE_CXX ( mkl_dfti.h HAVE_MKL_DFTI_H )
if ( HAVE_MKL_DFTI_H )
  if ( CMAKE_HOST_WIN32 AND COMPILER_VENDOR STREQUAL "Intel" )
    set ( Xyce_USE_INTEL_FFT ON )
  else ( CMAKE_HOST_WIN32 AND COMPILER_VENDOR STREQUAL "Intel" )
    CHECK_SYMBOL_EXISTS ( Dfti_CreateDescriptor mkl_dfti.h Xyce_USE_INTEL_FFT )
  endif ( CMAKE_HOST_WIN32 AND COMPILER_VENDOR STREQUAL "Intel" )
endif ( HAVE_MKL_DFTI_H )

if (Xyce_USE_INTEL_FFT OR Xyce_USE_FFTW)
  set (Xyce_USE_FFT ON)
else (Xyce_USE_INTEL_FFT OR Xyce_USE_FFTW)
  set (Xyce_USE_FFT OFF)
endif (Xyce_USE_INTEL_FFT OR Xyce_USE_FFTW)



# -- search for required functions

if ( Xyce_FORTRAN_TEST )
  CHECK_FORTRAN_FUNCTION_EXISTS ( dcopy HAVE_FORTRAN_DCOPY )  
  CHECK_FORTRAN_FUNCTION_EXISTS ( dgetrs HAVE_FORTRAN_DGETRS )  
else ( Xyce_FORTRAN_TEST )
  CHECK_FUNCTION_EXISTS ( dcopy HAVE_DCOPY )  
  CHECK_FUNCTION_EXISTS ( dgetrs HAVE_DGETRS )  
endif ( Xyce_FORTRAN_TEST )

CHECK_FUNCTION_EXISTS ( sqrt HAVE_SQRT )
CHECK_FUNCTION_EXISTS ( getopts HAVE_GETOPTS )
CHECK_FUNCTION_EXISTS ( drand48 HAVE_DRAND48 )

try_compile ( HAVE_IOTA           
              ${CMAKE_BINARY_DIR}
              ${CMAKE_SOURCE_DIR}/config/cmake_check_iota.cxx )



# -- check for Windows or Unix directory commands

if ( HAVE_UNISTD_H )
      CHECK_FUNCTION_EXISTS ( _chdir HAVE_WIN_CHDIR )
      CHECK_FUNCTION_EXISTS ( _getcwd HAVE_WIN_GETCWD )
      if ( HAVE_WIN_CHDIR AND HAVE_WIN_GETCWD )
	      add_definitions ( -DHAVE_WIN_DIRCOMMANDS )
      endif ( HAVE_WIN_CHDIR AND HAVE_WIN_GETCWD )
endif ( HAVE_UNISTD_H )


# -- Check if the NOX abstract solver is collecting statistics

CHECK_CXX_SOURCE_COMPILES ( "#include <NOX_SolverStats.hpp>
      int main(void)
      {
      return 0;
      }"
      USING_NOX_SOLVERSTATS )

if ( USING_NOX_SOLVERSTATS )
      add_definitions ( -DXyce_NOX_SOLVERSTATS )
endif ( USING_NOX_SOLVERSTATS )

# -- post-processing ----------------------------------------------------------


# use the config.h files for Trilinos
set ( HAVE_CONFIG_H ON )

# flag use of CMake
set ( USING_CMAKE ON )   

# add flag to config files
set ( Xyce_NOX_LOCA_SUPPORT ON )


# set lexer flags
if ( Xyce_REACTION_PARSER ) 
  if ( NOT HAVE_UNISTD_H )
    set ( YY_NO_UNISTD_H ON ) 
  endif ( NOT HAVE_UNISTD_H )
endif ( Xyce_REACTION_PARSER )

# DEBUG:  holdover to old version string approach
# generate the build timestamp and make available from timestamp.h
if ( ${CMAKE_VERSION} VERSION_GREATER 2.8.11 )
	string ( TIMESTAMP DATETIMESTAMP %Y%m%d%H%M )
	message ( STATUS "DATETIMESTAMP set as ${DATETIMESTAMP}" )
	if ( Xyce_REPO_ACCESS )
		execute_process ( COMMAND git log --format=format:%h -1 WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY} OUTPUT_VARIABLE SHA1STAMP )
		message ( STATUS "SHA1STAMP set as ${SHA1STAMP}" )
		set ( XYCEBUILDTIMESTAMP "${DATETIMESTAMP}-g-${SHA1STAMP}" CACHE INTERNAL "MUTABLE" FORCE )
	else ( Xyce_REPO_ACCESS )
		set ( XYCEBUILDTIMESTAMP "${DATETIMESTAMP}" CACHE INTERNAL "MUTABLE" FORCE )
	endif ( Xyce_REPO_ACCESS )
	message ( STATUS "XYCEBUILDTIMESTAMP set as ${XYCEBUILDTIMESTAMP}" )
else( ${CMAKE_VERSION} VERSION_GREATER 2.8.11 )
	message ( STATUS "CMake version is less than 2.8.11.  Reverting to old timestamp generation.  This will fail if you are configuring with Windows CMake and not under Cygwin." )
	execute_process ( COMMAND date "+%Y%m%d%H%M" 
		OUTPUT_VARIABLE DATETIMESTAMP
		OUTPUT_STRIP_TRAILING_WHITESPACE )
	set ( DATETIMESTAMP "${DATETIMESTAMP}" CACHE INTERNAL "MUTABLE" FORCE )
	add_custom_target ( TIMESTAMP_H ALL
		COMMAND /usr/bin/touch ${Xyce_SOURCE_DIR}/src/UtilityPKG/timestamp.h.cmake )
	message ( STATUS "DATETIMESTAMP set as ${DATETIMESTAMP}" )
	if ( Xyce_REPO_ACCESS )
		execute_process ( COMMAND git log --format=format:%h -1 WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY} OUTPUT_VARIABLE SHA1STAMP )
		message ( STATUS "SHA1STAMP set as ${SHA1STAMP}" )
		set ( XYCEBUILDTIMESTAMP "${DATETIMESTAMP}-g-${SHA1STAMP}" CACHE INTERNAL "MUTABLE" FORCE )
	else ( Xyce_REPO_ACCESS )
		set ( XYCEBUILDTIMESTAMP "${DATETIMESTAMP}" CACHE INTERNAL "MUTABLE" FORCE )
	endif ( Xyce_REPO_ACCESS )
	message ( STATUS "XYCEBUILDTIMESTAMP set as ${XYCEBUILDTIMESTAMP}" )
endif( ${CMAKE_VERSION} VERSION_GREATER 2.8.11 )

configure_file ( ${Xyce_SOURCE_DIR}/src/UtilityPKG/timestamp.h.cmake
                 ${Xyce_BINARY_DIR}/src/UtilityPKG/timestamp.h )
include_directories ( ${Xyce_BINARY_DIR}/src/UtilityPKG )



# generate the version strings
set ( Xyce_VERSION_STRING_SHORT "${Xyce_VERSION_MAJOR}.${Xyce_VERSION_MINOR}" )

if ( DEFINED Xyce_VERSION_PATCH AND NOT (Xyce_VERSION_PATCH STREQUAL "" OR Xyce_VERSION_PATCH STREQUAL "0"))
  set ( Xyce_VERSION_STRING_SHORT "${Xyce_VERSION_STRING_SHORT}.${Xyce_VERSION_PATCH}" )
endif ( DEFINED Xyce_VERSION_PATCH AND NOT (Xyce_VERSION_PATCH STREQUAL "" OR Xyce_VERSION_PATCH STREQUAL "0"))

set ( Xyce_VERSION_STRING_LONG "${Xyce_VERSION_STRING_SHORT}" )

if ( DEFINED Xyce_VERSION_EXTRA AND NOT(Xyce_VERSION_EXTRA STREQUAL ""))
    set ( Xyce_VERSION_STRING_LONG "${Xyce_VERSION_STRING_LONG} ${Xyce_VERSION_EXTRA}" )
endif ( DEFINED Xyce_VERSION_EXTRA AND NOT( Xyce_VERSION_EXTRA STREQUAL ""))


if ( Xyce_ATHENA )
  if ( Xyce_RAD_MODELS )
    if( EXISTS ${Xyce_SOURCE_DIR}/src/DeviceModelPKG/SandiaModels/ATHENA/ )
      find_package(Boost 1.41 REQUIRED COMPONENTS)
      include_directories(${Boost_INCLUDE_DIRS})
      set ( EXTERNAL_LIBRARIES ${EXTERNAL_LIBRARIES} ${Boost_LIBRARIES})
      # MacPorts adds "-mt" to boot library names for "multi-threaded" and the 
      # regular boost library names will not be present.  So look for the -mt 
      # library names if the expected ones are not found
      add_definitions ( -DBOOST_ALL_NO_LIB )
    else ( EXISTS ${Xyce_SOURCE_DIR}/src/DeviceModelPKG/SandiaModels/ATHENA/ )
      message ( STATUS "Xyce_RAD_MODELS is enabled, but the SandiaModels/ATHENA directory was not found." )
      message ( STATUS "Try updating your version of SandiaModels." )
      message ( STATUS "Disabling Xyce_ATHENA." )
      set ( Xyce_ATHENA OFF )
    endif ( EXISTS ${Xyce_SOURCE_DIR}/src/DeviceModelPKG/SandiaModels/ATHENA/ )
  else ( Xyce_RAD_MODELS )
    set ( Xyce_ATHENA OFF )
  endif ( Xyce_RAD_MODELS )
endif ( Xyce_ATHENA )

if ( Xyce_RELEASE_NOTES_VERSION )
  message ( STATUS "Release notes version manually specified as " ${Xyce_RELEASE_NOTES_VERSION} )
else ( Xyce_RELEASE_NOTES_VERSION )
  set ( Xyce_RELEASE_NOTES_VERSION ${Xyce_VERSION_STRING_SHORT} )
  message ( STATUS "Release notes version defaulting to " ${Xyce_RELEASE_NOTES_VERSION} )
endif ( Xyce_RELEASE_NOTES_VERSION )

#
# not great but need the boost libraries that ATHENA uses 
# for the ARAEnv library to be useful
#
if( Xyce_ATHENA AND Xyce_RAD_MODELS )
  XYCE_FIND_OPTIONAL_LIBRARY( ARAEnv )
  if( LIBRARY_ARAEnv )
    # have to set the compiler to just compile, not link as the tested header file has static 
    # member functions that don't get resolved by just this check.  Could avoid this by 
    # checking for simpler header like SRdata.h 
    if ( CMAKE_HOST_WIN32 )
      set( CMAKE_REQUIRED_FLAGS "/C")
    else ( CMAKE_HOST_WIN32 )
      set( CMAKE_REQUIRED_FLAGS "-c")
    endif ( CMAKE_HOST_WIN32 )
    CHECK_INCLUDE_FILE_CXX ( SRenvironment.h HAVE_SRENVIRONMENT_H)
  
    set ( HAVE_LIBARAENV ON  )
  endif( LIBRARY_ARAEnv )
endif( Xyce_ATHENA AND Xyce_RAD_MODELS )

# -- generate source files and update the include list -----------------------

configure_file ( "${Xyce_SOURCE_DIR}/src/Xyce_config.h.cmake"
                 "${Xyce_BINARY_DIR}/src/Xyce_config.h" )

set ( XYCE_INCLUDE_DIRECTORIES ${XYCE_INCLUDE_DIRECTORIES}
      ${Xyce_BINARY_DIR}/src 
      CACHE INTERNAL "X_I" )




# -- process subdirectories --------------------------------------------------

add_subdirectory ( src )
add_subdirectory ( utils )



# -- enable dashboard and CTest ----------------------------------------------

# DEBUG: should this be in Xyce_Regression as an include() instead?
set ( DROP_METHOD "http" )
set ( DROP_SITE "s858869.sandia.gov" )
set ( DROP_LOCATION "/CDash/submit.php?project=Xyce" )
include ( CTest )
