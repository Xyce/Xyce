# This is to define virtual target out of some of the Third Party Softwares
# that provide clumps of variables but not nice little virtual target
# keywords, like 'fftw3'

# The below text is copied directly from Bryan Hughes' tps.cmake file on another project. It's art.
# ---
# include third party info

# ---------------------------------------------------------------------
# ***** READ THE TEXT BELOW BEFORE YOU CHANGE ANYTHING IN THIS FILE *****

# You are modifying a file that affects the ENTIRE BUILD. Before making changes, understand what you are doing!

# DO NOT add "include_directories" or "link_directories" statements to this file. Yes, I am talking to YOU.

# Try the following instead:
# * If there is a Find<yourlibrary>.cmake file, add find_package(yourlibrary REQUIRED) to this file.
# * If there is NOT a Find<yourlibrary>.cmake file, and you are are unable to find one on the internet,
#   then "set" the following variables:
#     a variable named <yourlibrary>_INCLUDE_DIR (or something similiar)
#     a variable named <yourlibrary>_LIBRARIES (or something similar)

# After you have done one of the things above, change the leaf CMakeLists.txt file inside your target directory to
# include the directories or libraries, local to your target.

#   THE WHOLE WORLD DOES NOT NEED YOUR THIRD PARTY LIBRARY'S INCLUDE AND LINKER PATHS.

# If you do not understand these instructions or are otherwise unable to follow them, ASK FOR HELP.

# If you disregard this message, I reserve the right to spray you with a hose, because you are a bad person.

# ***** READ THE TEXT ABOVE BEFORE YOU CHANGE ANYTHING IN THIS FILE *****
# ---------------------------------------------------------------------

###################
## Trilinos
###################

# A starting reference for the Trilinos CMake package is:
#    `Trilinos/commonTools/importing/README.Finding_Trilinos`
# Another starting point is:
#    `Trilinos/cmake/doc/export_system/Finding_Trilinos.txt`

# Trilinos recommends it be found BEFORE project() is called.
# One might put the search here, BUT some important (and subtle?) aspects setup
# by project() would likely be missing.

# Fix the library lists, as they contain a lot of duplicates
LIST(REVERSE Trilinos_LIBRARIES)
LIST(REMOVE_DUPLICATES Trilinos_LIBRARIES)
LIST(REVERSE Trilinos_LIBRARIES)

LIST(REVERSE Trilinos_TPL_LIBRARIES)
LIST(REMOVE_DUPLICATES Trilinos_TPL_LIBRARIES)
LIST(REVERSE Trilinos_TPL_LIBRARIES)

#
# some libraries imported as dependencies from Trilinos
# may exist in other locations for this build.
# scan the list for libraries given with an absolute path
# and if found try to find them again (without the absolute path)
#
# foreach( TPL ${Trilinos_TPL_LIBRARIES})
#   message( "TPLs BEFORE abs path removal ${TPL}")
# endforeach()
# 
# foreach( TPL ${Trilinos_TPL_LIBRARIES})
#   if(IS_ABSOLUTE ${TPL})
#     #LIST(REMOVE_ITEM Trilinos_TPL_LIBRARIES ${TPL})
#     get_filename_component( LIBNAME ${TPL} NAME_WE)
#     string(REGEX REPLACE "^lib" "" LIBNAME ${LIBNAME} "" )
#     message( "extracted name is ${LIBNAME}" )
#     #string( CONCAT MKLLIBDIR $ENV{MKLROOT} "/lib/intel64")
#     #find_library( LIBFOUND NAMES ${LIBNAME} HINTS ${MKLLIBDIR} ) 
#     find_library( LIBFOUND NAMES ${LIBNAME} HINTS $ENV{MKLROOT} PATH_SUFFIXES "lib" "lib/intel64" ) 
#     if( LIBFOUND STREQUAL "LIBFOUND-NOTFOUND" )
#       message( "Cound not find library ${TPL} ${LIBNAME}" )
#     else( LIBFOUND STREQUAL "LIBFOUND-NOTFOUND" )
#       message( "Found absolute path library and it's real location ${TPL} ${LIBNAME} ${LIBFOUND}" )
#       LIST(REMOVE_ITEM Trilinos_TPL_LIBRARIES ${TPL})
#       LIST(APPEND Trilinos_TPL_LIBRARIES ${LIBNAME})
#       unset( LIBFOUND CACHE )
#     endif( LIBFOUND STREQUAL "LIBFOUND-NOTFOUND" )
#     #if( EXISTS ${TPL} )
#     #  message( "TPL found ${TPL}" )
#     #else( EXISTS ${TPL} )
#     #  message( "TPL NOT found  ${TPL}" )
#     #endif( EXISTS ${TPL} )
#   endif(IS_ABSOLUTE ${TPL})
##  endforeach()
# 
# foreach( TPL ${Trilinos_TPL_LIBRARIES})
#   message( "AFTER abs path removal ${TPL}")
# endforeach()
set( Xyce_TPL_LIBRARIES ${Trilinos_TPL_LIBRARIES} CACHE STRING "TPLs for Xyce")
add_library(trilinos INTERFACE IMPORTED GLOBAL)

# MPI check
message(STATUS "Checking if MPI is enabled in Trilinos")
list(FIND Trilinos_TPL_LIST MPI MPI_Enabled)
if (MPI_Enabled GREATER -1)
     message(STATUS "Checking if MPI is enabled in Trilinos - MPI enabled")
     set(Xyce_PARALLEL_MPI TRUE CACHE BOOL "Build Xyce with MPI enabled")

     # For MPI builds, Isorropia and Zoltan are REQUIRED
     message(STATUS "Looking for Isorropia in Trilinos")
     list(FIND Trilinos_PACKAGE_LIST Isorropia Isorropia_FOUND)
     if ( Isorropia_FOUND GREATER -1)
          set(Xyce_USE_ISORROPIA TRUE)
          message(STATUS "Looking for Isorropia in Trilinos - found")
     else ()
          message(STATUS "Looking for Isorropia in Trilinos - not found")
          message("Isorropia is required for MPI parallel builds.\n"
               "Enable the following in the Trilinos build:\n"
               "  -D Trilinos_ENABLE_Isorropia=ON")
          set(Trilinos_IS_MISSING_FEATURES TRUE)
     endif()

     message(STATUS "Looking for Zoltan in Trilinos")
     list(FIND Trilinos_PACKAGE_LIST Zoltan Zoltan_FOUND)
     if (Zoltan_FOUND LESS 0)
          message(STATUS "Looking for Zoltan in Trilinos - not found")
          message("Zoltan is required for MPI parallel builds.\n"
               "Enable the following in the Trilinos build:\n"
               "  -D Trilinos_ENABLE_Zoltan=ON")
          set(Trilinos_IS_MISSING_FEATURES TRUE)
     endif()
     message(STATUS "Looking for Zoltan in Trilinos - found")
else()
     message(STATUS "Checking if MPI is enabled in Trilinos - MPI not enabled")
     set(Xyce_PARALLEL_MPI FALSE CACHE BOOL "Build Xyce with MPI enabled")
     set(Xyce_USE_ISORROPIA FALSE)
endif()

# Search for required TPL packages

message(STATUS "Looking for BLAS and LAPACK in Trilinos")
list(FIND Trilinos_TPL_LIST BLAS BLAS_IN_Trilinos)
list(FIND Trilinos_TPL_LIST LAPACK LAPACK_IN_Trilinos)
if ((BLAS_IN_Trilinos GREATER -1) AND (LAPACK_IN_Trilinos GREATER -1))
     message(STATUS "Looking for BLAS and LAPACK in Trilinos - found")
else ()
     message(STATUS "Looking for BLAS and LAPACK in Trilinos - not found")
     message("BLAS and LAPACK are not available via Trilinos.\n"
          "Enable the following in the Trilinos build:\n"
          "  -D TPL_ENABLE_BLAS=ON\n"
          "  -D TPL_ENABLE_LAPACK=ON")
     set(Trilinos_IS_MISSING_FEATURES TRUE)
endif ()

# Search for required features

set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES} ${Trilinos_INCLUDE_DIRS} )
check_cxx_symbol_exists(HAVE_AMESOS_KLU Amesos_config.h KLU_IN_Trilinos)
if (NOT KLU_IN_Trilinos)
     message("Trilinos was not built with KLU support in Amesos.\n"
          "Enable the following in the Trilinos build:\n"
          "  -D Amesos_ENABLE_KLU=ON")
     set(Trilinos_IS_MISSING_FEATURES TRUE)
endif()

check_cxx_symbol_exists(HAVE_BTF EpetraExt_config.h Epetra_BTF_IN_Trilinos)
if (NOT Epetra_BTF_IN_Trilinos)
     message("Trilinos was not built with BTF support in EpetraExt.\n"
          "Enable the following in the Trilinos build:\n"
          "  -D EpetraExt_BUILD_BTF=ON")
     set(Trilinos_IS_MISSING_FEATURES TRUE)
endif()

check_cxx_symbol_exists(HAVE_EXPERIMENTAL EpetraExt_config.h Epetra_EXPERIMENTAL_IN_Trilinos)
if (NOT Epetra_EXPERIMENTAL_IN_Trilinos)
     message("Trilinos was not built with EXPERIMENTAL support in EpetraExt.\n"
          "Enable the following in the Trilinos build:\n"
          "  -D EpetraExt_BUILD_EXPERIMENTAL=ON")
     set(Trilinos_IS_MISSING_FEATURES TRUE)
endif()

check_cxx_symbol_exists(HAVE_GRAPH_REORDERINGS EpetraExt_config.h Epetra_GRAPH_REORD_IN_Trilinos)
if (NOT Epetra_GRAPH_REORD_IN_Trilinos)
     message("Trilinos was not built with GRAPH REORDERINGS support in EpetraExt.\n"
          "Enable the following in the Trilinos build:\n"
          "  -D EpetraExt_BUILD_GRAPH_REORDERINGS=ON")
     set(Trilinos_IS_MISSING_FEATURES TRUE)
endif()

# The following flag should also be enabled in the Trilinos build:
#    -D NOX_ENABLE_LOCA=ON
# However, I don't see a way to probe Trilinos' CMake for this.

check_cxx_symbol_exists(HAVE_TEUCHOS_COMPLEX Teuchos_config.h Teuchos_COMPLEX_IN_Trilinos)
if (NOT Teuchos_COMPLEX_IN_Trilinos)
     message("Trilinos was not built with COMPLEX support in Teuchos.\n"
          "Enable the following in the Trilinos build:\n"
          "  -D Teuchos_ENABLE_COMPLEX=ON")
     set(Trilinos_IS_MISSING_FEATURES TRUE)
endif()

# Should we be checking if Sacado has complex as well?
#    check_cxx_symbol_exists(HAVE_SACADO_COMPLEX Sacado_config.h Sacado_COMPLEX_IN_Trilinos)
#    if (NOT Sacado_COMPLEX_IN_Trilinos)
#         message("Trilinos was not built with COMPLEX support in Sacado.\n"
#              "Enable the following in the Trilinos build:\n"
#              "  -D Sacado_ENABLE_COMPLEX=ON")
#         set(Trilinos_IS_MISSING_FEATURES TRUE)
#    endif()

if (Trilinos_IS_MISSING_FEATURES)
     message(FATAL_ERROR "Halting the Xyce configure due to missing features in Trilinos.\n"
          "Rebuild Trilinos with the required features, and try again.")
endif()

# Search for optional Trilinos packages

# THIS IS A GUESS FOR THE PACKAGE NAME FOR ShyLU; IT NEEDS TO BE VERIFIED
message(STATUS "Looking for ShyLU in Trilinos")
list(FIND Trilinos_PACKAGE_LIST ShyLU ShyLU_IN_Trilinos)
if (ShyLU_IN_Trilinos GREATER -1)
     message(STATUS "Looking for ShyLU in Trilinos - found")
     set(Xyce_SHYLU TRUE CACHE BOOL "Enables the ShyLU linear solver")
else ()
     message(STATUS "Looking for ShyLU in Trilinos - not found")
     set(Xyce_SHYLU FALSE CACHE BOOL "Enables the ShyLU linear solver" FORCE)
endif ()

# THIS IS A GUESS FOR THE PACKAGE NAME FOR Basker; IT NEEDS TO BE VERIFIED
message(STATUS "Looking for Basker in Trilinos")
list(FIND Trilinos_PACKAGE_LIST Basker Basker_IN_Trilinos)
if (Basker_IN_Trilinos GREATER -1)
     message(STATUS "Looking for Basker in Trilinos - found")
     set(SHYLUBASKER TRUE CACHE BOOL "Enables the Basker linear solver")
else ()
     message(STATUS "Looking for Basker in Trilinos - not found")
     set(SHYLUBASKER FALSE CACHE BOOL "Enables the Basker linear solver" FORCE)
endif ()

message(STATUS "Looking for Amesos2 in Trilinos")
list(FIND Trilinos_PACKAGE_LIST Amesos2 Amesos2_IN_Trilinos)
if (Amesos2_IN_Trilinos GREATER -1)
     message(STATUS "Looking for Amesos2 in Trilinos - found")
     set(Xyce_AMESOS2 TRUE CACHE BOOL "Enables the Amesos2 linear solver")
else ()
     message(STATUS "Looking for Amesos2 in Trilinos - not found")
     set(Xyce_AMESOS2 FALSE CACHE BOOL "Enables the Amesos2 linear solver" FORCE)
endif ()

message(STATUS "Looking for Stokhos in Trilinos")
list(FIND Trilinos_PACKAGE_LIST Stokhos Stokhos_IN_Trilinos)
if (Stokhos_IN_Trilinos GREATER -1)
     message(STATUS "Looking for Stokhos in Trilinos - found")
     set(Xyce_STOKHOS_ENABLE TRUE CACHE BOOL "Enables the Stokhos linear solver")
else ()
     message(STATUS "Looking for Stokhos in Trilinos - not found")
     set(Xyce_STOKHOS_ENABLE FALSE CACHE BOOL "Enables the Stokhos linear solver" FORCE)
endif ()

message(STATUS "Looking for ROL in Trilinos")
list(FIND Trilinos_PACKAGE_LIST ROL ROL_IN_Trilinos)
if (ROL_IN_Trilinos GREATER -1)
     message(STATUS "Looking for ROL in Trilinos - found")
     set(Xyce_ROL TRUE CACHE BOOL "Enables the ROL linear solver")
else ()
     message(STATUS "Looking for ROL in Trilinos - not found")
     set(Xyce_ROL FALSE CACHE BOOL "Enables the ROL linear solver" FORCE)
endif ()

# Search for optional TPL packages in Trilinos (some of these are just informational)

message(STATUS "Looking for ParMETIS in Trilinos")
list(FIND Trilinos_TPL_LIST ParMETIS ParMETIS_IN_Trilinos)
if (ParMETIS_IN_Trilinos GREATER -1)
     message(STATUS "Looking for ParMETIS in Trilinos - found")
     set(Xyce_USE_PARMETIS TRUE CACHE BOOL "Enables the ParMETIS linear solver")
else()
     message(STATUS "Looking for ParMETIS in Trilinos - not found")
endif()

message(STATUS "Looking for AMD in Trilinos")
list(FIND Trilinos_TPL_LIST AMD AMD_IN_Trilinos)
if (AMD_IN_Trilinos GREATER -1)
     message(STATUS "Looking for AMD in Trilinos - found")
     set(Xyce_AMD TRUE CACHE BOOL "Enables the option of AMD ordering for the linear solver")
else()
     message(STATUS "Looking for AMD in Trilinos - not found")
endif()

message(STATUS "Looking for PARDISO_MKL in Trilinos")
list(FIND Trilinos_TPL_LIST PARDISO_MKL PARDISO_IN_Trilinos)
if (PARDISO_IN_Trilinos GREATER -1)
     message(STATUS "Looking for PARDISO_MKL in Trilinos - found")
else()
     message(STATUS "Looking for PARDISO_MKL in Trilinos - not found")
endif()

message(STATUS "Looking for SuperLU in Trilinos")
list(FIND Trilinos_TPL_LIST SuperLU SuperLU_IN_Trilinos)
if (SuperLU_IN_Trilinos GREATER -1)
     message(STATUS "Looking for SuperLU in Trilinos - found")
else()
     message(STATUS "Looking for SuperLU in Trilinos - not found")
endif()

message(STATUS "Looking for SuperLUDist in Trilinos")
list(FIND Trilinos_TPL_LIST SuperLUDist SuperLUDist_IN_Trilinos)
if (SuperLUDist_IN_Trilinos GREATER -1)
     message(STATUS "Looking for SuperLUDist in Trilinos - found")
else()
     message(STATUS "Looking for SuperLUDist in Trilinos - not found")
endif()

###################
## End Trilinos
###################


# Find a usable FFT library
message(STATUS "Looking for usable FFT libraries")
string( FIND "${Trilinos_TPL_LIBRARIES}" "mkl_core" MKL_LIBS_FOUND )
if (MKL_LIBS_FOUND GREATER -1)

     message(STATUS "Looking for usable FFT libraries - found the Intel Math Kernel Library")
     set(Xyce_USE_FFT TRUE CACHE BOOL "Enable the FFT capability")
     set(Xyce_USE_INTEL_FFT TRUE CACHE BOOL "Use the Intel Math Kernel Library FFT capability")
     
     set(FFT "")
     #So far when using Intel compilers we must line up with a Trilinos built with MKL,
     #and that then must be exposed as a Trilinos TPL, so specifing it in the link line
     #is unnecessary given that we already include Trilinos in the link line.
     #Because the IntelMKL is technically free to redistribute, I'm leaving this here, although
     #how much use the free Intel MKL is without the corresponding notsofree compilers is something
     #the community could help us with. 
else ()
     find_package(FFTW)
     if(FFTW_FOUND)
          set(Xyce_USE_FFT TRUE CACHE BOOL "Enable the FFT capability")
          add_library(FFTW::FFTW INTERFACE IMPORTED GLOBAL)
          set_target_properties(FFTW::FFTW PROPERTIES
               INTERFACE_INCLUDE_DIRECTORIES "${FFTW_INCLUDE_DIRS}"
               INTERFACE_LINK_LIBRARIES "${FFTW_DOUBLE_LIB}")
          set(Xyce_USE_FFTW TRUE CACHE BOOL "Use FFTW library")
	  set(FFT FFTW::FFTW) 
     else()
          message("Neither FFTW or Intel MKL found - disabling the FFT capability")
          set(Xyce_USE_FFT FALSE CACHE BOOL "Enable the FFT capability")
     endif ()
endif ()

# Find flex and Bison
message(STATUS "Looking for flex and Bison")
find_package(FLEX)
# The autotools probe also does the following:
#    "Define YYTEXT_POINTER if yytext defaults to 'char *' instead of 'char'"
# We do not appear to use YYTEXT_POINTER anywhere, so is it safe to not probe
# for that behavior?
find_package(BISON 2.4)
# The 2.4 specifies the minimum version.  That is ok at the moment, as Bison
# has been functional for many versions (through 3.6 at the time of this
# writing).  Historically, though, new versions have had backward
# incompatibility issues.  If that occurs again, the BISON_VERSION variable
# will have to be probed for a certain range.
if (FLEX_FOUND AND BISON_FOUND)
     set(Xyce_REACTION_PARSER TRUE CACHE BOOL "Enable the chemical reaction parsing capability")
else()
     set(Xyce_REACTION_PARSER FALSE CACHE BOOL "Enable the chemical reaction parsing capability")
endif ()

# Find CURL
if (Xyce_USE_CURL)
     if (Xyce_TRACKING_URL)
       message(STATUS "Looking for cURL")
       # We have been having trouble on Windows, where it seems that
       # CURLConfig.cmake is not getting the right variables set.
       #  Let's try telling FindCURL not to use it.
          set(CURL_NO_CURL_CMAKE ON)
          find_package(CURL REQUIRED)
          # Now, what did this tell us?
          message(STATUS "The usage tracking capability is enabled. Using: ${Xyce_TRACKING_URL}")
          message(STATUS "CURL_LIBRARIES is ${CURL_LIBRARIES} and CURL_INCLUDE_DIRS is ${CURL_INCLUDE_DIRS}")
     else()
          message("Xyce_USE_CURL is TRUE, but no URL is supplied in Xyce_TRACKING_URL.\n"
                  "Changing Xyce_USE_CURL to FALSE - disabling usage tracking")
          set(Xyce_USE_CURL FALSE CACHE BOOL "Enable the usage tracking capability using CURL" FORCE)
     endif()
else()
     message(STATUS "Usage tracking is not enabled")
endif()

