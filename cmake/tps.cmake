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
#    Trilinos/cmake/doc/export_system/Finding_Trilinos.txt
# Also see the directories:
#    Trilinos/commonTools/importing/
#    Trilinos/demos/buildAgainstTrilinos/
# Some of that information is out-of-date, though.

# The following is recommended to make the library lists easier to read, as
# they contain a lot of duplicates.  However, this is done automatically in
# recent versions of Trilinos, so we will comment it out here.
#    LIST(REVERSE Trilinos_LIBRARIES)
#    LIST(REMOVE_DUPLICATES Trilinos_LIBRARIES)
#    LIST(REVERSE Trilinos_LIBRARIES)

# It is not safe to do the above for the TPLs, because someone might specify
# the libraries explicitly using flags like -L, -l, etc. (despite that being
# "bad form" for CMake builds).  In that case, order and duplication within
# that order matter.  Nevertheless, cleaning up the lists is sometimes useful
# for troubleshooting, so we will leave this here for that reason.
#    LIST(REVERSE Trilinos_TPL_LIBRARIES)
#    LIST(REMOVE_DUPLICATES Trilinos_TPL_LIBRARIES)
#    LIST(REVERSE Trilinos_TPL_LIBRARIES)

message(STATUS "Looking for Trilinos\n"
  "   Required packages:\n"
  "        Amesos Epetra EpetraExt Ifpack NOX Teuchos Sacado\n"
  "        Triutils AztecOO Belos TrilinosCouplings\n"
  "   Optional packages:\n"
  "        Isorropia Zoltan ShyLU ShyLU_DDCore Amesos2 Stokhos ROL MKL")
find_package(Trilinos CONFIG
  REQUIRED Amesos Epetra EpetraExt Ifpack NOX Teuchos Sacado Triutils
  AztecOO Belos TrilinosCouplings
  OPTIONAL_COMPONENTS Isorropia Zoltan ShyLU ShyLU_DDCore
  Amesos2 Stokhos ROL MKL)
message(STATUS "Looking for Trilinos - found")

if(Trilinos_VERSION VERSION_LESS "13.5")
  message(FATAL_ERROR
    "ERROR: Trilinos version ${Trilinos_VERSION} is less than the required minimum of 13.5. Install a version of Trilinos of 13.5 or greater.")
endif()

# MPI check
message(STATUS "Checking if MPI is enabled in Trilinos")
if(TARGET MPI::all_libs)
     message(STATUS "Checking if MPI is enabled in Trilinos - MPI enabled")
     set(Xyce_PARALLEL_MPI TRUE CACHE BOOL "Build Xyce with MPI enabled" FORCE)

     # For MPI builds, Isorropia and Zoltan are REQUIRED
     message(STATUS "Looking for Isorropia in Trilinos")
     if(TARGET Isorropia::all_libs)
          set(Xyce_USE_ISORROPIA TRUE)
          message(STATUS "Looking for Isorropia in Trilinos - found")
     else()
          message(STATUS "Looking for Isorropia in Trilinos - not found")
          message("Isorropia is required for MPI parallel builds.\n"
               "Enable the following in the Trilinos build:\n"
               "  -D Trilinos_ENABLE_Isorropia=ON")
          set(Trilinos_IS_MISSING_FEATURES TRUE)
     endif()

     message(STATUS "Looking for Zoltan in Trilinos")
     if(NOT TARGET Zoltan::all_libs)
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

message(STATUS "Looking for BLAS via Trilinos")
if(TARGET BLAS::all_libs)
     message(STATUS "Looking for BLAS via Trilinos - found")
else()
     message(STATUS "Looking for BLAS via Trilinos - not found")
     message("BLAS is not available via Trilinos.\n"
          "Enable the following in the Trilinos build:\n"
          "  -D TPL_ENABLE_BLAS=ON")
     set(Trilinos_IS_MISSING_FEATURES TRUE)
endif()

message(STATUS "Looking for LAPACK via Trilinos")
if(TARGET LAPACK::all_libs)
     message(STATUS "Looking for LAPACK via Trilinos - found")
else()
     message(STATUS "Looking for LAPACK via Trilinos - not found")
     message("LAPACK is not available via Trilinos.\n"
          "Enable the following in the Trilinos build:\n"
          "  -D TPL_ENABLE_LAPACK=ON")
     set(Trilinos_IS_MISSING_FEATURES TRUE)
endif()

# Search for required features

# Since Trilinos depends on OpenMP, we look for it first, so it is available
# for some of the Trilinos feature probes.
list(FIND Kokkos_DEVICES OPENMP OpenMP_IN_Kokkos)
if (OpenMP_IN_Kokkos GREATER -1)
   find_package(OpenMP REQUIRED)
endif()

if(NOT TARGET NOX::loca)
     message("Trilinos was not built with LOCA support in NOX.\n"
          "Enable the following in the Trilinos build:\n"
          "  -D NOX_ENABLE_LOCA=ON")
     set(Trilinos_IS_MISSING_FEATURES TRUE)
endif()

get_target_property(CMAKE_REQUIRED_INCLUDES Amesos::all_libs INTERFACE_INCLUDE_DIRECTORIES)
check_cxx_symbol_exists(HAVE_AMESOS_KLU Amesos_config.h KLU_IN_Trilinos)
if (NOT KLU_IN_Trilinos)
     message("Trilinos was not built with KLU support in Amesos.\n"
          "Enable the following in the Trilinos build:\n"
          "  -D Amesos_ENABLE_KLU=ON")
     set(Trilinos_IS_MISSING_FEATURES TRUE)
endif()

get_target_property(CMAKE_REQUIRED_INCLUDES EpetraExt::all_libs INTERFACE_INCLUDE_DIRECTORIES)
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

get_target_property(CMAKE_REQUIRED_INCLUDES Teuchos::all_libs INTERFACE_INCLUDE_DIRECTORIES)
check_cxx_symbol_exists(HAVE_TEUCHOS_COMPLEX Teuchos_config.h Teuchos_COMPLEX_IN_Trilinos)
if (NOT Teuchos_COMPLEX_IN_Trilinos)
     message("Trilinos was not built with COMPLEX support in Teuchos.\n"
          "Enable the following in the Trilinos build:\n"
          "  -D Trilinos_ENABLE_COMPLEX_DOUBLE=ON")
     set(Trilinos_IS_MISSING_FEATURES TRUE)
endif()

# Compilers targeting the MSVC ABI have flags to select the MSVC runtime
# library.  This is because, with the MSVC, there is a separate runtime library
# for Debug builds.  (This is related to the CMAKE_BUILD_TYPE discussion in the
# root "CMakeLists.txt file.  When the minimum CMake version requirement is
# increased to 3.15 or above, see CMake policy CMP0091:
#    <https://cmake.org/cmake/help/latest/policy/CMP0091.html>)
# Teuchos links to the MSVC runtime library.  Since LOCA.H includes at least one
# Teuchos header, the test below has to point to the same library as Trilinos.
# (Note that the debug library is the default for the test.)  The best way to
# accomplish this is to use the same compiler flags as Trilinos.
# One might think to use:
#   set (CMAKE_REQUIRED_FLAGS ${Trilinos_CXX_COMPILER_FLAGS})
# or
#   check_include_file_cxx(LOCA.H LOCA_IN_Trilinos ${Trilinos_CXX_COMPILER_FLAGS})
# However, the test *prepends* those flags, so they get overridden by the
# debug flags of the check, itself.  The following is probably an abuse of
# CMAKE_REQUIRED_DEFINITIONS, but it seems to work.

if (MSVC)
     set(CMAKE_REQUIRED_DEFINITIONS "${Trilinos_CXX_COMPILER_FLAGS}")
endif()
get_target_property(CMAKE_REQUIRED_LIBRARIES Teuchos::all_libs INTERFACE_LINK_LIBRARIES)

# Perform an initial check to see if we can compile against Trilinos at all.
# This could reveal compiler setup problems and/or Trilinos setup problems.
check_include_file_cxx(Teuchos_SerialDenseMatrix.hpp Trilinos_COMPILE_SUCCESS ${OpenMP_CXX_FLAGS})
if (NOT Trilinos_COMPILE_SUCCESS)
     message(FATAL_ERROR "Unable to compile against Trilinos. It is possible\
     Trilinos was not properly configured, or the environment has changed since\
     Trilinos was installed. See the CMake log files for more information.")
endif()

# After the release of Trilinos 12.12.1, the abstract solver interface in NOX
# was changed to include a new method that returns solver statistics.  This
# test and the Xyce_NOX_SOLVERSTATS ifdefs can be removed if the minimum
# required version of Trilinos is raised.
# 8/24/2022 - Minimum version was raised to 13.5 for cmake, but not autotools.
#             This ifdef will be set to true for cmake builds.
get_target_property(CMAKE_REQUIRED_INCLUDES NOX::all_libs INTERFACE_INCLUDE_DIRECTORIES)
get_target_property(CMAKE_REQUIRED_LIBRARIES NOX::all_libs INTERFACE_LINK_LIBRARIES)
set(Xyce_NOX_SOLVERSTATS TRUE CACHE BOOL "Use new method to return NOX solver statistics.")

unset(CMAKE_REQUIRED_INCLUDES)
unset(CMAKE_REQUIRED_LIBRARIES)
if (MSVC)
     unset(CMAKE_REQUIRED_DEFINITIONS)
endif()

if (Trilinos_IS_MISSING_FEATURES)
     message(FATAL_ERROR "Halting the Xyce configure due to missing features in Trilinos.\n"
          "Rebuild Trilinos with the required features, and try again.")
endif()

# Search for optional Trilinos packages

# Hybrid-hybrid ShyLU linear solver 
if (DEFINED Xyce_SHYLU AND NOT Xyce_SHYLU)
   set(Xyce_SHYLU FALSE CACHE BOOL "Enables the ShyLU linear solver package")
else()
   message(STATUS "Looking for ShyLU in Trilinos")
   if(TARGET ShyLU_DDCore::all_libs)
        message(STATUS "Looking for ShyLU in Trilinos - found")
        set(Xyce_SHYLU TRUE CACHE BOOL "Enables the ShyLU linear solver package")
   else()
        message(STATUS "Looking for ShyLU in Trilinos - not found")
        set(Xyce_SHYLU FALSE CACHE BOOL "Enables the ShyLU linear solver package" FORCE)
   endif()
endif()

if (DEFINED Xyce_AMESOS2 AND NOT Xyce_AMESOS2)
     set(Xyce_AMESOS2 FALSE CACHE BOOL "Enables the Amesos2 linear solver")
else()
     message(STATUS "Looking for Amesos2 in Trilinos")
     if(TARGET Amesos2::all_libs)
          set(Xyce_AMESOS2 TRUE CACHE BOOL "Enables the Amesos2 linear solver")
          message(STATUS "Looking for Amesos2 in Trilinos - found")
     else()
          message(STATUS "Looking for Amesos2 in Trilinos - not found")
          set(Xyce_AMESOS2 FALSE CACHE BOOL "Enables the Amesos2 linear solver" FORCE)
     endif()
endif()

if (Xyce_AMESOS2)
     get_target_property(CMAKE_REQUIRED_INCLUDES Amesos2::all_libs INTERFACE_INCLUDE_DIRECTORIES)

     # KLU2
     if (DEFINED Xyce_AMESOS2_KLU2 AND NOT Xyce_AMESOS2_KLU2)
          set(Xyce_AMESOS2_KLU2 FALSE CACHE BOOL "Enables the templated KLU2 linear solver in Amesos2")
     else()
          check_cxx_symbol_exists(HAVE_AMESOS2_KLU2 Amesos2_config.h Amesos2_KLU2_IN_Trilinos)
          if (Amesos2_KLU2_IN_Trilinos)
               set(Xyce_AMESOS2_KLU2 TRUE CACHE BOOL "Enables the templated KLU2 linear solver in Amesos2")
          else()
               set(Xyce_AMESOS2_KLU2 FALSE CACHE BOOL "Enables the templated KLU2 linear solver in Amesos2" FORCE)
          endif()
     endif()

     # Templated Basker
     if (DEFINED Xyce_AMESOS2_BASKER AND NOT Xyce_AMESOS2_BASKER)
          set(Xyce_AMESOS2_BASKER FALSE CACHE BOOL "Enables the templated Basker linear solver in Amesos2")
     else()
          check_cxx_symbol_exists(HAVE_AMESOS2_BASKER Amesos2_config.h Amesos2_Basker_IN_Trilinos)
          if (Amesos2_Basker_IN_Trilinos)
               set(Xyce_AMESOS2_BASKER TRUE CACHE BOOL "Enables the templated Basker linear solver in Amesos2")
               set(Xyce_NEW_BASKER TRUE)
               # After the release of Trilinos 12.12 (maybe 12.14?), the Amesos2/Basker interface was changed.
               # Xyce's cmake builds only support the new basker interface, but the autotools build still
               # supports options for using both. The Xyce_NEW_BASKER ifdefs can be removed if the minimum 
               # required version of Trilinos is raised for both build systems.
          else()
               set(Xyce_AMESOS2_BASKER FALSE CACHE BOOL "Enables the templated Basker linear solver in Amesos2" FORCE)
          endif()
     endif()

     # ShyLU-Basker
     if (DEFINED Xyce_AMESOS2_SHYLUBASKER AND NOT Xyce_AMESOS2_SHYLUBASKER)
          set(Xyce_AMESOS2_SHYLUBASKER FALSE CACHE BOOL "Enables the multi-threaded ShyLU-Basker linear solver in Amesos2")
     else()
          check_cxx_symbol_exists(HAVE_AMESOS2_SHYLUBASKER Amesos2_config.h Amesos2_ShyLUBasker_IN_Trilinos)
          check_cxx_symbol_exists(HAVE_AMESOS2_SHYLU_NODEBASKER Amesos2_config.h Amesos2_ShyLUNodeBasker_IN_Trilinos)
          if (Amesos2_ShyLUBasker_IN_Trilinos OR Amesos2_ShyLUNodeBasker_IN_Trilinos)
               set(Xyce_AMESOS2_SHYLUBASKER TRUE CACHE BOOL "Enables the multi-threaded ShyLU-Basker linear solver in Amesos2")
          else()
               set(Xyce_AMESOS2_SHYLUBASKER FALSE CACHE BOOL "Enables the multi-threaded ShyLU-Basker linear solver in Amesos2" FORCE)
          endif()
     endif()

     unset(CMAKE_REQUIRED_FLAGS)
     unset(CMAKE_REQUIRED_LIBRARIES)
     unset(CMAKE_REQUIRED_INCLUDES)
endif()

if (DEFINED Xyce_STOKHOS_ENABLE AND NOT Xyce_STOKHOS_ENABLE)
     set(Xyce_STOKHOS_ENABLE FALSE CACHE BOOL "Enables the Stokhos linear solver")
else()
     message(STATUS "Looking for Stokhos in Trilinos")
     if(TARGET Stokhos::all_libs)
          message(STATUS "Looking for Stokhos in Trilinos - found")
          set(Xyce_STOKHOS_ENABLE TRUE CACHE BOOL "Enables the Stokhos linear solver")
     else()
          message(STATUS "Looking for Stokhos in Trilinos - not found")
          set(Xyce_STOKHOS_ENABLE FALSE CACHE BOOL "Enables the Stokhos linear solver" FORCE)
     endif()
endif()

if (DEFINED Xyce_ROL AND NOT Xyce_ROL)
     set(Xyce_ROL FALSE CACHE BOOL "Enables the ROL linear solver")
else()
     message(STATUS "Looking for ROL in Trilinos")
     if(TARGET ROL::all_libs)
          set(Xyce_ROL TRUE CACHE BOOL "Enables the ROL linear solver")
          message(STATUS "Looking for ROL in Trilinos - found")
     else()
          message(STATUS "Looking for ROL in Trilinos - not found")
          set(Xyce_ROL FALSE CACHE BOOL "Enables the ROL linear solver" FORCE)
     endif()
endif()

# Search for optional TPL packages in Trilinos

# Because of the way the autotools script works, it is not trivial to change
# HAVE_LIBPARMETIS to something like Xyce_USE_PARMETIS.  We will simply use the
# autotools variable for now (which comes from Trilinos), and change this some
# time in the future to help with consistency.
if (DEFINED HAVE_LIBPARMETIS AND NOT HAVE_LIBPARMETIS)
     set(HAVE_LIBPARMETIS FALSE CACHE BOOL "Enables the ParMETIS partitioning library")
else()
     message(STATUS "Looking for ParMETIS via Trilinos")
     if(TARGET ParMETIS::all_libs)
          message(STATUS "Looking for ParMETIS via Trilinos - found")
          set(HAVE_LIBPARMETIS TRUE CACHE BOOL "Enables the ParMETIS partitioning library")
     else()
          message(STATUS "Looking for ParMETIS via Trilinos - not found")
          set(HAVE_LIBPARMETIS FALSE CACHE BOOL "Enables the ParMETIS partitioning library" FORCE)
     endif()
endif()

if (DEFINED Xyce_AMD AND NOT Xyce_AMD)
     set(Xyce_AMD FALSE CACHE BOOL "Enables the option of AMD ordering for the linear solver")
else()
     message(STATUS "Looking for AMD via Trilinos")
     if(TARGET AMD::all_libs)
          message(STATUS "Looking for AMD via Trilinos - found")
          set(Xyce_AMD TRUE CACHE BOOL "Enables the option of AMD ordering for the linear solver")

          # libamd may not always depend on libsuitesparseconfig.  So this extra 
          # check is not safe.
          # libamd is dependent on libsuitesparseconfig. this should
          # be handled wherever AMD::all_libs is created but it
          # isn't. when building shared libraries it all works out but
          # when building static libraries you have to explicitly set
          # this dependency.
          #find_library(SUITESPARSECONFIG_LIB NAMES suitesparseconfig
          #  HINTS $ENV{Trilinos_DIR}/lib64 $ENV{Trilinos_DIR}/lib)
          #
          #if (SUITESPARSECONFIG_LIB)
          #  target_link_libraries(AMD::all_libs INTERFACE ${SUITESPARSECONFIG_LIB})
          #else()
          #  message(WARNING "Unable to find libsuitesparseconfig. Note that this may cause unresolved SuiteSParse*() functions during link when trilinos/amd/suitesparse are built statically")
          #endif()
     else()
          message(STATUS "Looking for AMD via Trilinos - not found")
          set(Xyce_AMD FALSE CACHE BOOL "Enables the option of AMD ordering for the linear solver" FORCE)
     endif()
endif()

# The following are simply for reporting purposes; no build variables need to
# be defined.

message(STATUS "Looking for PARDISO_MKL via Trilinos")
if (TARGET PARDISO_MKL::all_libs)
     message(STATUS "Looking for PARDISO_MKL via Trilinos - found")
else()
     message(STATUS "Looking for PARDISO_MKL via Trilinos - not found")
endif()

message(STATUS "Looking for SuperLU via Trilinos")
if (TARGET SuperLU::all_libs)
     message(STATUS "Looking for SuperLU via Trilinos - found")
else()
     message(STATUS "Looking for SuperLU via Trilinos - not found")
endif()

message(STATUS "Looking for SuperLUDist via Trilinos")
if (TARGET SuperLUDist::all_libs)
     message(STATUS "Looking for SuperLUDist via Trilinos - found")
else()
     message(STATUS "Looking for SuperLUDist via Trilinos - not found")
endif()

###################
## End Trilinos
###################

###################################
## Find a usable FFT library
###################################

# Address the pathological case first
if (Xyce_USE_FFT AND Xyce_USE_INTEL_FFT AND Xyce_USE_FFTW)
     message(FATAL_ERROR
          "Both \"Xyce_USE_INTEL_FFT\" and \"Xyce_USE_FFTW\" have been set as\
          TRUE.  It is recommended to delete all \"FFT\" variables and let\
          CMake determine the configuration.  You may also set either\
          \"Xyce_USE_INTEL_FFT\" or \"Xyce_USE_FFTW\" to TRUE.  However, if you\
          explicitly specify the use of the Intel MKL FFT capability, you must\
          also ensure the appropriate MKL flags are set in the Xyce CMake\
          invocation.")
endif()

# Both Xyce_USE_INTEL_FFT and Xyce_USE_FFT must be true to force the use of the
# Intel MKL (see below).  If Xyce_USE_INTEL_FFT is true and Xyce_USE_FFT is
# unspecified or explicitly false, then that is considered a contradiction.
if (Xyce_USE_INTEL_FFT AND NOT Xyce_USE_FFT)
     set(Xyce_USE_INTEL_FFT FALSE CACHE BOOL "Use the Intel Math Kernel Library FFT capability" FORCE)
endif()

# Note that the following isn't forced, so it does not override an explicit set
# of Xyce_USE_FFT=False.
set(Xyce_USE_FFT TRUE CACHE BOOL "Enable the FFT capability")

if (Xyce_USE_FFT AND Xyce_USE_INTEL_FFT)
     message("\"Xyce_USE_FFT\" and \"Xyce_USE_INTEL_FFT\" are TRUE.  It is assumed "
             "the Intel MKL was found in a previous configuration run, or the user "
             "is explicitly enabling it.")
elseif (Xyce_USE_FFT)
     message(STATUS "Looking for FFT libraries")
endif()

# If Xyce_USE_FFTW is true and Xyce_USE_FFT is unspecified or explicitly true,
# then it is taken that the user wants to use FFTW.  If Xyce_USE_FFTW is true and
# Xyce_USE_FFT is explicitly false, then that is considered a contradiction.
# Note that this is slightly different from the Intel MKL case, so this must be
# located after the "set(Xyce_USE_FFT..." command, above.
if (Xyce_USE_FFTW AND NOT Xyce_USE_FFT)
     set(Xyce_USE_FFTW FALSE CACHE BOOL "Use FFTW library" FORCE)
endif()

# Note that the following is not done on a re-run of the configuration, when
# Xyce_USE_INTEL_FFT=TRUE (as might happen with ccmake).  That should not be a
# problem, since we're relying on information from Trilinos, anyway.
if (Xyce_USE_FFT AND NOT Xyce_USE_INTEL_FFT AND NOT Xyce_USE_FFTW)
     if(TARGET MKL::all_libs)
          set(Tri_KNOWS_MKL TRUE)
          get_target_property(CMAKE_REQUIRED_INCLUDES MKL::all_libs INTERFACE_INCLUDE_DIRECTORIES)
          check_include_file_cxx("mkl_dfti.h" HAVE_MKL_FFT)
          unset(CMAKE_REQUIRED_INCLUDES)
     endif()
     if (Tri_KNOWS_MKL AND HAVE_MKL_FFT)
          message(STATUS "Looking for FFT libraries - found the Intel Math Kernel Library")
          set(Xyce_USE_INTEL_FFT TRUE CACHE BOOL "Use the Intel Math Kernel Library FFT capability" FORCE)
          # Since the MKL is included in the "Trilinos_TPL_LIBRARIES" list, or CMake
          # is already aware of it, specifing it on the link line again is unnecessary.
          set(FFT "")
     elseif (HAVE_MKL_FFT AND NOT Tri_KNOWS_MKL)
          message(STATUS "Looking for FFT libraries - found the Intel Math Kernel Library")
          message(WARNING "The Intel Math Kernel Library was found, but it is not known to Trilinos.\n"
                  "Disabling the Intel Math Kernel Library.\n"
                  "If you want to force the use of the Intel MKL FFT capability, set "
                  "\"Xyce_USE_FFT=TRUE\" and \"Xyce_USE_INTEL_FFT=TRUE\", and set "
                  "the appropriate MKL flags in the Xyce CMake invocation (see "
                  "https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor).")
     elseif (Tri_KNOWS_MKL AND NOT HAVE_MKL_FFT)
          message(STATUS "Looking for FFT libraries - Intel Math Kernel Library not found")
          message(WARNING "Trilinos seems to be aware of the Intel Math Kernel Library, but the FFT "
                  "capability was not found.  This may be an indication of a problem with "
                  "the current platform configuration.")
     else()
          message(STATUS "Looking for FFT libraries - Intel Math Kernel Library not found")
     endif()
endif()

# If the Intel MKL is not being used, try to find FFTW.
if (Xyce_USE_FFT AND NOT Xyce_USE_INTEL_FFT)
     if( NOT (Xyce_USE_APPLEFFT OR Xyce_USE_BASICFFT))
          find_package(FFTW)
          if(FFTW_FOUND)
               message(STATUS "Looking for FFT libraries - found FFTW")
               add_library(FFTW::FFTW INTERFACE IMPORTED GLOBAL)
               set_target_properties(FFTW::FFTW PROPERTIES
                    INTERFACE_INCLUDE_DIRECTORIES "${FFTW_INCLUDE_DIRS}"
                    INTERFACE_LINK_LIBRARIES "${FFTW_DOUBLE_LIB}")
               set(Xyce_USE_FFTW TRUE CACHE BOOL "Use FFTW library")
               set(FFT FFTW::FFTW)
          endif()
     endif()
     if( APPLE AND Xyce_USE_APPLEFFT)
          # try searching for Apple's FFT functions 
          find_library(Accelerate_Framework Accelerate)
          if(Accelerate_Framework)
               message(STATUS "Looking for Apple FFT libraries - found framework Accelerate")
               # note ${Accelerate_Framework} will point to Apple's library with full path.
               # save this in the Cache for later use by downstream products 
               #message(STATUS "${Accelerate_Framework}")
               #add_library(${Accelerate_Framework} INTERFACE IMPORTED GLOBAL)
               set(Xyce_AppleAccelerate ${Accelerate_Framework} CACHE STRING "Apple Accelerate library")
               set(Xyce_USE_APPLE_FFT TRUE CACHE BOOL "Use Apple FFT Library")
          endif()
     else()
       # fall back to the basic FFT implementation 
       message(STATUS "Using basic FFT class for FFT Support.")
       set(Xyce_USE_BASIC_FFT TRUE CACHE BOOL "Use Basic FFT Library")
     endif()
endif()
if (Xyce_USE_FFT AND NOT Xyce_USE_INTEL_FFT AND NOT Xyce_USE_FFTW AND NOT Xyce_USE_APPLEFFT AND NOT Xyce_USE_BASICFFT)
     message("Neither FFTW, Apple FFT or Intel MKL enabled - disabling the FFT capability")
     set(Xyce_USE_FFT FALSE CACHE BOOL "Enable the FFT capability" FORCE)
endif()

###################################
## End find a usable FFT library
###################################

# Find flex and Bison
message(STATUS "Looking for flex and Bison")
find_package(FLEX 2.6 REQUIRED)
find_package(BISON 3.3 REQUIRED)

# Find CURL
if (Xyce_USE_CURL)
     if (Xyce_TRACKING_URL)
          message(STATUS "Looking for cURL")
          find_package(CURL REQUIRED)
          message(STATUS "The usage tracking capability is enabled.  Using: ${Xyce_TRACKING_URL}")
     else()
          message("Xyce_USE_CURL is TRUE, but no URL is supplied in Xyce_TRACKING_URL.\n"
                  "Changing Xyce_USE_CURL to FALSE - disabling usage tracking")
          set(Xyce_USE_CURL FALSE CACHE BOOL "Enable the usage tracking capability using CURL" FORCE)
     endif()
else()
     message(STATUS "Usage tracking is not enabled")
endif()

find_package(Git)

#
# Look for optional Dakota libraries to build integrated Xyce-Dakota binary
#
if (Xyce_DAKOTA OR Xyce_Dakota)
     find_package(Dakota REQUIRED)
     if(Dakota_FOUND)
          message(STATUS "Looking for Dakota libraries - Dakota found")

          message(STATUS "Looking for Boost")
          find_package(Boost 1.69 REQUIRED
            COMPONENTS filesystem program_options regex serialization system)

          set(Xyce_Dakota TRUE CACHE BOOL "Use Dakota library")
     endif()

endif()

include(CTest)

if(Xyce_GTEST_UNIT_TESTS)
     message(STATUS "Looking for GTest")
     # If the wrong version of GTest is found, try setting GTest_DIR or GTEST_ROOT to the install
     # directory of the desired verstion of GTest; specify when invoking cmake to configure.
     find_package(GTest)
     include(GoogleTest)
endif()
