# This should always be true, so why is it optional?
set(Xyce_SPICE_NORMS TRUE)

# For communicating the exact version of Xyce to the code
set(Xyce_RELEASE                   TRUE CACHE BOOL "Set to TRUE to designate a release version")
set(Xyce_QUALIFICATION             TRUE CACHE BOOL "Set to TRUE to designate a qualification release")
if(Xyce_RELEASE)
     set(RELEASE_CHARACTER "R")
     if(Xyce_QUALIFICATION)
          set(QUALIFICATION_CHARACTER "D")
     endif()
else()
     set(RELEASE_CHARACTER "D")
endif()
set(Xyce_XyceVERSION "${RELEASE_CHARACTER}:${QUALIFICATION_CHARACTER}:${Xyce_VERSION_MAJOR}.${Xyce_VERSION_MINOR}.${Xyce_VERSION_PATCH}")

# generate the version strings
set ( Xyce_VERSION_STRING_SHORT "${Xyce_VERSION_MAJOR}.${Xyce_VERSION_MINOR}" )

if ( DEFINED Xyce_VERSION_PATCH AND NOT (Xyce_VERSION_PATCH STREQUAL "" OR Xyce_VERSION_PATCH STREQUAL "0"))
  set ( Xyce_VERSION_STRING_SHORT "${Xyce_VERSION_STRING_SHORT}.${Xyce_VERSION_PATCH}" )
endif ( DEFINED Xyce_VERSION_PATCH AND NOT (Xyce_VERSION_PATCH STREQUAL "" OR Xyce_VERSION_PATCH STREQUAL "0"))

set ( Xyce_VERSION_STRING_LONG "${Xyce_VERSION_STRING_SHORT}" )

if ( DEFINED Xyce_VERSION_EXTRA AND NOT(Xyce_VERSION_EXTRA STREQUAL ""))
  set ( Xyce_VERSION_STRING_LONG "${Xyce_VERSION_STRING_LONG} ${Xyce_VERSION_EXTRA}" )
endif ( DEFINED Xyce_VERSION_EXTRA AND NOT(Xyce_VERSION_EXTRA STREQUAL ""))

if ( DEFINED QUALIFICATION_CHARACTER AND NOT(QUALIFICATION_CHARACTER STREQUAL ""))
    set ( Xyce_VERSION_STRING_LONG "${Xyce_VERSION_STRING_LONG} ${QUALIFICATION_CHARACTER}" )
endif ( DEFINED QUALIFICATION_CHARACTER AND NOT(QUALIFICATION_CHARACTER STREQUAL ""))

# Enable the Unit Tests
option(Xyce_TEST_SUITE "Enables the unit tests" OFF)

# Enable the Plugin capability
option(Xyce_PLUGIN_SUPPORT "Install Xyce with plugin compatibility" OFF)

# Enable parallel
set(Xyce_PARALLEL_MPI              FALSE CACHE BOOL "Build Xyce with MPI enabled")

# Tracking capabilities
set(Xyce_USE_CURL                  FALSE CACHE BOOL "Enable the usage tracking capability using CURL")
set(Xyce_TRACKING_URL              ""  CACHE STRING "The URL for the usage tracking capability")

# Speed up compiles by disabling the analytic sensitivities in the
# ADMS-generated devices
set(Xyce_ADMS_SENSITIVITIES        TRUE CACHE BOOL "Enable analytic sensitivities in ADMS-generated devices")

# Self-explanatory: Enable the chemical reaction parsing capability
set(Xyce_REACTION_PARSER           TRUE CACHE BOOL "Enable the chemical reaction parsing capability")

# Support for Charon coupling
set(Xyce_CHARON                    FALSE CACHE BOOL "Enable Charon device support")

# Verbose output
set(Xyce_VERBOSE_CONDUCTANCE       FALSE CACHE BOOL "Enable verbose output for ???")
set(Xyce_VERBOSE_LINEAR            FALSE CACHE BOOL "Enable verbose output in the linear solver")
set(Xyce_VERBOSE_NONLINEAR         FALSE CACHE BOOL "Enable verbose output in the nonlinear solver")
set(Xyce_VERBOSE_NOX               FALSE CACHE BOOL "Enable verbose output in the NOX nonlinear solver library")
set(Xyce_VERBOSE_TIME              FALSE CACHE BOOL "Enable verbose output in the time integrator")

# Debug output
set(Xyce_DEBUG_ALL_PROCS_SAME_WD   FALSE CACHE BOOL "Enable debug output for the ")
set(Xyce_DEBUG_ANALYSIS            FALSE CACHE BOOL "Enable debug output for the Analysis package")
set(Xyce_DEBUG_CIRCUIT             FALSE CACHE BOOL "Enable debug output for the Circuit package")
set(Xyce_DEBUG_CONDUCTANCE         FALSE CACHE BOOL "Enable debug output for the ???")
set(Xyce_DEBUG_DEVICE              FALSE CACHE BOOL "Enable debug output for the Device package")
set(Xyce_DEBUG_DISTRIBUTION        FALSE CACHE BOOL "Enable debug output for the ???")
set(Xyce_DEBUG_EXCESS_PHASE        FALSE CACHE BOOL "Enable debug output for the ???")
set(Xyce_DEBUG_EXPRESSION          FALSE CACHE BOOL "Enable debug output for the Expression package")
set(Xyce_DEBUG_HB                  FALSE CACHE BOOL "Enable debug output for the Harmonic Balance package")
set(Xyce_DEBUG_IC                  FALSE CACHE BOOL "Enable debug output for the ???")
set(Xyce_DEBUG_IC_Gmin             FALSE CACHE BOOL "Enable debug output for the ???")
set(Xyce_DEBUG_IO                  FALSE CACHE BOOL "Enable debug output for the I/O Package")
set(Xyce_DEBUG_LINEAR              FALSE CACHE BOOL "Enable debug output for the Linear Solver package")
set(Xyce_DEBUG_MOR                 FALSE CACHE BOOL "Enable debug output for the Model Order Reduction package")
set(Xyce_DEBUG_MPDE                FALSE CACHE BOOL "Enable debug output for the MPDE Solver package")
set(Xyce_DEBUG_MPI                 FALSE CACHE BOOL "Enable debug output for the ???")
set(Xyce_DEBUG_NONLINEAR           FALSE CACHE BOOL "Enable debug output for the Nonlinear Solver package")
set(Xyce_DEBUG_OP_START            FALSE CACHE BOOL "Enable debug output for the ???")
set(Xyce_DEBUG_PARALLEL            FALSE CACHE BOOL "Enable debug output for the ???")
set(Xyce_DEBUG_RESTART             FALSE CACHE BOOL "Enable debug output for the Restart capability")
set(Xyce_DEBUG_TIME                FALSE CACHE BOOL "Enable debug output for the Time Integration package")
set(Xyce_DEBUG_TOPOLOGY            FALSE CACHE BOOL "Enable debug output for the Topology package")
set(Xyce_DEBUG_UNITS               FALSE CACHE BOOL "Enable debug output for the ???")
set(Xyce_DEBUG_VOLTLIM             FALSE CACHE BOOL "Enable debug output for the Voltage Limiting capability")
set(Xyce_Dakota_Debug              FALSE CACHE BOOL "Enable debug output for the Dakota interface")
set(Xyce_Dakota_Parallel_Debug     FALSE CACHE BOOL "Enable debug output for the Dakota interface in parallel")
set(Xyce_GRAPH_DEBUG               FALSE CACHE BOOL "Enable debug output for the ???")

# Troubleshooting
set(Xyce_DEBUG_TESTJAC             FALSE CACHE BOOL "Enable debug output for the ???")
set(Xyce_TEST_SOLN_VAR_MAP         FALSE CACHE BOOL "Enable debug output for the ???")


# Include the optional directories

if ((Xyce_ADMS_MODELS OR NOT DEFINED Xyce_ADMS_MODELS) AND EXISTS "${PROJECT_SOURCE_DIR}/src/DeviceModelPKG/ADMS")
     set (Xyce_ADMS_MODELS TRUE CACHE BOOL "Include the ADMS directory, if it exists")
     message(STATUS "Including the src/DeviceModelPKG/ADMS model directory")
elseif (NOT DEFINED Xyce_ADMS_MODELS)
     # The flag was not set, and the directory does not exist
     # Silently add the flag
     set (Xyce_ADMS_MODELS FALSE CACHE BOOL "Include the ADMS directory, if it exists")
elseif (NOT Xyce_ADMS_MODELS)
     # The flag was set to FALSE; the directory may or may not exist
     message(STATUS "NOT including the src/DeviceModelPKG/ADMS model directory")
     set (Xyce_ADMS_MODELS FALSE CACHE BOOL "Include the ADMS directory, if it exists")
else ()
     # The flag was set to TRUE, but the directory doesn't exist
     set (Xyce_ADMS_MODELS FALSE CACHE BOOL "Include the ADMS directory, if it exists" FORCE)
     message("The src/DeviceModelPKG/ADMS directory is not in place - "
          "changing Xyce_ADMS_MODELS to FALSE")
endif ()

if ((Xyce_NEURON_MODELS OR NOT DEFINED Xyce_NEURON_MODELS) AND EXISTS "${PROJECT_SOURCE_DIR}/src/DeviceModelPKG/NeuronModels")
     set (Xyce_NEURON_MODELS TRUE CACHE BOOL "Include the Neuron directory, if it exists")
     message(STATUS "Including the src/DeviceModelPKG/NeuronModels model directory")
elseif (NOT DEFINED Xyce_NEURON_MODELS)
     # The flag was not set, and the directory does not exist
     # Silently add the flag
     set (Xyce_NEURON_MODELS FALSE CACHE BOOL "Include the Neuron directory, if it exists")
elseif (NOT Xyce_NEURON_MODELS)
     # The flag was set to FALSE; the directory may or may not exist
     set (Xyce_NEURON_MODELS FALSE CACHE BOOL "Include the Neuron directory, if it exists")
     message(STATUS "NOT including the src/DeviceModelPKG/NeuronModels model directory")
else ()
     # The flag was set to TRUE, but the directory doesn't exist
     set (Xyce_NEURON_MODELS FALSE CACHE BOOL "Include the Neuron directory, if it exists" FORCE)
     message("The src/DeviceModelPKG/NeuronModels directory is not in place - "
          "changing Xyce_NEURON_MODELS to FALSE.")
endif ()

if ((Xyce_NONFREE_MODELS OR NOT DEFINED Xyce_NONFREE_MODELS) AND EXISTS "${PROJECT_SOURCE_DIR}/src/DeviceModelPKG/Xyce_NonFree")
     set (Xyce_NONFREE_MODELS TRUE CACHE BOOL "Include the NonFree directory, if it exists")
     message(STATUS "Including the src/DeviceModelPKG/Xyce_NonFree model directory")
elseif (NOT DEFINED Xyce_NONFREE_MODELS)
     # The flag was not set, and the directory does not exist
     # Silently add the flag (is this even necessary?)
     set (Xyce_NONFREE_MODELS FALSE CACHE BOOL "Include the NonFree directory, if it exists")
elseif (NOT Xyce_NONFREE_MODELS)
     # The flag was set to FALSE; the directory may or may not exist
     set (Xyce_NONFREE_MODELS FALSE CACHE BOOL "Include the NonFree directory, if it exists")
     message(STATUS "NOT including the src/DeviceModelPKG/Xyce_NonFree model directory")
else ()
     # The flag was set to TRUE, but the directory doesn't exist
     set (Xyce_NONFREE_MODELS FALSE CACHE BOOL "Include the NonFree directory, if it exists" FORCE)
     message("The src/DeviceModelPKG/NonFree directory is not in place - "
          "changing Xyce_NONFREE_MODELS to FALSE.")
endif ()

if ((Xyce_RAD_MODELS OR NOT DEFINED Xyce_RAD_MODELS) AND EXISTS "${PROJECT_SOURCE_DIR}/src/DeviceModelPKG/SandiaModels")
     set (Xyce_RAD_MODELS TRUE CACHE BOOL "Include the SandiaModels directory, if it exists")
     message(STATUS "Including the src/DeviceModelPKG/SandiaModels model directory")
elseif (NOT DEFINED Xyce_RAD_MODELS)
     # The flag was not set, and the directory does not exist
     # Silently add the flag (is this even necessary?)
     set (Xyce_RAD_MODELS FALSE CACHE BOOL "Include the SandiaModels directory, if it exists")
elseif (NOT Xyce_RAD_MODELS)
     # The flag was set to FALSE; the directory may or may not exist
     message(STATUS "NOT including the src/DeviceModelPKG/SandiaModels model directory")
     set (Xyce_RAD_MODELS FALSE CACHE BOOL "Include the SandiaModels directory, if it exists")
else ()
     # The flag was set to TRUE, but the directory doesn't exist
     set (Xyce_RAD_MODELS FALSE CACHE BOOL "Include the SandiaModels directory, if it exists" FORCE)
     message("The src/DeviceModelPKG/SandiaModels directory is not in place - "
          "changing Xyce_RAD_MODELS to FALSE.")
endif ()

# This logic could be in tps.cmake, since it's looking for Boost; but this is
# fine for now
if (Xyce_RAD_MODELS)
     if ((Xyce_ATHENA OR NOT DEFINED Xyce_ATHENA) AND EXISTS "${PROJECT_SOURCE_DIR}/src/DeviceModelPKG/SandiaModels/ATHENA")
          message(STATUS "Looking for Boost")
          find_package(Boost)
          if (Boost_FOUND)
               message(STATUS "Enabling the ATHENA model")
               set (Xyce_ATHENA TRUE CACHE BOOL "Include the ATHENA model, if it exists")
          else ()
               message("Boost was not found - disabling the ATHENA model")
               set (Xyce_ATHENA FALSE CACHE BOOL "Include the ATHENA model, if it exists" FORCE)
               set (Boost_INCLUDE_DIRS "")
          endif ()
     elseif (NOT DEFINED Xyce_ATHENA)
          # The flag was not set, and the directory does not exist
          # Silently add the flag
          set (Xyce_ATHENA FALSE CACHE BOOL "Include the ATHENA model, if it exists")
     elseif (NOT Xyce_ATHENA)
          # The flag was set to FALSE; the directory may or may not exist
          message(STATUS "NOT including the src/DeviceModelPKG/SandiaModels/ATHENA model directory")
          set (Xyce_ATHENA FALSE CACHE BOOL "Include the ATHENA model, if it exists")
     else ()
          # The flag was set to TRUE, but the directory doesn't exist
          set (Xyce_ATHENA FALSE CACHE BOOL "Include the ATHENA model, if it exists" FORCE)
          message("The src/DeviceModelPKG/SandiaModels/ATHENA directory is not in place - "
               "changing Xyce_ATHENA to FALSE.")
     endif ()
endif ()

