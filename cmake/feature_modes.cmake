
# These are all C++11 features that have ifdef's in the code.
# Setting them true, always
set (HAVE_IOTA TRUE)
set (HAVE_UNORDERED_MAP TRUE)
set (HAVE_UNORDERED_SET TRUE)
set (HAVE_ERF TRUE)
set (HAVE_ERFC TRUE)
# These may be problematic on Windows????????
set (HAVE_ISNAN TRUE)
set (HAVE_ISINF TRUE)

# END C++11 features

# This should always be true, so why is it optional?
set(Xyce_BELOS TRUE)

# For communicating the exact version of Xyce to the code
set(Xyce_RELEASE                   FALSE CACHE BOOL "Set to TRUE to designate a release version")
set(Xyce_QUALIFICATION             FALSE CACHE BOOL "Set to TRUE to designate a qualification release")
if(Xyce_RELEASE)
     set(RELEASE_CHARACTER "R")
     if(Xyce_QUALIFICATION)
          set(QUALIFICATION_CHARACTER "Q")
     endif()
else()
     set(RELEASE_CHARACTER "D")
endif()
set(Xyce_XyceVERSION "${RELEASE_CHARACTER}:${QUALIFICATION_CHARACTER}:${Xyce_VERSION_MAJOR}.${Xyce_VERSION_MINOR}.${Xyce_VERSION_PATCH}")

# Enable parallel
set(Xyce_PARALLEL_MPI              FALSE CACHE BOOL "Build Xyce with MPI enabled")

# Tracking capabilities
set(Xyce_USE_CURL                  FALSE CACHE BOOL "Enable the usage tracking capability using CURL")
set(Xyce_TRACKING_URL              ""  CACHE STRING "The URL for the usage tracking capability")

# Speed up compiles by disabling the analytic sensitivities in the
# ADMS-generated devices
set(Xyce_ADMS_SENSITIVITIES        TRUE CACHE BOOL "Enable analytic sensitivities in ADMS-generated devices")

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

if (EXISTS "${PROJECT_SOURCE_DIR}/src/DeviceModelPKG/ADMS" AND ( (NOT DEFINED Xyce_ADMS_MODELS) OR Xyce_ADMS_MODELS) )
     set (Xyce_ADMS_MODELS TRUE CACHE BOOL "Include the ADMS directory, if it exists")
     message(STATUS "Including the src/DeviceModelPKG/ADMS model directory")
elseif (NOT DEFINED Xyce_ADMS_MODELS)
     set (Xyce_ADMS_MODELS FALSE CACHE BOOL "Include the ADMS directory, if it exists")
elseif (NOT Xyce_ADMS_MODELS)
     message(STATUS "NOT including the src/DeviceModelPKG/ADMS model directory")
     set (Xyce_ADMS_MODELS FALSE CACHE BOOL "Include the ADMS directory, if it exists")
else ()
     set (Xyce_ADMS_MODELS FALSE CACHE BOOL "Include the ADMS directory, if it exists" FORCE)
     message("The src/DeviceModelPKG/ADMS directory is not in place - "
          "changing Xyce_ADMS_MODELS to FALSE")
endif ()

if (EXISTS "${PROJECT_SOURCE_DIR}/src/DeviceModelPKG/NeuronModels" AND ( (NOT DEFINED Xyce_NEURON_MODELS) OR Xyce_NEURON_MODELS) )
     set (Xyce_NEURON_MODELS TRUE CACHE BOOL "Include the Neuron directory, if it exists")
     message(STATUS "Including the src/DeviceModelPKG/NeuronModels model directory")
elseif (NOT DEFINED Xyce_NEURON_MODELS)
     set (Xyce_NEURON_MODELS FALSE CACHE BOOL "Include the Neuron directory, if it exists")
elseif (NOT Xyce_NEURON_MODELS)
     set (Xyce_NEURON_MODELS FALSE CACHE BOOL "Include the Neuron directory, if it exists")
     message(STATUS "NOT including the src/DeviceModelPKG/NeuronModels model directory")
else ()
     set (Xyce_NEURON_MODELS FALSE CACHE BOOL "Include the Neuron directory, if it exists" FORCE)
     message("The src/DeviceModelPKG/NeuronModels directory is not in place - "
          "changing Xyce_NEURON_MODELS to FALSE.")
endif ()

if (EXISTS "${PROJECT_SOURCE_DIR}/src/DeviceModelPKG/Xyce_NonFree" AND ( (NOT DEFINED Xyce_NONFREE_MODELS) OR Xyce_NONFREE_MODELS) )
     set (Xyce_NONFREE_MODELS TRUE CACHE BOOL "Include the NonFree directory, if it exists")
     message(STATUS "Including the src/DeviceModelPKG/Xyce_NonFree model directory")
elseif (NOT DEFINED Xyce_NONFREE_MODELS)
     set (Xyce_NONFREE_MODELS FALSE CACHE BOOL "Include the NonFree directory, if it exists")
elseif (NOT Xyce_NONFREE_MODELS)
     set (Xyce_NONFREE_MODELS FALSE CACHE BOOL "Include the NonFree directory, if it exists")
     message(STATUS "NOT including the src/DeviceModelPKG/Xyce_NonFree model directory")
else ()
     set (Xyce_NONFREE_MODELS FALSE CACHE BOOL "Include the NonFree directory, if it exists" FORCE)
     message("The src/DeviceModelPKG/NonFree directory is not in place - "
          "changing Xyce_NONFREE_MODELS to FALSE.")
endif ()

if (EXISTS "${PROJECT_SOURCE_DIR}/src/DeviceModelPKG/SandiaModels" AND ( (NOT DEFINED Xyce_RAD_MODELS) OR Xyce_RAD_MODELS) )
     set (Xyce_RAD_MODELS TRUE CACHE BOOL "Include the SandiaModels directory, if it exists")
     message(STATUS "Including the src/DeviceModelPKG/SandiaModels model directory")
elseif (NOT DEFINED Xyce_RAD_MODELS)
     set (Xyce_RAD_MODELS FALSE CACHE BOOL "Include the SandiaModels directory, if it exists")
elseif (NOT Xyce_RAD_MODELS)
     message(STATUS "NOT including the src/DeviceModelPKG/SandiaModels model directory")
     set (Xyce_RAD_MODELS FALSE CACHE BOOL "Include the SandiaModels directory, if it exists")
else ()
     set (Xyce_RAD_MODELS FALSE CACHE BOOL "Include the SandiaModels directory, if it exists" FORCE)
     message("The src/DeviceModelPKG/SandiaModels directory is not in place - "
          "changing Xyce_RAD_MODELS to FALSE.")
endif ()


