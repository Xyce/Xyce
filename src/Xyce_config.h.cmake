#ifndef  Xyce_CONFIG_H
#define  Xyce_CONFIG_H

// Compile Xyce for MPI parallelism
#cmakedefine Xyce_PARALLEL_MPI

// Trilinos
#cmakedefine Xyce_SHYLU
#cmakedefine Xyce_AMESOS2
#cmakedefine Xyce_AMESOS2_KLU2
#cmakedefine Xyce_AMESOS2_BASKER
#cmakedefine Xyce_AMESOS2_SHYLUBASKER
#cmakedefine Xyce_NEW_BASKER
#cmakedefine Xyce_STOKHOS_ENABLE
#cmakedefine Xyce_ROL
#cmakedefine Xyce_USE_ISORROPIA

// Trilinos TPLs
#cmakedefine Xyce_AMD

// FFT
#cmakedefine Xyce_USE_FFT
#cmakedefine Xyce_USE_INTEL_FFT
#cmakedefine Xyce_USE_FFTW

// Binary Restart files
#cmakedefine Xyce_RESTART_NOPACK

// Enable usage tracking
#cmakedefine Xyce_USE_CURL
#cmakedefine Xyce_TRACKING_URL "@Xyce_TRACKING_URL@"

// Dakota coupling
#cmakedefine Xyce_Dakota

// Windows features
#cmakedefine HAVE_WINDOWS_H
#cmakedefine HAVE_WIN_DIRCOMMANDS
#cmakedefine _CRT_SECURE_NO_WARNINGS
#cmakedefine HAVE_DIRECT_H

// POSIX features
#cmakedefine HAVE_DLFCN_H
#cmakedefine HAVE_DRAND48
#cmakedefine HAVE_UNISTD_H
#cmakedefine HAVE_SYS_RESOURCE_H
#cmakedefine HAVE_SYS_STAT_H
#cmakedefine HAVE_MALLOC_H
#cmakedefine HAVE_MALLINFO
#cmakedefine HAVE_PWD_H
#cmakedefine HAVE_GETPWUID
#cmakedefine HAVE_GETHOSTNAME
#cmakedefine HAVE_GETDOMAINNAME
#cmakedefine HAVE_SYS_UTSNAME_H
#cmakedefine HAVE_UNAME

// Reaction parser
#cmakedefine Xyce_REACTION_PARSER
#cmakedefine YY_NO_UNISTD_H

// Build directories
#cmakedefine Xyce_ADMS_MODELS
#cmakedefine Xyce_NONFREE_MODELS
#cmakedefine Xyce_NEURON_MODELS
#cmakedefine Xyce_RAD_MODELS
#cmakedefine Xyce_ATHENA
// This is not in the CMake at the moment
#cmakedefine Xyce_MODSPEC_MODELS

// Add the NOX solver statistics, if available
#cmakedefine Xyce_NOX_SOLVERSTATS

// Enable analytic sensitivities in ADMS-generated devices
#cmakedefine Xyce_ADMS_SENSITIVITIES

// Support for Charon coupling
#cmakedefine Xyce_CHARON

// Verbose output
#cmakedefine Xyce_VERBOSE_CONDUCTANCE
#cmakedefine Xyce_VERBOSE_LINEAR
#cmakedefine Xyce_VERBOSE_NONLINEAR
#cmakedefine Xyce_VERBOSE_NOX
#cmakedefine Xyce_VERBOSE_TIME

// Debug output
#cmakedefine Xyce_DEBUG_ALL_PROCS_SAME_WD
#cmakedefine Xyce_DEBUG_ANALYSIS
#cmakedefine Xyce_DEBUG_CIRCUIT
#cmakedefine Xyce_DEBUG_CONDUCTANCE
#cmakedefine Xyce_DEBUG_DEVICE
#cmakedefine Xyce_DEBUG_DISTRIBUTION
#cmakedefine Xyce_DEBUG_EXCESS_PHASE
#cmakedefine Xyce_DEBUG_EXPRESSION
#cmakedefine Xyce_DEBUG_HB
#cmakedefine Xyce_DEBUG_IC
#cmakedefine Xyce_DEBUG_IC_Gmin
#cmakedefine Xyce_DEBUG_IO
#cmakedefine Xyce_DEBUG_LINEAR
#cmakedefine Xyce_DEBUG_MOR
#cmakedefine Xyce_DEBUG_MPDE
#cmakedefine Xyce_DEBUG_MPI
#cmakedefine Xyce_DEBUG_NONLINEAR
#cmakedefine Xyce_DEBUG_OP_START
#cmakedefine Xyce_DEBUG_PARALLEL
#cmakedefine Xyce_DEBUG_RESTART
#cmakedefine Xyce_DEBUG_TIME
#cmakedefine Xyce_DEBUG_TOPOLOGY
#cmakedefine Xyce_DEBUG_UNITS
#cmakedefine Xyce_DEBUG_VOLTLIM
#cmakedefine Xyce_Dakota_Debug
#cmakedefine Xyce_Dakota_Parallel_Debug
#cmakedefine Xyce_GRAPH_DEBUG

// Troubleshooting
#cmakedefine Xyce_DEBUG_TESTJAC

// Use constants from BSIM3 device 
#cmakedefine Xyce_USE_BSIM3_CONST

// Set the Xyce version for the code
#define VERSION "@Xyce_XyceVERSION@"

#endif 
