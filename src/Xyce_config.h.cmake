#ifndef  Xyce_CONFIG_H
#define  Xyce_CONFIG_H

// C++11 Features; these should go away ASAP

#cmakedefine HAVE_IOTA
#cmakedefine HAVE_UNORDERED_MAP
#cmakedefine HAVE_UNORDERED_SET
#cmakedefine HAVE_ERF
#cmakedefine HAVE_ERFC
#cmakedefine HAVE_ISNAN
#cmakedefine HAVE_ISINF

// This is required (so why is it ifdef'd?)
#cmakedefine Xyce_BELOS

// Compile Xyce for MPI parallelism
#cmakedefine Xyce_PARALLEL_MPI

// Trilinos
#cmakedefine Xyce_SHYLU
#cmakedefine SHYLUBASKER
#cmakedefine Xyce_AMESOS2
#cmakedefine Xyce_STOKHOS_ENABLE
#cmakedefine Xyce_ROL

#cmakedefine Xyce_USE_ISORROPIA

// Trilinos TPLs
#cmakedefine Xyce_USE_PARMETIS
#cmakedefine Xyce_AMD

// FFT
#cmakedefine Xyce_USE_FFT
#cmakedefine Xyce_USE_INTEL_FFT
#cmakedefine Xyce_USE_FFTW

// Reaction parser
#cmakedefine Xyce_REACTION_PARSER

// Enable usage tracking
#cmakedefine Xyce_USE_CURL

// Dakota coupling
#cmakedefine Xyce_Dakota

// Windows features
#cmakedefine HAVE_WINDOWS_H
#cmakedefine HAVE_WIN_DIRCOMMANDS
#cmakedefine YY_NO_UNISTD_H
#cmakedefine _CRT_SECURE_NO_WARNINGS



// POSIX features
#cmakedefine HAVE_DLFCN_H
#cmakedefine HAVE_DRAND48
#cmakedefine HAVE_UNISTD_H
#cmakedefine HAVE_SYS_RESOURCE_H
#cmakedefine HAVE_SYS_STAT_H

// Build directories
#cmakedefine Xyce_ADMS_MODELS
#cmakedefine Xyce_NONFREE_MODELS
#cmakedefine Xyce_NEURON_MODELS
#cmakedefine Xyce_RAD_MODELS
#cmakedefine Xyce_ATHENA
// This is not in the CMake at the moment
#cmakedefine Xyce_MODSPEC_MODELS

// Enable analytic sensitivities in ADMS-generated devices
#cmakedefine Xyce_ADMS_SENSITIVITIES
#cmakedefine01 Xyce_STOKHOS_ENABLE

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
#cmakedefine Xyce_TEST_SOLN_VAR_MAP

// Set the Xyce version for the code
#define VERSION "@Xyce_XyceVERSION@"

#endif 
