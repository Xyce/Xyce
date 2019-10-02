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


// ALL OF THE FOLLOWING IS OLD, BUT I'M NOT READY TO REMOVE IT, YET


/* Define to 1 if `lex' declares `yytext' as a `char *' by default, not a
   `char[]'. */
#cmakedefine YYTEXT_POINTER

/* Define to 1 if SRenvironment.h header is found */
#cmakedefine HAVE_SRENVIRONMENT_H


/* -- Xyce configuration ---------------------------------------------------- */


// -- version string ----------------
// This doesn't work when values are 0...
#cmakedefine Xyce_VERSION_MAJOR @Xyce_VERSION_MAJOR@
#cmakedefine Xyce_VERSION_MINOR @Xyce_VERSION_MINOR@
#cmakedefine Xyce_VERSION_PATCH @Xyce_VERSION_PATCH@
#cmakedefine Xyce_VERSION_TWEAK @Xyce_VERSION_TWEAK@
// So maybe use something like this?
#define X_V_MAJOR @Xyce_VERSION_MAJOR@
#define X_V_MINOR @Xyce_VERSION_MINOR@
#define X_V_PATCH @Xyce_VERSION_PATCH@

#endif 
