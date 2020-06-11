
#ifndef  Xyce_CONFIG_H
#define  Xyce_CONFIG_H




/* 
  This flag is used to signal some code in N_UTL_Version.C which
  expcects that we're using autotools and not cmake 
 */
#cmakedefine USING_CMAKE


/* 
  Define this to trick trilinos include files into using package_config.h 
  files.  Without this a lot of trilinos headers will faile. 
*/
#cmakedefine HAVE_CONFIG_H


/*
  These are conditionally handled by cmake.  I've added the ifndef's to avoid
  warning about redefinging names as some of these may be defined in the
  package_include.h files from Trilinos
*/

/* Define to 1 if you have the <algorithm> header file. */
#ifndef HAVE_ALGORITHM
#cmakedefine HAVE_ALGORITHM
#endif 

/* Define to 1 if you have the <algo.h> header file. */
#ifndef HAVE_ALGO_H
#cmakedefine HAVE_ALGO_H
#endif

/* define if you have the bison parser generator */
#ifndef HAVE_BISON
#cmakedefine HAVE_BISON
#endif

/* Define to 1 if you have the <dlfcn.h> header file. */
#ifndef HAVE_DLFCN_H
#cmakedefine HAVE_DLFCN_H
#endif

/* define if you have the flex lexical scanner generator */
#ifndef HAVE_FLEX
#cmakedefine HAVE_FLEX
#endif

/* define if the Standard Template Library has flexible insert */
#ifndef HAVE_FLEXIBLE_INSERT
#cmakedefine HAVE_FLEXIBLE_INSERT
#endif

/* Define to 1 if you have the <float.h> header file. */
#ifndef HAVE_FLOAT_H
#cmakedefine HAVE_FLOAT_H
#endif

/* Define to 1 if you have the `getopt' function. */
#ifndef HAVE_GETOPT
#cmakedefine HAVE_GETOPT
#endif

/* Define to 1 if you have the `DRAND48' function. */
#ifndef HAVE_DRAND48
#cmakedefine HAVE_DRAND48
#endif

/* Define to 1 if you have the <inttypes.h> header file. */
#ifndef HAVE_INTTYPES_H
#cmakedefine HAVE_INTTYPES_H
#endif

/* define if the Standard Template Library algorithms include iota */
#ifndef HAVE_IOTA
#cmakedefine HAVE_IOTA
#endif

/* Define to 1 if you have the `amd' library (-lamd). */
#ifndef HAVE_LIBAMD
#cmakedefine HAVE_LIBAMD
#endif

/* Define to 1 if you have the `expr' library (-lexpr). */
#ifndef HAVE_LIBEXPR
#cmakedefine HAVE_LIBEXPR
#endif

/* Define to 1 if you have the `metis' library (-lmetis). */
#ifndef HAVE_LIBMETIS
#cmakedefine HAVE_LIBMETIS
#endif

/* Define to 1 if you have the `parmetis' library (-lparmetis). */
#ifndef HAVE_LIBPARMETIS
#cmakedefine HAVE_LIBPARMETIS
#endif

/* Define to 1 if you have the `superlu' library (-lsuperlu). */
#ifndef HAVE_LIBSUPERLU
#cmakedefine HAVE_LIBSUPERLU
#endif

/* Define to 1 if you have the `superludist' library (-lsuperludist). */
#ifndef HAVE_LIBSUPERLUDIST
#cmakedefine HAVE_LIBSUPERLUDIST
#endif

/* Define to 1 if you have the `umfpack' library (-lumfpack). */
#ifndef HAVE_LIBUMFPACK
#cmakedefine HAVE_LIBUMFPACK
#endif

/* Define to 1 if you have the `y12m' library (-ly12m). */
#ifndef HAVE_LIBY12M
#cmakedefine HAVE_LIBY12M
#endif

/* Define to 1 if you have the `zoltan' library (-lzoltan). */
#ifndef HAVE_LIBZOLTAN
#cmakedefine HAVE_LIBZOLTAN
#endif

/* Define to 1 if you have the <limits.h> header file. */
#ifndef HAVE_LIMITS_H 
#cmakedefine HAVE_LIMITS_H
#endif

/* Define to 1 if you have the 'mallinfo' function. */
#ifndef HAVE_MALLINFO
#cmakedefine HAVE_MALLINFO
#endif

/* Define to 1 if you have the <malloc.h> header file. */
#ifndef HAVE_MALLOC_H 
#cmakedefine HAVE_MALLOC_H
#endif

/* Define to 1 if you have the <mathimf.h> header file. */
#ifndef HAVE_MATHIMF_H
#cmakedefine HAVE_MATHIMF_H
#endif

/* Define to 1 if you have the <memory.h> header file. */
#ifndef HAVE_MEMORY_H
#cmakedefine HAVE_MEMORY_H
#endif

/* define if the compiler implements namespaces */
#ifndef HAVE_NAMESPACES
#cmakedefine HAVE_NAMESPACES
#endif

/* define if the compiler supports isnan and isinf checks */
#ifndef HAVE_NAN_INF_SUPPORT
#cmakedefine HAVE_NAN_INF_SUPPORT
#endif

/* define if the compiler supports _isnan and _finite checks */
#ifndef HAVE__ISNAN_AND__FINITE_SUPPORT
#cmakedefine HAVE__ISNAN_AND__FINITE_SUPPORT
#endif

/* Define to 1 if you have the <pwd.h> header file. */
#ifndef HAVE_PWD_H 
#cmakedefine HAVE_PWD_H
#endif

/* Define to 1 if you have the <stdint.h> header file. */
#ifndef HAVE_STDINT_H
#cmakedefine HAVE_STDINT_H
#endif

/* Define to 1 if you have the <stdlib.h> header file. */
#ifndef HAVE_STDLIB_H
#cmakedefine HAVE_STDLIB_H
#endif

/* define if the compiler supports Standard Template Library */
#ifndef HAVE_STL
#cmakedefine HAVE_STL
#endif

/* Define to 1 if you have the <strings.h> header file. */
#ifndef HAVE_STRINGS_H
#cmakedefine HAVE_STRINGS_H
#endif

/* Define to 1 if you have the <string.h> header file. */
#ifndef HAVE_STRING_H
#cmakedefine HAVE_STRING_H
#endif

/* Define to 1 if you have the <Windows.h> header file. */
#ifndef HAVE_WINDOWS_H
#cmakedefine HAVE_WINDOWS_H
#endif

/* Define to 1 if you have the <sys/resource.h> header file. */
#ifndef HAVE_SYS_RESOURCE_H
#cmakedefine HAVE_SYS_RESOURCE_H
#endif

/* Define to 1 if you have the <functional> header file. */
#ifndef HAVE_FUNCTIONAL
#cmakedefine HAVE_FUNCTIONAL
#endif

/* Define to 1 if you have the 'getdomainname' function. */
#ifndef HAVE_GETDOMAINNAME
#cmakedefine HAVE_GETDOMAINNAME
#endif

/* Define to 1 if you have the 'gethostname' function. */
#ifndef HAVE_GETHOSTNAME
#cmakedefine HAVE_GETHOSTNAME
#endif

/* Define to 1 if you have the 'getpwuid' function. */
#ifndef HAVE_GETPWUID
#cmakedefine HAVE_GETPWUID
#endif

/* Define to 1 if you have the 'erf' function. */
#ifndef HAVE_ERF
#cmakedefine HAVE_ERF
#endif

/* Define to 1 if you have the 'erf' function. */
#ifndef HAVE_ERFC
#cmakedefine HAVE_ERFC
#endif

/* Define to 1 if you have the <tr1/functional> header file. */
#ifndef HAVE_TR1_FUNCTIONAL
#cmakedefine HAVE_TR1_FUNCTIONAL
#endif

/* Define to 1 if you have the <unordered_map> header file. */
#ifndef HAVE_UNORDERED_MAP
#cmakedefine HAVE_UNORDERED_MAP
#endif

/* Define to 1 if you have the <tr1/unordered_map> header file. */
#ifndef HAVE_TR1_UNORDERED_MAP
#cmakedefine HAVE_TR1_UNORDERED_MAP
#endif

/* Define to 1 if you have the <unordered_set> header file. */
#ifndef HAVE_UNORDERED_SET
#cmakedefine HAVE_UNORDERED_SET
#endif

/* Define to 1 if you have the <tr1/unordered_set> header file. */
#ifndef HAVE_TR1_UNORDERED_SET
#cmakedefine HAVE_TR1_UNORDERED_SET
#endif

/* define if the compiler has strcasecmp function */
#ifndef HAVE_STRCASECMP
#cmakedefine HAVE_STRCASECMP
#endif

/* Define to 1 if you have the <sys/stat.h> header file. */
#ifndef HAVE_SYS_STAT_H
#cmakedefine HAVE_SYS_STAT_H
#endif

/* Define to 1 if you have the <sys/types.h> header file. */
#ifndef HAVE_SYS_TYPES_H
#cmakedefine HAVE_SYS_TYPES_H
#endif

/* Define to 1 if you have the <sys/utsname.h> header file. */
#ifndef HAVE_SYS_UTSNAME_H 
#cmakedefine HAVE_SYS_UTSNAME_H
#endif

/* Define to 1 if you have the <unistd.h> header file. */
#ifndef HAVE_UNISTD_H
#cmakedefine HAVE_UNISTD_H
#endif

/* Define to 1 if you have the <values.h> header file. */
#ifndef HAVE_VALUES_H
#cmakedefine HAVE_VALUES_H
#endif

/* Define to 1 if you have the <Epetra_MultiVector.h> header file. */
#ifndef HAVE_EPETRA_MULTIVECTOR_H
#cmakedefine HAVE_EPETRA_MULTIVECTOR_H
#endif

/* Define to 1 if you have the <Ifpack_CrsRiluk.h> header file. */
#ifndef HAVE_IFPACK_CRSRILUK_H
#cmakedefine HAVE_IFPACK_CRSRILUK_H
#endif

/* Define to 1 if you have the <Amesos_Klu.h> header file. */
#ifndef HAVE_AMESOS_KLU_H
#cmakedefine HAVE_AMESOS_KLU_H
#endif

/* Define to 1 if you have the <Amesos_Umfpack.h> header file. */
#ifndef HAVE_AMESOS_UMFPACK_H
#cmakedefine HAVE_AMESOS_UMFPACK_H
#endif

/* Define to 1 if you have the <Amesos_Superlu.h> header file. */
#ifndef HAVE_AMESOS_SUPERLU_H
#cmakedefine HAVE_AMESOS_SUPERLU_H
#endif

/* Define to 1 if you have the <NOX_Abstract_Vector.h> header file. */
#ifndef HAVE_NOX_ABSTRACT_VECTOR_H
#cmakedefine HAVE_NOX_ABSTRACT_VECTOR_H
#endif

/* Define to 1 if you have the <LOCA_Parameter_Vector.h> header file. */
#ifndef HAVE_LOCA_PARAMETER_VECTOR_H
#cmakedefine HAVE_LOCA_PARAMETER_VECTOR_H
#endif

/* define if stl string header does not cause pair.h to be included */
#ifndef NEED_PAIR_H
#cmakedefine NEED_PAIR_H
#endif


/*
  These are carried over from autotools.  Don't really need them except
  that N_UTL_Version.C references them
*/
/* Name of package */
#ifndef PACKAGE
#cmakedefine PACKAGE
#endif

/* Define to the address where bug reports for this package should be sent. */
#ifndef PACKAGE_BUGREPORT
#cmakedefine PACKAGE_BUGREPORT
#endif

/* Define to the full name of this package. */
#ifndef PACKAGE_NAME
#cmakedefine PACKAGE_NAME
#endif

/* Define to the full name and version of this package. */
#ifndef PACKAGE_STRING
#cmakedefine PACKAGE_STRING
#endif

/* Define to the one symbol short name of this package. */
#ifndef PACKAGE_TARNAME
#cmakedefine PACKAGE_TARNAME
#endif



/* DEBUG:  Define to the version of this package. */
#ifndef PACKAGE_VERSION
#cmakedefine PACKAGE_VERSION
#endif

/* Define to 1 if `lex' declares `yytext' as a `char *' by default, not a
   `char[]'. */
#ifndef YYTEXT_POINTER
#cmakedefine YYTEXT_POINTER
#endif

/* Define to 1 if unistd.h header is not found */
#ifndef YY_NO_UNISTD_H 
#cmakedefine YY_NO_UNISTD_H
#endif

/* Define to 1 if SRenvironment.h header is found */
#ifndef HAVE_SRENVIRONMENT_H 
#cmakedefine HAVE_SRENVIRONMENT_H 
#endif



/* -- Xyce configuration ---------------------------------------------------- */


/* -- version string ---------------- */
#cmakedefine Xyce_VERSION_MAJOR @Xyce_VERSION_MAJOR@
#cmakedefine Xyce_VERSION_MINOR @Xyce_VERSION_MINOR@
#cmakedefine Xyce_VERSION_PATCH @Xyce_VERSION_PATCH@
#cmakedefine Xyce_VERSION_EXTRA @Xyce_VERSION_EXTRA@


/* -- build options ----------------- */
#cmakedefine Xyce_AMD
#cmakedefine Xyce_BELOS
#cmakedefine Xyce_BTF
#cmakedefine Xyce_CHARON
#cmakedefine Xyce_Dakota
#cmakedefine Xyce_Dakota50
#cmakedefine Xyce_DEPENDENCY_TRACKING
#cmakedefine Xyce_EXTDEV 
#cmakedefine Xyce_NEW_DAE_FORMULATION 

// The autotools build system only passes this in on the commandline.
#ifndef Xyce_PARALLEL_MPI
#cmakedefine Xyce_PARALLEL_MPI
#endif

#cmakedefine Xyce_PARDISO_MKL 
#cmakedefine Xyce_RAD_MODELS
#cmakedefine Xyce_NONFREE_MODELS
#cmakedefine Xyce_ADMS_MODELS
#cmakedefine Xyce_NEURON_MODELS
#cmakedefine Xyce_REACTION_PARSER
#cmakedefine Xyce_SENSITIVITY_ENABLE
#cmakedefine Xyce_SHYLU
#cmakedefine Xyce_SPICE_NORMS
#cmakedefine Xyce_SUPERLU   
#cmakedefine Xyce_SUPERLUDIST 
#cmakedefine Xyce_TRILINOS_DEV
#cmakedefine Xyce_UMFPACK         
#cmakedefine Xyce_USE_BSIM3_CONST
#cmakedefine Xyce_USE_ISORROPIA
#cmakedefine Xyce_USE_ZOLTAN
#cmakedefine Xyce_USE_PARMETIS
#cmakedefine Xyce_USE_HDF5
#cmakedefine Xyce_USE_CURL
#cmakedefine Xyce_ATHENA
#cmakedefine HAVE_LIBARAENV
#cmakedefine Xyce_ADMS_SENSITIVITIES
#cmakedefine01 Xyce_STOKHOS_ENABLE

/* -- verbosity options ------------- */
#cmakedefine Xyce_VERBOSE_LINEAR
#cmakedefine Xyce_VERBOSE_NONLINEAR
#cmakedefine Xyce_VERBOSE_NOX
#cmakedefine Xyce_VERBOSE_TIME


/* -- debugging options ------------- */
#cmakedefine Xyce_DEBUG
#cmakedefine Xyce_DEBUG_ANALYSIS
#cmakedefine Xyce_DEBUG_CIRCUIT
#cmakedefine Xyce_DEBUG_DEVICE
#cmakedefine Xyce_DEBUG_DIRECTSOLVE
#cmakedefine Xyce_DEBUG_DISTRIBUTION
#cmakedefine Xyce_DEBUG_EXPRESSION
#cmakedefine Xyce_DEBUG_IO    
#cmakedefine Xyce_DEBUG_LINEAR
#cmakedefine Xyce_DEBUG_NONLINEAR
#cmakedefine Xyce_DEBUG_PARALLEL
#cmakedefine Xyce_DEBUG_RESTART
#cmakedefine Xyce_DEBUG_TIME  
#cmakedefine Xyce_DEBUG_TOPOLOGY
#cmakedefine Xyce_TEST_SOLN_VAR_MAP


/* -- internally set options -------- */

#cmakedefine Xyce_NOX_LOCA_SUPPORT  
#cmakedefine Xyce_RESTART_NOPACK
#cmakedefine Xyce_USE_INTEL_FFT
#cmakedefine Xyce_USE_FFTW
#cmakedefine Xyce_USE_FFT 




#endif 
