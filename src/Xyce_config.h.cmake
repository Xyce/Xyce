
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
#cmakedefine HAVE_ALGORITHM

/* Define to 1 if you have the <algo.h> header file. */
#cmakedefine HAVE_ALGO_H

/* define if you have the bison parser generator */
#cmakedefine HAVE_BISON

/* Define to 1 if you have the <cctype> header file. */
#cmakedefine HAVE_CCTYPE

/* Define to 1 if you have the <climits> header file. */
#cmakedefine HAVE_CLIMITS

/* Define to 1 if you have the <cmath> header file. */
#cmakedefine HAVE_CMATH

/* Define to 1 if you have the <cstdio> header file. */
#cmakedefine HAVE_CSTDIO

/* Define to 1 if you have the <cstdlib> header file. */
#cmakedefine HAVE_CSTDLIB

/* Define to 1 if you have the <cstring> header file. */
#cmakedefine HAVE_CSTRING

/* Define to 1 if you have the <strings.h> header file. */
#cmakedefine HAVE_STRINGS_H

/* Define to 1 if you have the <string.h> header file. */
#cmakedefine HAVE_STRING_H

/* Define to 1 if you have the <dlfcn.h> header file. */
#cmakedefine HAVE_DLFCN_H

/* define if you have the flex lexical scanner generator */
#cmakedefine HAVE_FLEX

/* define if the Standard Template Library has flexible insert */
#cmakedefine HAVE_FLEXIBLE_INSERT

/* Define to 1 if you have the <float.h> header file. */
#cmakedefine HAVE_FLOAT_H

/* Define to 1 if you have the <fstream> header file. */
#cmakedefine HAVE_FSTREAM

/* Define to 1 if you have the `getopt' function. */
#cmakedefine HAVE_GETOPT

/* Define to 1 if you have the `DRAND48' function. */
#cmakedefine HAVE_DRAND48

/* Define to 1 if you have the <inttypes.h> header file. */
#cmakedefine HAVE_INTTYPES_H

/* Define to 1 if you have the <iostream> header file. */
#cmakedefine HAVE_IOSTREAM

/* define if the Standard Template Library algorithms include iota */
#cmakedefine HAVE_IOTA

/* Define to 1 if you have the `amd' library (-lamd). */
#cmakedefine HAVE_LIBAMD

/* Define to 1 if you have the `expr' library (-lexpr). */
#cmakedefine HAVE_LIBEXPR

/* Define to 1 if you have the `metis' library (-lmetis). */
#cmakedefine HAVE_LIBMETIS

/* Define to 1 if you have the `parmetis' library (-lparmetis). */
#cmakedefine HAVE_LIBPARMETIS

/* Define to 1 if you have the `superlu' library (-lsuperlu). */
#cmakedefine HAVE_LIBSUPERLU

/* Define to 1 if you have the `superludist' library (-lsuperludist). */
#cmakedefine HAVE_LIBSUPERLUDIST

/* Define to 1 if you have the `umfpack' library (-lumfpack). */
#cmakedefine HAVE_LIBUMFPACK

/* Define to 1 if you have the `y12m' library (-ly12m). */
#cmakedefine HAVE_LIBY12M

/* Define to 1 if you have the `zoltan' library (-lzoltan). */
#cmakedefine HAVE_LIBZOLTAN

/* Define to 1 if you have the <limits.h> header file. */
#cmakedefine HAVE_LIMITS_H

/* Define to 1 if you have the 'mallinfo' function. */
#cmakedefine HAVE_MALLINFO

/* Define to 1 if you have the <malloc.h> header file. */
#cmakedefine HAVE_MALLOC_H

/* Define to 1 if you have the <math.h> header file. */
#cmakedefine HAVE_MATH_H

/* Define to 1 if you have the <mathimf.h> header file. */
#cmakedefine HAVE_MATHIMF_H

/* Define to 1 if you have the <memory.h> header file. */
#cmakedefine HAVE_MEMORY_H

/* define if the compiler implements namespaces */
#cmakedefine HAVE_NAMESPACES

/* define if the compiler supports isnan and isinf checks */
#cmakedefine HAVE_NAN_INF_SUPPORT

/* define if the compiler supports _isnan and _finite checks */
#cmakedefine HAVE__ISNAN_AND__FINITE_SUPPORT


/* Define to 1 if you have the <ostream> header file. */
#cmakedefine HAVE_OSTREAM

/* Define to 1 if you have the <pwd.h> header file. */
#cmakedefine HAVE_PWD_H

/* Define to 1 if you have the <ostream.h> header file. */
#cmakedefine HAVE_OSTREAM_H

/* Define to 1 if you have the <stdint.h> header file. */
#cmakedefine HAVE_STDINT_H

/* Define to 1 if you have the <stdlib.h> header file. */
#cmakedefine HAVE_STDLIB_H

/* define if the compiler supports Standard Template Library */
#cmakedefine HAVE_STL

/* Define to 1 if you have the <strings.h> header file. */
#cmakedefine HAVE_STRINGS_H

/* Define to 1 if you have the <string.h> header file. */
#cmakedefine HAVE_STRING_H

/* Define to 1 if you have the <Windows.h> header file. */
#cmakedefine HAVE_WINDOWS_H

/* Define to 1 if you have the <sys/resource.h> header file. */
#cmakedefine HAVE_SYS_RESOURCE_H

/* Define to 1 if you have the <functional> header file. */
#cmakedefine HAVE_FUNCTIONAL

/* Define to 1 if you have the 'getdomainname' function. */
#cmakedefine HAVE_GETDOMAINNAME

/* Define to 1 if you have the 'gethostname' function. */
#cmakedefine HAVE_GETHOSTNAME

/* Define to 1 if you have the 'getpwuid' function. */
#cmakedefine HAVE_GETPWUID

/* Define to 1 if you have the 'erf' function. */
#cmakedefine HAVE_ERF

/* Define to 1 if you have the 'erf' function. */
#cmakedefine HAVE_ERFC

/* Define to 1 if you have the <tr1/functional> header file. */
#cmakedefine HAVE_TR1_FUNCTIONAL

/* Define to 1 if you have the <unordered_map> header file. */
#cmakedefine HAVE_UNORDERED_MAP

/* Define to 1 if you have the <tr1/unordered_map> header file. */
#cmakedefine HAVE_TR1_UNORDERED_MAP

/* Define to 1 if you have the <unordered_set> header file. */
#cmakedefine HAVE_UNORDERED_SET

/* Define to 1 if you have the <tr1/unordered_set> header file. */
#cmakedefine HAVE_TR1_UNORDERED_SET

/* define if the compiler has strcasecmp function */
#cmakedefine HAVE_STRCASECMP

/* Define to 1 if you have the <sys/stat.h> header file. */
#cmakedefine HAVE_SYS_STAT_H

/* Define to 1 if you have the <sys/types.h> header file. */
#cmakedefine HAVE_SYS_TYPES_H

/* Define to 1 if you have the <sys/utsname.h> header file. */
#cmakedefine HAVE_SYS_UTSNAME_H

/* Define to 1 if you have the <unistd.h> header file. */
#cmakedefine HAVE_UNISTD_H

/* Define to 1 if you have the <values.h> header file. */
#cmakedefine HAVE_VALUES_H

/* Define to 1 if you have the <Epetra_MultiVector.h> header file. */
#cmakedefine HAVE_EPETRA_MULTIVECTOR_H

/* Define to 1 if you have the <Ifpack_CrsRiluk.h> header file. */
#cmakedefine HAVE_IFPACK_CRSRILUK_H

/* Define to 1 if you have the <Amesos_Klu.h> header file. */
#cmakedefine HAVE_AMESOS_KLU_H

/* Define to 1 if you have the <Amesos_Umfpack.h> header file. */
#cmakedefine HAVE_AMESOS_UMFPACK_H

/* Define to 1 if you have the <Amesos_Superlu.h> header file. */
#cmakedefine HAVE_AMESOS_SUPERLU_H

/* Define to 1 if you have the <NOX_Abstract_Vector.h> header file. */
#cmakedefine HAVE_NOX_ABSTRACT_VECTOR_H

/* Define to 1 if you have the <LOCA_Parameter_Vector.h> header file. */
#cmakedefine HAVE_LOCA_PARAMETER_VECTOR_H

/* define if math.h defines M_PI */
#cmakedefine MATH_H_HAS_M_PI

/* define if stl string header does not cause pair.h to be included */
#cmakedefine NEED_PAIR_H

/*
  These are carried over from autotools.  Don't really need them except
  that N_UTL_Version.C references them
*/
/* Name of package */
#cmakedefine PACKAGE

/* Define to the address where bug reports for this package should be sent. */
#cmakedefine PACKAGE_BUGREPORT

/* Define to the full name of this package. */
#cmakedefine PACKAGE_NAME

/* Define to the full name and version of this package. */
#cmakedefine PACKAGE_STRING

/* Define to the one symbol short name of this package. */
#cmakedefine PACKAGE_TARNAME



/* DEBUG:  Define to the version of this package. */
#cmakedefine PACKAGE_VERSION

/* Define to 1 if you have the ANSI C header files. */
#cmakedefine STDC_HEADERS


/* Define to 1 if `lex' declares `yytext' as a `char *' by default, not a
   `char[]'. */
#cmakedefine YYTEXT_POINTER

/* Define to 1 if unistd.h header is not found */
#cmakedefine YY_NO_UNISTD_H

/* Define to 1 if SRenvironment.h header is found */
#cmakedefine HAVE_SRENVIRONMENT_H 

#cmakedefine HAVE_LINUX_EXCEPTIONS

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
#cmakedefine Xyce_PARALLEL_MPI

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
