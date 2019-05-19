dnl Process this file with autoconf to produce a configure script.
dnl NEVER edit the configure script itself, as it will be overwritten
dnl when autoconf processes this file.

dnl ********************************************************************
dnl Autoconf boilerplate
dnl Initialize the macros, define the package name, version number and 
dnl bug reporting address
AC_INIT(Xyce,D::7.0,xyce@sandia.gov)

dnl This is a FILTHY HACK designed to get rid of the troublesome space that
dnl AC_INIT stuffs into PACKAGE_STRING, and which buggers MPILAM's "hcp" 
dnl script.

sed -e 's/PACKAGE_STRING\(.*\) \(.*\)/PACKAGE_STRING\1_\2/g' < confdefs.h > foobie.h && mv foobie.h confdefs.h
PACKAGE_STRING=`echo $PACKAGE_STRING|sed -e 's/ /_/g'`

dnl Identify a single file within the source tree to make sure we're 
dnl configuring the right source
AC_CONFIG_SRCDIR(src/CircuitPKG/N_CIR_Xyce.C)

dnl Auxillary Build Tools
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR(config)

# Configure should create Xyce_config.h from Xyce_config.h.in
AC_CONFIG_HEADERS(src/Xyce_config.h)

AM_SILENT_RULES(yes)

dnl Require version 2.59 or higher of autoconf
AC_PREREQ(2.59)

dnl Determine a canonicalized system type (e.g. i386-unknown-freebsd4.4)
AC_CANONICAL_TARGET

dnl Initialize Automake
AM_INIT_AUTOMAKE([1.11 foreign])

dnl ********************************************************************
dnl the autoconfing proper 

dnl *****************************************************************
dnl Before we can get to the meat, we need to do some options
dnl *****************************************************************


dnl ***********************************************************************
dnl The Xyce project has a directory full of pre-built libraries for
dnl each architecture.  Use it if we can.
dnl If the environment variable ARCH_ROOT is set, use it, 
dnl otherwise pick a safe default
dnl Allow user to specify an "ARCH_ROOT", and kludge a -L for it
if test "x$ARCH_ROOT" = x; then
 if test -d /Net/Proj/Xyce/arch; then
   AC_MSG_NOTICE([Using /Net/Proj/Xyce/arch as default ARCH_ROOT])
   ARCH_ROOT=/Net/Proj/Xyce/arch
 fi
fi

dnl ********************************************
dnl Check to see if the executable is disabled.
dnl RPP: 03/15/2006 - For Xyce/Charon coupling,
dnl we need to turn off building of the executable,
dnl autotools is screwing the link line in Nevada's TPL
AC_ARG_ENABLE(xyce-binary,
[AS_HELP_STRING([--disable-xyce-binary],[Disables building of the xyce executable. (Default is to build the binary)])],
[ENABLE_XYCE_BINARY=$enableval],
[ENABLE_XYCE_BINARY=yes]
)

AC_MSG_CHECKING(whether to build the xyce binary)
if test "X$ENABLE_XYCE_BINARY" != "Xno"; then
  AC_MSG_RESULT(yes)
else
  CPPFLAGS="-DDINT $CPPFLAGS"
  AC_MSG_RESULT(no)
fi

AC_ARG_ENABLE(xyce-shareable,
[AS_HELP_STRING([--enable-xyce-shareable],[Enables building of the shareable xyce executable. (Default is to build the non-shareable)])],
[ENABLE_XYCE_SHAREABLE=$enableval],
[ENABLE_XYCE_SHAREABLE=no]
)

AC_MSG_CHECKING(whether to build xyce shareable)
if test "X$ENABLE_XYCE_SHAREABLE" != "Xno"; then
  AC_MSG_RESULT(yes)
  MAYBE_SCRIPTS="buildxyceplugin"
  MAYBE_LIBEXEC_SCRIPTS="libtool"
  MAYBE_UTILS="utils"
else
  AC_MSG_RESULT(no)
  MAYBE_SCRIPTS=""
  MAYBE_UTILS=""
  MAYBE_LIBEXEC_SCRIPTS=""
fi

AM_CONDITIONAL(BUILD_XYCE_SHAREABLE, test "X$ENABLE_XYCE_SHAREABLE" != "Xno")
AC_SUBST([MAYBE_SCRIPTS])
AC_SUBST([MAYBE_LIBEXEC_SCRIPTS])
AC_SUBST([MAYBE_UTILS])

# --enable-user-plugin means that "make all" will build the plugin if
# we have --enable-xyce-shareable, and make install will install it.
# Without it, the user has to do "make plugin" and "make install-plugin"

MAYBE_PLUGIN=""
AC_ARG_ENABLE(user-plugin,
[AS_HELP_STRING([--enable-user-plugin],[enables building of the demo verilog plug-in as part of "make all" when doing shared executable build. You must have ADMS installed for this to work. (Default is only to build the plugin with "make plugin")])],
[MAYBE_PLUGIN="user_plugin"],
[MAYBE_PLUGIN=""]
)
AC_SUBST([MAYBE_PLUGIN])



AM_CONDITIONAL(BUILD_XYCE_BINARY, test "X$ENABLE_XYCE_BINARY" != "Xno")
dnl ********************************************

dnl ********************************************
dnl Here we determine whether we're doing a threaded build
AC_ARG_ENABLE(thread-safe-logging, AC_HELP_STRING([--enable-thread-safe-logging], 
                                      [enable thread safe logging for one Xyce per thread runs (default is no)]), 
                                      ac_enable_threads=$enableval, 
                                      ac_enable_threads=no) 
if test "$ac_enable_threads" = yes; then 
  AC_DEFINE(USE_THREADS, 1, [use threads]) 
fi

AM_CONDITIONAL(THREADS, test $ac_enable_threads = yes) 
dnl ********************************************

dnl ********************************************
dnl Here we determine whether we're doing a parallel build or not
XYCE_DEBUG_OPTION(mpi,no,[parallel build with MPI],Xyce_PARALLEL_MPI,USE_MPI)

if test "$USE_MPI" = yes; then
   AC_DEFINE([Xyce_PARALLEL_MPI], [1], [define if running with MPI parallel]) 
fi

#Lastly, we allow the user to completely override our choice of Xyce_ARCH
# if need be.

dnl And now we know precisely what build we're doing, so maybe we can use
dnl ARCH_ROOT
dnl but if user has specified ARCHDIR use it instead of $ARCH_ROOT/$Xyce_ARCH
if test "x$ARCHDIR" = "x"; then
 if test "x$ARCH_ROOT" != "x"; then
  if test "x$Xyce_ARCH" != "x"; then
    if test -d $ARCH_ROOT/$Xyce_ARCH; then
      ARCHDIR=$ARCH_ROOT/$Xyce_ARCH
    fi
  fi  
 fi
fi


dnl see comments at header of XYCE_DEBUG_OPTION in acinclude.m4 for
dnl usage
XYCE_DEBUG_OPTION(verbose_time,no,[verbosity in time integrator])
XYCE_DEBUG_OPTION(verbose_nonlinear,no,[verbosity in nonlinear solver])
XYCE_DEBUG_OPTION(verbose_nox,no,[verbosity in NOX nonlinear solver library])
XYCE_DEBUG_OPTION(verbose_linear,no,[verbosity in linear solver])
XYCE_DEBUG_OPTION(spice_norms,yes,[using SPICE type norms])
XYCE_DEBUG_OPTION(bsim3_const,no,[using constants from BSIM3 instead of more precise constants consistent with other devices],Xyce_USE_BSIM3_CONST)
XYCE_DEBUG_OPTION(nox_use_vector_copy,yes,[using vector copies in NOX])

dnl optionally compile in support for rad models if the directory is there and
dnl the user hasn't specified disable-radmodels
MAYBE_RAD=""
if test -d $srcdir/src/DeviceModelPKG/SandiaModels; then
  XYCE_DEBUG_OPTION(radmodels,yes,[Radiation Model Library],Xyce_RAD_MODELS,RAD_MODELS)

  if test "x$RAD_MODELS" = "xyes"
  then
    MAYBE_RAD="SandiaModels"
  fi
else
  RAD_MODELS=no
fi
AM_CONDITIONAL(RADMODELS,test $RAD_MODELS = yes)
AC_SUBST([MAYBE_RAD])

dnl optionally compile in support for adms models if the directory is there and
dnl the user hasn't specified disable-admsmodels
MAYBE_ADMS=""
if test -d $srcdir/src/DeviceModelPKG/ADMS; then
  XYCE_DEBUG_OPTION(admsmodels,yes,[ADMS Model Library],Xyce_ADMS_MODELS,ADMS_MODELS)

  if test "x$ADMS_MODELS" = "xyes"
  then
    MAYBE_ADMS="ADMS"
  fi
else
  ADMS_MODELS=no
fi
AM_CONDITIONAL(ADMSMODELS,test $ADMS_MODELS = yes)
AC_SUBST([MAYBE_ADMS])

dnl optionally compile in support for neuron models if the directory is there and
dnl the user hasn't specified disable-neuronmodels
MAYBE_NEURON=""
if test -d $srcdir/src/DeviceModelPKG/NeuronModels; then
  XYCE_DEBUG_OPTION(neuronmodels,yes,[NEURON Model Library],Xyce_NEURON_MODELS,NEURON_MODELS)

  if test "x$NEURON_MODELS" = "xyes"
  then
    MAYBE_NEURON="NeuronModels"
  fi
else
  NEURON_MODELS=no
fi
AM_CONDITIONAL(NEURONMODELS,test $NEURON_MODELS = yes)
AC_SUBST([MAYBE_NEURON])

dnl optionally compile in support for non-free models if the directory is there and
dnl the user hasn't specified disable-nonfree
MAYBE_NONFREE=""    
if test -d $srcdir/src/DeviceModelPKG/Xyce_NonFree; then
  XYCE_DEBUG_OPTION(nonfree,yes,[Non-Free Models Library],Xyce_NONFREE_MODELS,NONFREE_MODELS)
  if test "x$NONFREE_MODELS" = "xyes" 
  then
    MAYBE_NONFREE="Xyce_NonFree"
  fi
else
  NONFREE_MODELS=no
fi
AM_CONDITIONAL(NONFREEMODELS,test $NONFREE_MODELS = yes)
AC_SUBST([MAYBE_NONFREE])

dnl optionally compile in support for ModSpec models if the directory is there and
dnl the user hasn't specified disable-modspec
MAYBE_BERKELEY=""    
if test -d $srcdir/src/DeviceModelPKG/BERKELEY; then
  XYCE_DEBUG_OPTION(modspec,yes,[ModSpec Interface],Xyce_MODSPEC_MODELS,MODSPEC_MODELS)
  if test "x$MODSPEC_MODELS" = "xyes" 
  then
    MAYBE_BERKELEY="BERKELEY"
  fi
else
  MODSPEC_MODELS=no
fi

AM_CONDITIONAL(MODSPECMODELS,test $MODSPEC_MODELS = yes)
AC_SUBST([MAYBE_BERKELEY])

dnl optionally compile in support for vpi tests if user has specified enable-vpitests
dnl the user hasn't specified disable-vpitests
MAYBE_VPITESTS=""
VPI_BASE=""
AC_CHECK_PROG([HAS_IVERILOG],[iverilog],[yes],[no])
if test "x$HAS_IVERILOG" = "xyes"; then
  AC_MSG_NOTICE([Enabling VPI tests.])
  VPI_BASE=$(dirname $(dirname $(which iverilog)))
  MAYBE_VPITESTS="VPITests"
fi
AC_SUBST([MAYBE_VPITESTS])
AC_SUBST([VPI_BASE])

dnl optionally compile in support for Dakota linkage
XYCE_DEBUG_OPTION(dakota,no,[Dakota direct linkage support],Xyce_Dakota, USE_DAKOTA)
AM_CONDITIONAL(DAKOTA_LINK, test $USE_DAKOTA = yes)
dnl For the LHS package in Dakota we need to link in a Dakota object file as well as
dnl the dakota libraries.  
AM_CONDITIONAL(DAKOTA_OBJ_NEEDED, test $USE_DAKOTA = yes )

dnl optionally compile in support for Charon linkage
XYCE_DEBUG_OPTION(charon,no,[Charon device support],Xyce_CHARON, USE_CHARON)
AM_CONDITIONAL(BUILD_CHARON_SUPPORT, test $USE_CHARON = yes)

dnl optionally compile in support for Trilinos Development
XYCE_DEBUG_OPTION(trilinos_dev,no,[support for development (pre-release) versions of Trilinos],Xyce_TRILINOS_DEV, USE_TRILINOS_DEV)
AM_CONDITIONAL(TRILINOS_DEV, test $USE_TRILINOS_DEV = yes)

dnl compile in support for partitioning
XYCE_DEBUG_OPTION(isorropia,yes,[isorropia partitioning option],Xyce_USE_ISORROPIA,USE_ISORROPIA)
XYCE_DEBUG_OPTION(zoltan,yes,[zoltan partitioning option],Xyce_USE_ZOLTAN,USE_ZOLTAN)

dnl optionally compile in support for SuperLU and UmfPack.
XYCE_DEBUG_OPTION(superlu,no,[SuperLU direct solver],Xyce_SUPERLU,SUPERLU)
XYCE_DEBUG_OPTION(pardiso_mkl,no,[Pardiso direct solver (through MKL)],Xyce_PARDISO_MKL,PARDISO_MKL)

dnl optionally compile in ksparse solver
XYCE_DEBUG_OPTION(ksparse,yes,[Kundert's Sparse direct solver],Xyce_KSPARSE,USE_KSPARSE)
AM_CONDITIONAL(BUILD_KSPARSE, test "x$USE_KSPARSE" = "xyes")

dnl optionally compile in support for Belos (as of Xyce 6.0, this is "yes" by default)
XYCE_DEBUG_OPTION(belos,yes,[Belos iterative solver support],Xyce_BELOS, USE_BELOS)
AM_CONDITIONAL(BUILD_BELOS_SUPPORT, test "x$USE_BELOS" = "xyes")

dnl optionally compile in support for ShyLU (we need Trilinos development branch at this time)
XYCE_DEBUG_OPTION(shylu,no,[ShyLU hybrid solver support],Xyce_SHYLU, USE_SHYLU)
AM_CONDITIONAL(BUILD_SHYLU_SUPPORT, test "x$USE_SHYLU" = "xyes")

dnl optionally compile in support for Amesos2 (we need Trilinos development branch at this time)
XYCE_DEBUG_OPTION(amesos2,no,[Amesos2 direct solver support],Xyce_AMESOS2, USE_AMESOS2)
AM_CONDITIONAL(BUILD_AMESOS2_SUPPORT, test "x$USE_AMESOS2" = "xyes")

dnl optionally compile in support for ROL (we need Trilinos development branch at this time)
XYCE_DEBUG_OPTION(rol,no,[ROL optimization library support],Xyce_ROL, USE_ROL)
AM_CONDITIONAL(BUILD_ROL_SUPPORT, test "x$USE_ROL" = "xyes")

dnl optionally compile in support for Stokhos
XYCE_DEBUG_OPTION(stokhos,no,[Stokhos UQ library support],Xyce_STOKHOS_ENABLE, USE_STOKHOS)
AM_CONDITIONAL(BUILD_STOKHOS, test "x$USE_STOKHOS" = "xyes")

dnl Must include AMD
XYCE_DEBUG_OPTION(amd,yes,[AMD reordering],Xyce_AMD,AMD)

if test "$USE_MPI" = "yes"; then
  XYCE_DEBUG_OPTION(superludist,no,[SuperLU_Dist linear solver],Xyce_SUPERLUDIST,SUPERLUDIST)
else
  SUPERLUDIST=no
fi

XYCE_DEBUG_OPTION(test_soln_var_map,no,[mapping of internal variables to external for debugging])
XYCE_DEBUG_OPTION(debug_device,no,[device package debugging output])
XYCE_DEBUG_OPTION(debug_io,no,[I/O package debugging])
XYCE_DEBUG_OPTION(debug_expression,no,[Expression package debugging output])
XYCE_DEBUG_OPTION(debug_time,no,[time integrator package debugging output])
XYCE_DEBUG_OPTION(debug_analysis,no,[analysis package debugging output])
XYCE_DEBUG_OPTION(debug_hb,no,[Harmonic balance analysis package debugging output])
XYCE_DEBUG_OPTION(debug_es,no,[embedded sampling analysis package debugging output])
XYCE_DEBUG_OPTION(debug_sampling,no,[traditional sampling analysis package debugging output])
XYCE_DEBUG_OPTION(debug_mor,no,[Model order reduction debugging output])
XYCE_DEBUG_OPTION(debug_parallel,no,[parallel distribution package debugging output])
XYCE_DEBUG_OPTION(debug_distribution,no,[distribution package debugging output])
XYCE_DEBUG_OPTION(debug_topology,no,[topology package debugging output])
XYCE_DEBUG_OPTION(debug_linear,no,[linear solver package debugging output])
XYCE_DEBUG_OPTION(debug_nonlinear,no,[nonlinear solver package debugging output])
XYCE_DEBUG_OPTION(debug_circuit,no,[circuit package debugging output])
XYCE_DEBUG_OPTION(debug_directsolve,no,[direct solver package debugging output])
XYCE_DEBUG_OPTION(debug_restart,no,[restart debugging output])

XYCE_DEBUG_OPTION(op_start,yes,[allow startup from previous operating point])

dnl check if the developer wants to skip the probe of the fortran compiler
XYCE_DEBUG_OPTION(fortran_blas,yes,[test for fortran blas using fortran compiler --- disable this if you are using C blas],Xyce_FORTRAN_BLAS,FORTRAN_BLAS)

# Allow disabling of reaction parser
XYCE_DEBUG_OPTION(reaction_parser,yes,[reaction parser], Xyce_REACTION_PARSER,REACTION_PARSER)

# Allow disabling of FFT libraries
XYCE_DEBUG_OPTION(fft,yes,[support for FFT libraries],Xyce_FFT, USE_FFT)
if test "$USE_FFT" = "yes"; then
  XYCE_DEBUG_OPTION(intel_fft,yes,[support for Intel MKL FFT library],Xyce_IntelFFT, USE_IntelFFT)
  XYCE_DEBUG_OPTION(fftw,yes,[support for FFTW library],Xyce_FFTW, USE_FFTW)
else
  USE_FFTW=no
  USE_IntelFFT=no
fi

XYCE_DEBUG_OPTION([athena], [yes], [building the ATHENA device],Xyce_ATHENA_REQUSTED,USE_ATHENA)


XYCE_DEBUG_OPTION([adms_sensitivities], [yes], [analytic sensitivities in ADMS-generated device models],Xyce_ADMS_SENSITIVITIES)


# Default to HDF5 support off 
# Note:  a fourth argument *must* be given, because XYCE_DEBUG_OPTION 
# forces a "-D" into CPPFLAGS.  But because this particular option may be
# disabled later when we check if HDF5 libraries and headers are there, the
# fourth argument must not be the same one that Xyce uses to test for HDF5.
# That is, the fourth argument is just "user has asked for this feature"
# not "this feature works"
XYCE_DEBUG_OPTION(hdf5,no,[support for HDF5 output library],Xyce_HDF5, USE_HDF5)

# Default to curl support off
# see commentary above re: HDF5 and use of multiple symbols here.
# do NOT use "Xyce_CURL" in code, because it will be defined if cURL requested
# even if cURL isn't usable.
XYCE_DEBUG_OPTION(curl,no,[support for curl for metrics reporting],Xyce_CURL, USE_CURL)

# If cURL, set the default TRACKING_URL
if test "X$USE_CURL" = "Xyes"; then
   if test "X$TRACKING_URL" = "X"; then
      TRACKING_URL=http://127.0.0.1/
   fi
   AC_DEFINE_UNQUOTED([Xyce_TRACKING_URL], ["$TRACKING_URL"], [Audit Xyce runs.])
fi

dnl Checks for programs.
# Checks for programs.
AC_PROG_CXX(CC c++ g++ gcc cxx icpc)
AC_PROG_CC(cc gcc icc)
if test "$FORTRAN_BLAS" != "no"; then
  AC_PROG_F77(f77 g77 g95 gfortran ifort)
fi


if test "$REACTION_PARSER" = "yes"; then
 # check for flex --- we use flex-only features, and can't use vanilla lex
 # if don't find flex, don't build parser
 # Same goes for bison, as yacc won't do the job for us, and AC_PROG_YACC
 # will use bison in yacc compatibility mode if we use it.  WE also
 # cannot just check for "all versions of bison newer than X.Y" because
 # bison's C++ interface has a habit of changing incompatibly from time
 # to time.  We therefore enumerate actual versions of Bison our interface
 # is *known* to work with and against which Xyce has been tested. 
 AC_PROG_LEX
 XYCE_CHECK_BISON([3.3, 3.2, 3.1, 3.0, 2.7, 2.6, 2.5, 2.4], , [BISON=no])
 if test `basename $LEX` = "flex";
 then
    AC_DEFINE(HAVE_FLEX,1,[define if you have the flex lexical scanner generator])
 else
    AC_MSG_ERROR(["Requested reaction parser, but no FLEX."])
 fi
 if test "$BISON" != "no";
 then
   AC_DEFINE(HAVE_BISON,1,[define if you have the bison parser generator])
 else
   AC_MSG_ERROR(["Requested reaction parser, but no bison of correct version number."])
 fi
else
 LEX=no
 BISON=no
fi

AM_CONDITIONAL(BUILD_REACTION_PARSER, test `basename $LEX` = "flex" -a "$BISON" != "no")

dnl Check to see if valid Fortran compiler was found (HKT; 07/11/06; bug #1091).
if test "$FORTRAN_BLAS" != "no" -a "x$F77" = "x" ; then
  AC_MSG_ERROR([Fortran BLAS testing requested but no valid Fortran compiler was found in your path.])
fi

AC_PROG_MAKE_SET
# static exec and lib built by default; use --enable-shared to get shared lib
LT_INIT([disable-shared])
AC_SUBST([LIBTOOL_DEPS])
AC_PROG_GREP

dnl *********************************************************************
dnl KLUDGES:
dnl The following kludges add defined constants that are highly machine
dnl and compiler dependent --- so let's do them after we figure out what 
dnl machine, build and compiler we're using

dnl needed to move this thing to AFTER AC_PROG_CXX because putting it before
dnl that test disables the automatic insertion of -g -O2 in the default build.

USING_ICCLT11_CXX=no
XYCE_TEST_USING_ICCLT11_CXX

# Let the user force ADMS_CXXFLAGS.  If not forcing, check 'em ourselves
if test "x$ADMS_CXXFLAGS" = "x"
then
  ADMS_CXXFLAGS="$CXXFLAGS"
  # special CXX flags for ADMS directory on the ICC version <=11 compilers
  if test "$USING_ICCLT11_CXX" != "no"; then
    ADMS_CXXFLAGS=`echo $CXXFLAGS | sed -e 's/-O3/-O1/' -e 's/-O2/-O1/'`
  fi
fi

dnl Include the top-level of the include and lib directories in the search
dnl Add to loader flags so comes before all -l's
if test -d $ARCHDIR/lib; then
  LDFLAGS="$LDFLAGS -L$ARCHDIR/lib"
fi
if test -d $ARCHDIR/include; then
  CPPFLAGS="-I$ARCHDIR/include $CPPFLAGS"
fi

dnl Include the trilinos include directory in search
if test -d $ARCHDIR/include/trilinos; then
  CPPFLAGS="-I$ARCHDIR/include/trilinos $CPPFLAGS"
fi

dnl Include Dakota headers and library paths if asked for Dakota build
if test "x$USE_DAKOTA" = "xyes"; then
  if test "x$DAKOTADIR" = "x"; then
    if test -d $ARCHDIR/dakota; then
      DAKOTA_ARCHDIR="$ARCHDIR/dakota"
    fi
  else
    DAKOTA_ARCHDIR="$DAKOTADIR"
  fi
  if test -d $DAKOTA_ARCHDIR/include; then
    default_INCS="$default_INCS -I$DAKOTA_ARCHDIR/include"
    CPPFLAGS="-I$DAKOTA_ARCHDIR/include $CPPFLAGS"
  fi
  dnl if test -d $DAKOTA_ARCHDIR/reference/lib; then
  dnl   LDFLAGS="$LDFLAGS -L$DAKOTA_ARCHDIR/reference/lib"
  dnl fi
  if test -d $DAKOTA_ARCHDIR/lib; then
    LDFLAGS="$LDFLAGS -L$DAKOTA_ARCHDIR/lib"
  fi
fi


dnl Add some FLIBS and LDFLAGS for g95 based on Eric Keiter's Intel-Mac (HKT: 1/24/2007; bug# 1098)
if test "x$F77" = "xg95"; then
  FLIBS="-L/sw/lib -lg95 $FLIBS"
  LDFLAGS="-Wl,-framework,Accelerate $LDFLAGS"
fi

dnl *********************************************************************
dnl checks for libraries now
dnl Replace `main' with a function in -lm:
if test "$FORTRAN_BLAS" != "no"; then
  AC_F77_LIBRARY_LDFLAGS
fi

AC_SEARCH_LIBS(sqrt, [m])

dnl find a valid blas and lapack library

if test "$FORTRAN_BLAS" != "no"; then
  AC_LANG_PUSH(Fortran 77)
  AC_SEARCH_LIBS(dcopy,[essl cxml blas_slu blas])
  AC_SEARCH_LIBS(dgetrs,[essl cxml scs complib.sgimath lapack cxml])
  AC_LANG_POP(Fortran 77)
else
  AC_SEARCH_LIBS(dcopy,[essl cxml blas_slu blas])
  AC_SEARCH_LIBS(dgetrs,[essl cxml scs complib.sgimath lapack cxml])
fi

if  test "$ac_cv_search_dcopy" = "no"
then
  AC_MSG_ERROR([Configure has not found BLAS libraries for your system.  Please install BLAS or provide appropriate -l and -L flags to LIBS and LDFLAGS to allow configure to locate the dcopy function.])
fi
if  test "$ac_cv_search_dgetrs" = "no"
then
  AC_MSG_ERROR([Configure has not found LAPACK libraries for your system.  Please install LAPACK or provide appropriate -l and -L flags to LIBS and LDFLAGS to allow configure to locate the dgetrs function.])
fi

#It should now be safe to add FLIBS to LIBS
LIBS="${LIBS} ${FLIBS}"

dnl *********************************************************************
dnl *****DO NOT ADD ANY TRILINOS PROBES ABOVE THIS LINE! ***************
dnl *** These basic C++ checks need to be very early, because we need them
dnl *** in order some of the probes to work on all platforms
dnl *********************************************************************
dnl This is a showstopper, and must be done first!
AC_CXX_HAVE_STL

AC_LANG_PUSH(C++)
AC_CHECK_HEADERS(unistd.h values.h limits.h float.h unordered_map unordered_set functional tr1/unordered_map tr1/unordered_set sys/resource.h)

dnl *********************************************************************
dnl dlopen interface
AC_CHECK_HEADERS([stdlib.h dlfcn.h])
AC_SEARCH_LIBS([dlopen], [dl])

dnl *********************************************************************
dnl stats information
AC_CHECK_HEADERS(unistd.h malloc.h pwd.h sys/utsname.h)
AC_CHECK_FUNCS(mallinfo)
AC_CHECK_FUNCS(getpwuid)
AC_CHECK_FUNCS(gethostname)
AC_CHECK_FUNCS(getdomainname)
AC_CHECK_FUNCS(uname)
AC_CHECK_FUNCS(sysconf)
AC_CHECK_FILES(/proc/self/stat)
AC_CHECK_HEADERS(mathimf.h)

dnl *********************************************************************
dnl See if compiler supplies a nan/inf checker.  This must be done after the
dnl cmath probe, because some platforms have wonky cmath/math.h pairs.
AC_FINITE_NUMBER_CHECK

dnl CAREFUL --- if parallel is enabled, epetra needs MPI, so better have
dnl it later in the libs path! (which means detecting it earlier, since
dnl LIBS is assembled right to left)

dnl and if this is a parallel build we need all the Zoltan stuff
if test "x$USE_MPI" = "xyes"; then
  AC_SEARCH_LIBS(MPI_Bcast,[mpi])
fi

dnl If the user enabled AMD, let's try to add it here:
if test "x$AMD" = "xyes"; then
  AC_CHECK_LIB(amd,amd_valid,[],[AC_MSG_ERROR([FATAL: AMD  requested but AMD library not found or does not contain amd_valid])])
fi

# Make sure we have a version of Trilinos newer than 11.2.3
XYCE_LIBRARY_OPTION(teuchoscore,[],[],[XYCE_NEW_TEUCHOS_FOUND=1],[])

# Now look for the new style if we haven't found the old.  Any one of these
# not being found is a fatal configure error.
if test $XYCE_NEW_TEUCHOS_FOUND = 1; then
  XYCE_TRILINOS_LIBRARY_OPTION(Teuchos, USE_TEUCHOS_MAKEFILE, yes) 
  if test "x$USE_TEUCHOS_MAKEFILE" = "xno"; then
    XYCE_LIBRARY_OPTION(teuchoscore)
    XYCE_LIBRARY_OPTION(teuchosremainder)
    XYCE_LIBRARY_OPTION(teuchosnumerics)
    XYCE_LIBRARY_OPTION(teuchoscomm)
    XYCE_LIBRARY_OPTION(teuchosparameterlist)
    # Newer versions of trilinos need this, if we find it, use it
    AC_CHECK_LIB(trilinosss,trilinos_klu_defaults,[],[])
  fi
else
  AC_MSG_ERROR(
[===========================================================])

  AC_MSG_ERROR(
[Configure detected a version of Trilinos prior to 11.2.3.  Update Trilinos libraries and configure again."])

  AC_MSG_ERROR(
[===========================================================])
fi

XYCE_TRILINOS_LIBRARY_OPTION(Triutils, USE_TRIUTILS_MAKEFILE)
XYCE_TRILINOS_LIBRARY_OPTION(Epetra, USE_EPETRA_MAKEFILE)
XYCE_TRILINOS_LIBRARY_OPTION(EpetraExt, USE_EPETRAEXT_MAKEFILE)

dnl Check whether an essential AMD-related Trilinos header exists
if test "x$AMD" = "xyes"; then
  AC_CHECK_HEADER(EpetraExt_AMD_CrsGraph.h,,[AC_MSG_ERROR([Sorry, you have requested AMD, but EpetraExt_AMD_CrsGraph.h was not found.  Either recompile Trilinos with AMD support (recommended), or use '--disable-amd' (not recommended).])],)
fi

dnl Xyce depends on Teuchos having been built with COMPLEX enabled.  In
dnl old days, this happened automatically, but now you have to specify
dnl special cmake options to make it happen.  When users reuse an old
dnl trilinos configuration line or fail to read the building guide, they
dnl sometimes miss this.  And when they do, they get cryptic compilation
dnl errors.  Let's check that NOW instead of at compile time.

AC_CHECK_DECL([HAVE_TEUCHOS_COMPLEX],,[AC_MSG_ERROR([Your Trilinos was not built with the correct options and does not have COMPLEX support built into Teuchos.  Please rebuild Trilinos using the options specified in the Xyce building guide.])],[[#include <Teuchos_config.h>]])

dnl Now we must check for whether EpetraExt defines HAVE_PARMETIS, in which
dnl case we must also link to parmetis and metis libraries if we're building
dnl parallel
if test "x$USE_MPI" = "xyes"; then
   AC_CHECK_DECL([HAVE_LIBPARMETIS],
                 [
                  AC_CHECK_LIB(metis,METIS_PartGraphKway,[],[AC_MSG_ERROR([FATAL: METIS library not found or not valid])])
                  AC_CHECK_LIB(parmetis,ParMETIS_PartKway,[],
                               [AC_CHECK_LIB(parmetis,ParMETIS_V3_PartKway,[],
                               [AC_MSG_ERROR([FATAL: ParMETIS library not found or not valid])])])
                  ],[],[[#include <EpetraExt_config.h>]])
  if test "x$USE_ZOLTAN" = "xyes" -o "x$USE_ISORROPIA" = "xyes"; then
    XYCE_LIBRARY_INCLUDE(zoltan)
    AC_CHECK_LIB(zoltan,Zoltan_LB_Balance,[],[AC_MSG_ERROR([FATAL: zoltan library not found or not valid])])
  fi
fi

dnl and if this is a parallel build we (may) need all the Isorropia stuff
if test "x$USE_MPI" = "xyes" -a "x$USE_ISORROPIA" = "xyes"; then
  XYCE_TRILINOS_LIBRARY_OPTION(Isorropia, USE_ISORROPIA_MAKEFILE)
fi

dnl If the user enabled superludist, let's try to add it here:
if test "x$SUPERLUDIST" = "xyes"; then
    XYCE_LIBRARY_OPTION(superludist)
fi

dnl If the user enabled superlu, let's try to add it here:
if test "x$SUPERLU" = "xyes"; then
  AC_CHECK_LIB(superlu,dgssvx,[],[AC_MSG_ERROR([FATAL: SuperLU library does not contain dgssvx and this is required by Trilinos.])])
fi

#This is obnoxious.  If we do this probe *AFTER* the code that adds Amesos
# SOME platforms don't get the amesos linkage right.
AC_CHECK_HEADER(Amesos_Umfpack.h,
 [

  # This may already have been done if AMD enabled, but if it wasn't,
  # we require it now.  If it was already done, this is fast because
  # the results were cached
  AC_CHECK_LIB([amd],[amd_valid],[],[AC_MSG_ERROR([FATAL: Umfpack requires AMD but AMD library not found or does not contain amd_valid])])

  # we found that Trilinos has been built with Umfpack, we must now
  # find the libraries for umfpack
  # Newer cholmods may depend on suitesparseconfig.  If we find it, add it.
  AC_CHECK_LIB(suitesparseconfig,SuiteSparse_time,[],[])
  # cholmod may depend on colamd...
  AC_CHECK_LIB(colamd,colamd,[],[])
  # Newer umfpacks depend on cholmod.  If we find it, add it
  AC_CHECK_LIB(cholmod,cholmod_start,[],[])
  AC_CHECK_LIB(umfpack,umfpack_di_solve,[],[AC_MSG_ERROR([FATAL: UmfPack library does not contain umfpack_di_solve and this is required by Trilinos.])])
 ],[])

XYCE_TRILINOS_LIBRARY_OPTION(Amesos, USE_AMESOS_MAKEFILE)
XYCE_TRILINOS_LIBRARY_OPTION(TrilinosCouplings, USE_TRILINOSCOUPLINGS_MAKEFILE)

XYCE_TRILINOS_LIBRARY_OPTION(Ifpack, USE_IFPACK_MAKEFILE)
XYCE_TRILINOS_LIBRARY_OPTION(AztecOO, USE_AZTECOO_MAKEFILE)

dnl If the user enabled Belos, let's add it here:
if test "x$USE_BELOS" = "xyes"; then
    XYCE_TRILINOS_LIBRARY_OPTION(Belos, USE_BELOS_MAKEFILE)
    if test "x$USE_BELOS_MAKEFILE" = "xno"; then
      XYCE_LIBRARY_OPTION(belosepetra)
    fi
fi

dnl If the user enabled ShyLU, let's add it here:
if test "x$USE_SHYLU" = "xyes"; then
    XYCE_TRILINOS_LIBRARY_OPTION(ShyLU, USE_SHYLU_MAKEFILE)
fi

dnl If the user enabled Amesos2, let's add it here:
if test "x$USE_AMESOS2" = "xyes"; then
    XYCE_TRILINOS_LIBRARY_OPTION(Amesos2, USE_AMESOS2_MAKEFILE)
fi

dnl If the user enabled ROL, let's add it here:
if test "x$USE_ROL" = "xyes"; then
    XYCE_TRILINOS_LIBRARY_OPTION(ROL, USE_ROL_MAKEFILE)
fi

dnl NOX brings in the LOCA libraries as well using the Makefile.export system
XYCE_TRILINOS_LIBRARY_OPTION(NOX, USE_NOX_MAKEFILE)
if test "x$USE_NOX_MAKEFILE" = "xno"; then
  XYCE_LIBRARY_OPTION(loca)
fi
CPPFLAGS="-DXyce_NOX_LOCA_SUPPORT $CPPFLAGS"

dnl Check if NOX is using NOX::size_type or int as the return value for Vector::length()
USING_NOX_SIZETYPE=no
XYCE_TEST_USING_NOXSIZETYPE
if test "$USING_NOX_SIZETYPE" != "no"; then
  CPPFLAGS="-DXyce_NOX_SIZETYPE $CPPFLAGS"
fi

dnl Check if NOX generic solver interface has getSolverStatistics() method
USING_NOX_SOLVERSTATS=no
XYCE_TEST_USING_NOXSOLVERSTATS
if test "$USING_NOX_SOLVERSTATS" != "no"; then
  CPPFLAGS="-DXyce_NOX_SOLVERSTATS $CPPFLAGS"
fi


dnl check for the dakota library if the user wants to link to it
if test "x$USE_DAKOTA" = "xyes"; then

  echo "checking for Dakota libraries."

  XYCE_LIBRARY_OPTION(boost_serialization)
  XYCE_LIBRARY_OPTION(boost_signals)
  XYCE_LIBRARY_OPTION(boost_system)
  XYCE_LIBRARY_OPTION(boost_filesystem)
  XYCE_LIBRARY_OPTION(boost_regex)

  XYCE_DAKOTA_LIBRARY_OPTION(Dakota, USE_DAKOTA_MAKEFILE, yes)
  if test "x$USE_DAKOTA_MAKEFILE" = "xno"; then

    AC_MSG_ERROR([Sorry, due to circular dependencies in Dakota libraries, we cannot use the normal fall-back technique when Makefile.export.Dakota is not found.  Please set DAKOTADIR to the path to your Dakota install and rerun configure.])
  fi
fi

dnl *****************************************************************
dnl Now include the include directories in CPPFLAGS
dnl We don't need to check directory existence, because we cleared the
dnl variables above if the directory didn't exist.
RAW_CPPFLAGS=${CPPFLAGS}
AC_SUBST([RAW_CPPFLAGS])

CPPFLAGS="${CPPFLAGS} ${INCDIRS}"

dnl Lastly, add the current directory to CPPFLAGS
CPPFLAGS="-I. $CPPFLAGS"

dnl *********************************************************************
dnl Add all the extra things we offered the user the option of including
dnl --------- Set CPPFLAGS and bindir ---------

CPPFLAGS="${CPPFLAGS} ${XTRAFLAGS} ${XTRAINC}"

dnl --------- Set and Check Libs ---------
dnl add the ones we requested on the command line
LIBS="${LIBS} ${XTRALIBS}"

#dnl put original LDFLAGS back and keep -Ldirs on left
#if test "$Xyce_ARCH" = "solaris"; then
#    LDFLAGS="$TMP_LDFLAGS $LDFLAGS"
#fi

dnl Need to define HAVE_CONFIG_H because AC_CHECK_HEADER doesn't include the
dnl CXXFLAGS on the command line when it tries to compile the test program.
dnl Oddly, AC_CHECK_HEADERS does.
AC_CHECK_HEADER(Epetra_MultiVector.h,,[AC_MSG_ERROR([Sorry, Petra headers not found])],)
AC_CHECK_HEADER(Ifpack_CrsRiluk.h,,[AC_MSG_ERROR([Sorry, Ifpack headers not found])],)

if test "x$USE_MPI" = "xyes"; then
if test "x$USE_ZOLTAN" = "xyes"; then
if test "x$USE_HMETIS" = "xno"; then
  AC_CHECK_HEADER(Zoltan_QueryObject.h,,[AC_MSG_ERROR([Sorry, Zoltan headers not found])])
fi
fi
fi
AC_CHECK_HEADER(Amesos_Klu.h,,[AC_MSG_ERROR([Sorry, KLU headers not found])])


if test "x$SUPERLU" = "xyes"; then
  AC_CHECK_HEADER(Amesos_Superlu.h,,[AC_MSG_ERROR([Trilinos was not built with SuperLU enabled, please run configure with --disable-superlu])])
fi

AC_CHECK_HEADER(NOX_Abstract_Vector.H,,[AC_MSG_ERROR([Sorry, NOX headers not found])])
AC_CHECK_HEADER(LOCA_Parameter_Vector.H,,[AC_MSG_ERROR([Sorry, LOCA headers not found])])
AC_CHECK_HEADER([Sacado.hpp], [], [AC_MSG_ERROR([Sorry, Sacado headers not found])])
AC_LANG_POP(C++)

dnl Checks for typedefs, structures, and compiler characteristics.

#check if we can use floating point exceptions on linux
##AC_CXX_HAVE_FENV

dnl Checks for library functions.
AC_CHECK_FUNCS(strcasecmp)

dnl Checks for library functions.
AC_CHECK_FUNCS(getopt)

AC_CHECK_FUNCS(drand48)

dnl Checks for library function erf 
AC_CHECK_FUNCS(erf)

dnl Checks for library function erfc
AC_CHECK_FUNCS(erfc)

# check if FFT library support is enabled. If so, check for Intel FFT or FFTW.
if test "X$USE_FFT" = "Xyes"; then

AC_DEFINE([Xyce_USE_FFT],[1], [Set to use either FFTW or Intel's FFT library.])
#
# check if the Intel FFT functions are available 
# this should only pass if 
#  1) the right headers are around (mkl_dfti.h)
#  2) the right math libraries have been specified
if test "X$USE_IntelFFT" = "Xyes" ; then
  AC_CHECK_HEADER([mkl_dfti.h], [AC_SEARCH_LIBS( [DftiCreateDescriptor], [], [AC_DEFINE([Xyce_USE_INTEL_FFT],[1], [Set to use Intel's FFT library])],[AC_MSG_WARN([Could not find Intel FFT Library.  Disabling Intel FFT support]); ac_cv_header_mkl_dfti_h=no])])
else
  ac_cv_header_mkl_dfti_h=no
fi

# check if the FFTW functions are available 
#
# check if the FFTW library header is around (version 3+), and if it is,
# check for the library
if test "X$USE_FFTW" = "Xyes"; then
  AC_CHECK_HEADER([fftw3.h], [AC_SEARCH_LIBS( [fftw_plan_dft_1d], [fftw3], [AC_DEFINE([Xyce_USE_FFTW],[1], [Set to use FFTW library])],[AC_MSG_WARN([Could not find FFTW Library.  Disabling FFTW support]); ac_cv_header_fftw3_h=no])])
else
 ac_cv_header_fftw3_h=no
fi  

# ERROR, neither FFT library was found and FFT was enabled.
if test "X$ac_cv_header_mkl_dfti_h" != "Xyes" -a "X$ac_cv_header_fftw3_h" != "Xyes"; then
  AC_MSG_ERROR(
                  [Configure has not found an FFT library (FFTW or Intel FFT) in your path.
                  An FFT library is required for frequency-domain analysis capabilities.
                  Reconfigure Xyce with either FFTW or Intel FFT or disable FFT support
                  with --disable-fft.])
fi

else
 ac_cv_header_fftw3_h=no
 ac_cv_header_mkl_dfti_h=no
fi

AM_CONDITIONAL(BUILD_INTEL_FFT_SUPPORT, test "X$ac_cv_header_mkl_dfti_h" = "Xyes")
AM_CONDITIONAL(BUILD_FFTW_SUPPORT, test "X$ac_cv_header_fftw3_h" = "Xyes")

# check if the HDF5 functions are available 
#
# check if the HDF5 library header is around (version 1.8+), and if it is,
# check for the library (Note. HDF5 can depend on libsz.a and libz.a,
# which are available on their web site.  The libz.a used by hdf5 IS NOT the 
# same as libzc.a used for compressing files.)

if test "X$USE_HDF5" = "Xyes"; then
  AC_SEARCH_LIBS( [inflate], [z] )
  AC_SEARCH_LIBS( [SZ_Decompress], [sz])
  if test "$USE_MPI" = yes; then
    # look for an HDF5 function that is only defined if HDF is also
    # built with MPI support.  Currently Xyce in parallel needs 
    # Parallel HDF5 to do IO in parallel.
    AC_CHECK_HEADER( [hdf5.h],
                     [AC_SEARCH_LIBS( [H5Pset_fapl_mpio], [hdf5],
                                      [AC_DEFINE([Xyce_USE_HDF5],[1], [Set to use HDF5 library])],
                                      [AC_MSG_WARN([Could not find HDF5 Library with MPI support.  Disabling parallel HDF5 support])
                                       ac_cv_header_hdf5_h=no
                                      ])]
                   )
  else
    # configuring Xyce to build in serial, so only look for serial functions 
    # to check HDF5 availabilit 
    AC_CHECK_HEADER( [hdf5.h],
                     [AC_SEARCH_LIBS( [H5Fclose], [hdf5],
                                      [AC_DEFINE([Xyce_USE_HDF5],[1], [Set to use HDF5 library])],
                                      [AC_MSG_WARN([Could not find HDF5 Library.  Disabling HDF5 support])
                                       ac_cv_header_hdf5_h=no
                                      ])]
                   )
  fi      
else
 ac_cv_header_hdf5_h=no
fi  
AM_CONDITIONAL(BUILD_HDF5_SUPPORT, test "X$ac_cv_header_hdf5_h" = "Xyes")

dnl *********************************************************************
dnl cURL support
dnl if requested, check for curl.h header, if exists check for library, if 
dnl both exist, enable cURL.  Disable if any of the the conditions fails

if test "X$USE_CURL" = "Xyes"; then
  AC_CHECK_HEADER([curl/curl.h],
                  [AC_SEARCH_LIBS([curl_easy_init],[curl],
                                  [AC_DEFINE([Xyce_USE_CURL],[1],[Set to use cURL library])],
                                  [AC_MSG_WARN([Could not find cURL library.  Disabling cURL support])
                                  ac_cv_header_curl_curl_h=no
                                  ])])
else
  ac_cv_header_curl_curl_h=no
fi

XYCE_CHECK_IOTA

dnl Set the ATHENA requirements.
AS_IF([test x$enable_athena = xyes],
      [AS_IF([test x$RAD_MODELS = xyes],
             [AS_IF([test -d "${srcdir}/src/DeviceModelPKG/SandiaModels/ATHENA"],
                    [AC_DEFINE([Xyce_ATHENA], [1], [Build the ATHENA device.])],
                    [AC_MSG_ERROR([Your SandiaModels directory is too old.  ATHENA is not present.])])],
             [enable_athena=no])
      ],[enable_athena=no])



AM_CONDITIONAL([BUILD_ATHENA], [test x$enable_athena = xyes])

dnl Check for the ARAEnv library
dnl do it here as it requires the same boost libraries that ATHENA does.
dnl should be smarter about managing boost dependencies
AC_LANG_PUSH(C++)
AC_CHECK_LIB([ARAEnv], [main], [], [], [])
if test "x$ac_cv_lib_ARAEnv_main" = xyes; then
     AC_CHECK_HEADER([SRenvironment.h], [AC_DEFINE([HAVE_SRENVIRONMENT_H], [1],[Set to have SRenvironment.h])])
fi
AC_LANG_POP(C++)

dnl LASTLY:  add -nofor_main to FFLAGS for DEC.  Can't do it earlier than this
if test "$USING_DEC_CXX" != "no"; then
  FFLAGS="$FFLAGS -nofor_main"
fi

# Only after ALL libraries are added to the link line, check if it'll work:
echo $LIBS

AC_MSG_CHECKING(whether libraries (${LIBS}) will link)
AC_LANG_PUSH(C++)
AC_LINK_IFELSE(
             [AC_LANG_PROGRAM([[]],[[]])],
[AC_MSG_RESULT([yes])],
[
AC_MSG_RESULT([no])
AC_MSG_ERROR([FATAL: cannot link with ${LIBS}... check paths, fix with --with-package options.   See config.log for failure message.])
]
)

dnl necessary Xyce substitution so all the makefiles get it.
dnl This makes @Xyce_INCS@ be replaced with the full string of -I's needed
dnl for all packages to compile without problem.
dnl It is used by setting "INCLUDES = @Xyce_INCS@" in the Makefile.am for
dnl a package.
default_INCS="$default_INCS -I\$(top_srcdir)/src/AnalysisPKG -I\$(top_srcdir)/src/CircuitPKG -I\$(top_srcdir)/src/IOInterfacePKG -I\$(top_srcdir)/src/IOInterfacePKG/Output -I\$(top_srcdir)/src/DeviceModelPKG/Core -I\$(top_srcdir)/src/DeviceModelPKG/OpenModels -I\$(top_srcdir)/src/DeviceModelPKG/TCADModels -I\$(top_srcdir)/src/DeviceModelPKG/IBISModels -I\$(top_srcdir)/src/DeviceModelPKG/EXTSC -I\$(top_srcdir)/src/DeviceModelPKG/NeuronModels -I\$(top_srcdir)/src/DeviceModelPKG/ADMS -I\$(top_srcdir)/src/ErrorHandlingPKG -I\$(top_srcdir)/src/LinearAlgebraServicesPKG -I\$(top_srcdir)/src/LoaderServicesPKG -I\$(top_srcdir)/src/NonlinearSolverPKG -I\$(top_srcdir)/src/ParallelDistPKG -I\$(top_srcdir)/src/TimeIntegrationPKG -I\$(top_srcdir)/src/TopoManagerPKG -I\$(top_srcdir)/src/UtilityPKG -I\$(top_srcdir)/src/MultiTimePDEPKG -I\$(top_srcdir)/src/DakotaLinkPKG -I\$(top_srcdir)/utils/XyceCInterface -I\$(top_builddir)/src"
if test "x$RAD_MODELS" = "xyes"; then
  default_INCS="$default_INCS -I\$(top_srcdir)/src/DeviceModelPKG/SandiaModels"
fi
if test "x$NONFREE_MODELS" = "xyes"; then
  default_INCS="$default_INCS -I\$(top_srcdir)/src/DeviceModelPKG/Xyce_NonFree"
fi

if test "x$USE_KSPARSE" = "xyes"; then
  default_INCS="$default_INCS -I\$(top_srcdir)/src/LinearAlgebraServicesPKG/ksparse"
fi
if test "x$MODSPEC_MODELS" = "xyes"; then
   default_INCS="$default_INCS -I\$(top_srcdir)/src/DeviceModelPKG/BERKELEY/include -I\$(top_srcdir)/src/DeviceModelPKG/BERKELEY/Xyce_ModSpec_Interface/include"
fi
AC_SUBST(Xyce_INCS,[$default_INCS])

XYCE_SRCDIR="${srcdir}"
AC_SUBST(XYCE_SRCDIR)

XYCELIBS="${LIBS}"
LIBS=" "
AC_SUBST(XYCELIBS)
AC_SUBST(LIBTOOLFLAGS)

# shouldn't need to export this, but Dakota requires this for one of its object files
AC_SUBST(DAKOTA_ARCHDIR)

# Ugly stuff:  We must use different CXXFLAGS for Intel compilers of the wrong
# version in the ADMS directory:
AC_SUBST(ADMS_CXXFLAGS,[$ADMS_CXXFLAGS])

#############################################################################
# CPACK STUFF
# FINALLY:  Now that we know everything about this build, let's construct
# a set of variables that we can use to generate a CPack config file, so
# that we can use CPack to make a native installer.
# Sorry, but "make dist" just ain't good enough for us, and we do want
# the ability to create native installers.

# Make sure we at least have a generator defined.  TGZ is just a tarball.  We
# will certainly want to override this on the configure command line instead.
# our "reconfigure.cfg" file has settings for our supported platforms, so 
# this variable is always set on the command line in the build/test cycle.
if test "x$CPACK_GENERATOR" = "x" 
then
   CPACK_GENERATOR="TGZ"
fi
AC_SUBST([CPACK_GENERATOR])

# Now, let's decode the VERSION string that was defined in AC_INIT above.
DEVEL_OR_RELEASE=`echo $VERSION | awk -F: '{print $1}'`
if test $DEVEL_OR_RELEASE = "D" 
then
   DEVEL_OR_RELEASE="DEVELOPMENT"
else
   DEVEL_OR_RELEASE="Release"
fi

DEVEL_PHASE=`echo $VERSION | awk -F: '{print $2}'`
if test "x$DEVEL_PHASE" != "x" 
then
   DEVEL_PHASE="-$DEVEL_PHASE"
fi

VERSION_NUMBER=`echo $VERSION | awk -F: '{print $3}'`
Xyce_VERSION_MAJOR=`echo $VERSION_NUMBER | awk -F. '{print $1}'`
Xyce_VERSION_MINOR=`echo $VERSION_NUMBER | awk -F. '{print $2}'`
if test "x$Xyce_VERSION_MINOR" = "x"
then
  Xyce_VERSION_MINOR=0
fi
Xyce_VERSION_PATCH=`echo $VERSION_NUMBER | awk -F. '{print $3}'`
if test "x$Xyce_VERSION_PATCH" = "x"
then
  Xyce_VERSION_PATCH=0
fi

Xyce_VERSION_LONG="${Xyce_VERSION_MAJOR}.${Xyce_VERSION_MINOR}.${Xyce_VERSION_PATCH}"
if test "x$Xyce_VERSION_EXTRA" != "x"
then
  Xyce_VERSION_LONG="${Xyce_VERSION_MAJOR}.${Xyce_VERSION_MINOR}.${Xyce_VERSION_PATCH}-${Xyce_VERSION_EXTRA}"
fi

Xyce_VERSION="${Xyce_VERSION_MAJOR}.${Xyce_VERSION_MINOR}.${Xyce_VERSION_PATCH}"

AC_SUBST([Xyce_VERSION])
AC_SUBST([Xyce_VERSION_MAJOR])
AC_SUBST([Xyce_VERSION_MINOR])
AC_SUBST([Xyce_VERSION_PATCH])

CPACK_PACKAGE_NAME="Xyce-$DEVEL_OR_RELEASE-$Xyce_VERSION_LONG$DEVEL_PHASE"
AC_SUBST([CPACK_PACKAGE_NAME])
if test "x$CPACK_PACKAGE_FILE_NAME" = "x"
then
  CPACK_PACKAGE_FILE_NAME="Install_Xyce_$Xyce_VERSION_MAJOR.$Xyce_VERSION_MINOR.$Xyce_VERSION_PATCH"
fi
AC_SUBST([CPACK_PACKAGE_FILE_NAME])

CPACK_RPM_PACKAGE_NAME=$CPACK_PACKAGE_NAME
AC_SUBST([CPACK_RPM_PACKAGE_NAME])

# Set CPACK_SYSTEM_NAME and default packaging prefix
# This is a little annoying, but if you package on OSX with anything other than
# / for the install prefix
# you get something like /usr/local/Xyce.../usr/local/bin, which is just
# nuts.
# NOTE:  OS X, starting with El Capitan, refuses to respect the
# DYLD_LIBRARY_PATH environment variable in subshells, so have to force an
# RPATH to get our release builds to work.
# NOTE:  We're also adding and RPATH to Linux builds, so that we can someday
# deprecate and do away with the worthless "runxyce" and "xmpirun" scripts.
#
CPACK_PACKAGING_INSTALL_PREFIX="/usr/local"
case "$target" in
  *-freebsd*) 
    CPACK_SYSTEM_NAME="FreeBSD"
    ;;
  *-linux*)
    CPACK_SYSTEM_NAME="Linux"
    case "$CXX" in
    icpc|*-icc-*mpic++|*-icc-*mpicxx)
        CXXFLAGS="'-Wl,-rpath,\$\$ORIGIN/../lib' $CXXFLAGS"
        ;;
    mpicxx|mpic++)
        # Handle case where we're relying on PATH and OMPI_CXX to use intel
        if test "x$OMPI_CXX" = "xicpc" 
        then
          CXXFLAGS="'-Wl,-rpath,\$\$ORIGIN/../lib' $CXXFLAGS"
        fi
        ;;
    esac
    ;;
  *-apple-darwin*)
    CPACK_SYSTEM_NAME="Darwin"
    CPACK_PACKAGING_INSTALL_PREFIX="/"
    case "$CXX" in
    icpc|*-icc-*mpic++|*-icc-*mpicxx)
        CXXFLAGS="-Wl,-rpath,/usr/local/$CPACK_PACKAGE_NAME/lib/ $CXXFLAGS"
        ;;
    esac
    ;;
  *)
    CPACK_SYSTEM_NAME=$target
    ;;
esac
AC_SUBST([CPACK_SYSTEM_NAME])
AC_SUBST([CPACK_PACKAGING_INSTALL_PREFIX])

# OK, with this one we're going to try to force generated RPMs for OPENMPI
# builds to pull in the openmpi executable package.  Will it work?
if test "x$USE_MPI" = "xyes"  -a "$CPACK_GENERATOR" = "RPM"
then
  CPACK_RPM_DEPENDS_SET_LINE="SET (CPACK_RPM_PACKAGE_REQUIRES \"openmpi >= 1.4\")"
else
  CPACK_RPM_DEPENDS_SET_LINE=""
fi
AC_SUBST([CPACK_RPM_DEPENDS_SET_LINE])

# Finally, figure out what license we're using.

XYCE_LICENSE=ECI
CPACK_RPM_PACKAGE_LICENSE="Export Controlled"
if test $RAD_MODELS = "no"
then
  if test $NONFREE_MODELS = "no"
  then
     XYCE_LICENSE=OS
     CPACK_RPM_PACKAGE_LICENSE="GPLv3"
  else
     XYCE_LICENSE=NonFree
     CPACK_RPM_PACKAGE_LICENSE="NDA Required"
  fi
fi
AC_SUBST([XYCE_LICENSE])
AC_SUBST([CPACK_RPM_PACKAGE_LICENSE])

#
# END CPACK STUFF
############################################################################

AC_CONFIG_FILES([Makefile 
 src/Makefile 
 src/AnalysisPKG/Makefile 
 src/CircuitPKG/Makefile
 src/IOInterfacePKG/Makefile
 src/IOInterfacePKG/Output/Makefile
 src/DakotaLinkPKG/Makefile 
 src/DeviceModelPKG/Makefile 
 src/DeviceModelPKG/Core/Makefile 
 src/DeviceModelPKG/OpenModels/Makefile 
 src/DeviceModelPKG/TCADModels/Makefile
 src/DeviceModelPKG/IBISModels/Makefile 
 src/DeviceModelPKG/NeuronModels/Makefile 
 src/DeviceModelPKG/ADMS/Makefile 
 src/DeviceModelPKG/EXTSC/Makefile 
 src/ErrorHandlingPKG/Makefile 
 src/LinearAlgebraServicesPKG/Makefile 
 src/LinearAlgebraServicesPKG/ksparse/Makefile 
 src/LoaderServicesPKG/Makefile 
 src/NonlinearSolverPKG/Makefile 
 src/ParallelDistPKG/Makefile 
 src/TimeIntegrationPKG/Makefile 
 src/TopoManagerPKG/Makefile 
 src/UtilityPKG/Makefile 
 src/MultiTimePDEPKG/Makefile 
 utils/Makefile
 utils/XyceCInterface/Makefile 
 src/test/Makefile
 src/test/XyceAsLibrary/Makefile
 src/test/DeviceInterface/Makefile
 src/test/RandomTest/Makefile
 src/test/FFTInterface/Makefile
 src/test/LinearAlgebraTest/Makefile
 src/test/XygraTestHarnesses/Makefile
 src/test/GenExtTestHarnesses/Makefile
 src/test/MPITest/Makefile
 src/test/TwoLevelNewton/Makefile
 src/test/VPITests/Makefile
 user_plugin/Makefile
 distribution/CPackConfig.cmake])

#NOTE:
# By putting the makefiles into a shell variable like this we are in fact
# deliberately ignoring the guidance in Automake documentation for a specific
# reason.
# 
# When automake encounters the AC_CONFIG_FILES line here, it cannot tell
# that it needs to create the Makefile.in that AC_CONFIG_FILES will need.
# Thus, it is necessary to tell automake on its command line that it needs
# to do that generation.
#
# Why?  Because these two directories are optional, and are not generally 
# distributed with Xyce --- they are not open source directories, and we can't
# ship them.  But at automake/bootstrap time, automake can't tell that these
# AC_CONFIG_FILES are conditionalized, and if the directories don't exist
# then it pukes with an error that it can't find required files.
#
# So the solution is to confuse automake into thinking it has nothing to do
# automatically, ever, and explicitly tell it what to do in bootstrap if
# the directories actually exist.

if test "x$RAD_MODELS" = "xyes"; then
  rm_makefile="src/DeviceModelPKG/SandiaModels/Makefile src/DeviceModelPKG/SandiaModels/ATHENA/Makefile"
  AC_CONFIG_FILES([$rm_makefile])
fi
if test "x$NONFREE_MODELS" = "xyes"; then
  nf_makefile=src/DeviceModelPKG/Xyce_NonFree/Makefile
  AC_CONFIG_FILES([$nf_makefile])
fi
if test "x$MODSPEC_MODELS" = "xyes"; then
  ms_makefile1=src/DeviceModelPKG/BERKELEY/Makefile 
  ms_makefile2=src/DeviceModelPKG/BERKELEY/Xyce_ModSpec_Interface/Makefile 
  AC_CONFIG_FILES([$ms_makefile1])
  AC_CONFIG_FILES([$ms_makefile2])
fi


AC_OUTPUT
