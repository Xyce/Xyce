dnl @synopsis XYCE_ARG_WITH_Y12MLIB
dnl
dnl Test for --with-y12mlib="name".
dnl 
dnl Prepends the specified name to the list of files to check for BLAS
dnl routines.  
dnl
dnl @author Robert Hoekstra <rjhoeks@sandia.gov>
dnl
AC_DEFUN([XYCE_ARG_WITH_Y12MLIB],
[
AC_ARG_WITH(y12mlib,
AC_HELP_STRING([--with-y12mlib], 
[name of library containing Y12M: will search lib directories for
-lname]),
[
USE_Y12MLIB=yes
NEWY12MLIB=${withval}
]
)

   Y12MLIBS="y12m"

if test "X${USE_Y12MLIB}" = "Xyes"; then

   Y12MLIBS="${NEWY12MLIB} ${Y12MLIBS}"

fi
])

