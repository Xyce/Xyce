dnl XYCE_DEBUG_OPTION(mpi,no,[parallel build with MPI],Xyce_PARALLEL_MPI,USE_MPI)
dnl create a "--enable-mpi" configure option with help text
dnl   "parallel build with MPI".
dnl using the configure option will add the Xyce_PARALLEL_MPI symbol to the
dnl   Xyce_config.h file, and define the "USE_MPI" variable.
AC_DEFUN([XYCE_DEBUG_OPTION],
[dnl
dnl get the arguments
m4_pushdef([xyce_name], [$1]) dnl
m4_pushdef([xyce_name_upcase], [m4_toupper(xyce_name)]) dnl
dnl
dnl default
m4_pushdef([xyce_default], [$2]) dnl
dnl
dnl symbol
m4_pushdef([xyce_symbol], m4_default($4, [Xyce_]xyce_name_upcase))
dnl
AC_ARG_ENABLE([$1], [AS_HELP_STRING([--enable-$1], [enable $3.])], [],
              [enable_[]xyce_name[]=xyce_default; enableval=$enable_[]xyce_name[]])
if test "x$enable_[]xyce_name[]" != "xno"; then
   AC_DEFINE(xyce_symbol,[1],[Set to use $1.])
fi
dnl
m4_ifval($5, [$5=$enable_[]xyce_name[]])
dnl
m4_popdef([xyce_name], [xyce_name_upcase], [xyce_default]) dnl
m4_popdef([xyce_symbol]) dnl
])
