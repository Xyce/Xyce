# NOX has changed the abstract solver interface to include a new method that returns
# solver statistics.  We are testing for the existence of the new struct that contains
# those statistics, since it did not exist before this change.
AC_DEFUN([XYCE_TEST_USING_NOXSOLVERSTATS],
[AC_CACHE_CHECK([whether the NOX abstract solver has getSolverStatistics()],
ac_cv_xyce_using_noxsolverstats,
[AC_LANG_PUSH(C++)
 AC_COMPILE_IFELSE(
[AC_LANG_PROGRAM([[#include "NOX_SolverStats.hpp"]],[[
]])],ac_cv_xyce_using_noxsolverstats=yes, ac_cv_xyce_using_noxsolverstats=no)
AC_LANG_POP(C++)
])
if test "$ac_cv_xyce_using_noxsolverstats" = yes; then
  USING_NOX_SOLVERSTATS=yes
fi
])
