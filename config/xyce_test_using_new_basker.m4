# Basker has changed the class name for the abstract solver to deconflict the
# templated Basker (used in block analysis methods) and ShyLU-Basker (multi-threaded
# nodal solver).  This allows them to be buildable in the same executable. 
AC_DEFUN([XYCE_TEST_USING_NEW_BASKER],
[AC_CACHE_CHECK([whether the Basker abstract solver in Amesos2 has the old interface],
ac_cv_xyce_old_basker,
[AC_LANG_PUSH(C++)
 AC_COMPILE_IFELSE(
[AC_LANG_PROGRAM([[#include "Amesos2_Basker.hpp"]],
[[Basker::Basker<int, double> basker_;]])],ac_cv_xyce_old_basker=yes, ac_cv_xyce_old_basker=no)
AC_LANG_POP(C++)
])
AC_CACHE_CHECK([whether the Basker abstract solver in Amesos2 has the new interface],
ac_cv_xyce_new_basker,
[AC_LANG_PUSH(C++)
 AC_COMPILE_IFELSE(
[AC_LANG_PROGRAM([[#include "Amesos2_Basker.hpp"]],
[[BaskerClassicNS::BaskerClassic<int, double> basker_;]])],ac_cv_xyce_new_basker=yes, ac_cv_xyce_new_basker=no)
AC_LANG_POP(C++)
])
if test "$ac_cv_xyce_new_basker" = yes; then
  USING_NEW_BASKER=yes
fi
if test "$ac_cv_xyce_old_basker" = yes; then
  USING_OLD_BASKER=yes
fi
])
