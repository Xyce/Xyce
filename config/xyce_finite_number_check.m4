dnl
dnl Check to make sure nan and inf detection is supported.
dnl
dnl Author: Roger Pawlowski
dnl
dnl
AC_DEFUN([AC_FINITE_NUMBER_CHECK],
[AC_CACHE_CHECK(whether the compiler supports isnan and isinf,
ac_cv_xyce_have_finite_number,
[AC_LANG_PUSH(C++)
 ac_save_LIBS="$LIBS"
 LIBS="$LIBS -lm"
AC_LINK_IFELSE(
[AC_LANG_PROGRAM([[
#ifndef _ALL_SOURCE
 #define _ALL_SOURCE
#endif
#ifndef _XOPEN_SOURCE
 #define _XOPEN_SOURCE
#endif
#ifndef _XOPEN_SOURCE_EXTENDED
 #define _XOPEN_SOURCE_EXTENDED 1
#endif
#include <cmath>
]],[[

double x = 1.0; 
std::isnan(x); std::isinf(x);
return 0;]])],
 ac_cv_xyce_have_finite_number=yes, ac_cv_xyce_have_finite_number=no)
 LIBS="$ac_save_LIBS"
 AC_LANG_POP(C++)
])
if test "$ac_cv_xyce_have_finite_number" = yes; then
  AC_DEFINE(HAVE_NAN_INF_SUPPORT,1,[define if the compiler supports isnan and isinf checks])
else
  echo "****************************************************"
  echo "** Warning: Your compiler doesn't support isnan() and "
  echo "** isinf().  We will supply a default checker but it "
  echo "** is *NOT* guaranteed to work on your platform!"
  echo "****************************************************"
fi

])

