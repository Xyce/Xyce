# NOX has changed the return type on the length() method of their abstract vector class.
# We need to detect for older versions of NOX, where the return type is still an int.
AC_DEFUN([XYCE_TEST_USING_NOXSIZETYPE],
[AC_CACHE_CHECK([whether the NOX abstract vector length is a NOX::size_type],
ac_cv_xyce_using_noxsizetype,
[AC_LANG_PUSH(C++)
 AC_COMPILE_IFELSE(
[AC_LANG_PROGRAM([[#include "NOX_Common.H"]],[[
NOX::size_type nox_integer;
]])],ac_cv_xyce_using_noxsizetype=yes, ac_cv_xyce_using_noxsizetype=no)
AC_LANG_POP(C++)
])
if test "$ac_cv_xyce_using_noxsizetype" = yes; then
  USING_NOX_SIZETYPE=yes
fi
])
