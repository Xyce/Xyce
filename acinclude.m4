##########Does the stuff I wrote for Xyce work?
######################################################################
# These macros from the GNU Autoconf Macro Archive at 
# http://www.gnu.org/software/ac-archive/
######################################################################
AC_DEFUN([AC_CXX_NAMESPACES],
[AC_CACHE_CHECK(whether the compiler implements namespaces,
ac_cv_cxx_namespaces,
[AC_LANG_PUSH(C++)
 AC_COMPILE_IFELSE(
[AC_LANG_PROGRAM([[namespace Outer { namespace Inner { int i = 0; }}]],
                [[using namespace Outer::Inner; return i;]])],
 ac_cv_cxx_namespaces=yes, ac_cv_cxx_namespaces=no)
 AC_LANG_POP(C++)
])
if test "$ac_cv_cxx_namespaces" = yes; then
  AC_DEFINE(HAVE_NAMESPACES,,[define if the compiler implements namespaces])
fi
])

AC_DEFUN([AC_CXX_HAVE_NUMERIC_LIMITS],
[AC_CACHE_CHECK(whether the compiler has numeric_limits<T>,
ac_cv_cxx_have_numeric_limits,
[AC_REQUIRE([AC_CXX_NAMESPACES])
 AC_LANG_PUSH(C++)
 AC_COMPILE_IFELSE(
[AC_LANG_PROGRAM([[#include <limits>
]],[[double e = std::numeric_limits<double>::epsilon(); return 0;]])],
 ac_cv_cxx_have_numeric_limits=yes, ac_cv_cxx_have_numeric_limits=no)
 AC_LANG_POP(C++)
])
if test "$ac_cv_cxx_have_numeric_limits" = yes; then
  AC_DEFINE(HAVE_NUMERIC_LIMITS,1,[define if the compiler has numeric_limits<T>])
fi
])

########################################################################
# Exception handling for floating point exceptions can be enabled on
# linux.  As this is a common development and deployment platform, such
# handling on this platform should prevent this type of error from
# creeping into the code.
########################################################################
AC_DEFUN([AC_CXX_HAVE_FENV],
[AC_CACHE_CHECK(whether the compiler has fenv.h,
ac_cv_cxx_have_fenv_h,
[AC_REQUIRE([AC_CXX_NAMESPACES])
 AC_LANG_PUSH(C++)
 AC_COMPILE_IFELSE(
[AC_LANG_PROGRAM([[#include <fenv.h>
]],[[feenableexcept(FE_DIVBYZERO); return 0;]])],
 ac_cv_cxx_have_fenv_h=yes, ac_cv_cxx_have_fenv_h=no)
 AC_LANG_POP(C++)
])
if test "$ac_cv_cxx_have_fenv_h" = yes; then
  AC_DEFINE(HAVE_LINUX_EXCEPTIONS, 1,[define if the compiler has fenv include])
fi
])

AC_DEFUN([AC_CXX_HAVE_STL],
[AC_CACHE_CHECK(whether the compiler supports Standard Template Library,
ac_cv_cxx_have_stl,
[AC_REQUIRE([AC_CXX_NAMESPACES])
 AC_LANG_PUSH(C++)
 AC_COMPILE_IFELSE(
[AC_LANG_PROGRAM([[#include <list>
#include <deque>
]],[[std::list<int> x; x.push_back(5);
std::list<int>::iterator iter = x.begin(); if (iter != x.end()) ++iter; return 0;]])],
 ac_cv_cxx_have_stl=yes, ac_cv_cxx_have_stl=no)
 AC_LANG_POP(C++)
])
if test "$ac_cv_cxx_have_stl" = yes; then
  AC_DEFINE(HAVE_STL,,[define if the compiler supports Standard Template Library])
fi
])

