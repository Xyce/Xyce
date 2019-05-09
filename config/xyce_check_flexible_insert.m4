AC_DEFUN([XYCE_CHECK_FLEXIBLE_INSERT],
[AC_CACHE_CHECK(whether your STL supports flexible insertions,
ac_cv_xyce_flexible_insert,
[AC_LANG_PUSH(C++)
 AC_COMPILE_IFELSE(
[AC_LANG_PROGRAM([[#include <vector>
#include <list>
]],
[[
  std::list<int> ipList(1,0);
  std::vector<int> ipVec;
  
  ipVec.insert(ipVec.begin(),ipList.begin(),ipList.end());

]])], ac_cv_xyce_flexible_insert=yes , ac_cv_xyce_flexible_insert=no)
AC_LANG_POP(C++)
])
if test "$ac_cv_xyce_flexible_insert" = yes; then
  AC_DEFINE(HAVE_FLEXIBLE_INSERT,1,[define if the Standard Template Library has flexible insert])
fi
])
