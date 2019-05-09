# we are seeing some issues with Intel compilers below version 11 and 
# full optimization of Verilog-generated codes in the ADMS directory.
# So let's check whether we're using Intel <12
AC_DEFUN([XYCE_TEST_USING_ICCLT11_CXX],
[AC_CACHE_CHECK([whether your C++ compiler is the Intel compiler version less than 12],
ac_cv_cxx_xyce_using_icclt11_cxx,
[AC_LANG_PUSH(C++)
 AC_COMPILE_IFELSE(
[AC_LANG_PROGRAM([[]],[[
#ifdef __ICC
#if __ICC < 1200
#define A 1
#else
error intel but later than version 12
#endif
#else
error not the intel compiler
#endif
int i = 1;
]])],ac_cv_cxx_xyce_using_icclt11_cxx=yes, ac_cv_cxx_xyce_using_icclt11_cxx=no)
AC_LANG_POP(C++)
])
if test "$ac_cv_cxx_xyce_using_icclt11_cxx" = yes; then
  USING_ICCLT11_CXX=yes
fi
])
