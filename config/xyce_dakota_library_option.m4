
# XYCE_DAKOTA_LIBRARY_OPTION(LIBRARY, USE_EXPORT_MAKEFILE, no)
# Attempts to link with the library, and includes it in LIBS if its found.
# The default options mearly determine if the library can be linked against.
# If INCLUDE-DIRECTORY is specified, add it to the INCDIRS variable as "-I${DAKOTA_ARCHDIR}/include/INCLUDE-DIRECTORY"
# Call FUNCTION in library test.  Defaults to "main".
# Run ACTION-IF-FOUND if successful.  Defaults to AC_CHECK_LIB's actions.
# Run ACTION-IF-NOT-FOUND if unsuccessful.  Defaults to calling AC_MSG_ERROR.

AC_DEFUN([XYCE_DAKOTA_LIBRARY_OPTION],
[
# By default we check for the Makefile.export, then check for the library if it doesn't exist.
dnl the library name, like "petra"
m4_pushdef([dak_lib], m4_tolower([$1])) dnl
m4_pushdef([func], [main]) dnl
m4_pushdef([only_makefile], [$3]) dnl

# If we can find the Makefile.export for Dakota, then add all the libraries to LIBS.
# Otherwise, we're doing this the old fashioned way with a static list of dependencies.
if test -e "${DAKOTA_ARCHDIR}/include/Makefile.export.$1"; then
  m4_ifval($2, [$2=yes])
  ALL_LIBS=`$GREP $1_LIBRARIES ${DAKOTA_ARCHDIR}/include/Makefile.export.$1 | sed -e 's/$1_LIBRARIES//' -e 's/=//'`
#  TPL_LIBS=`$GREP $1_TPL_LIBRARIES ${DAKOTA_ARCHDIR}/include/Makefile.export.$1 | sed -e 's/$1_TPL_LIBRARIES//' -e 's/=//'`
  AC_MSG_NOTICE([Found library and dependencies for $1 in Makefile.export.$1])
  LIBS="${ALL_LIBS} ${TPL_LIBS} ${LIBS}"
else
  m4_ifval($2, [$2=no])
  AC_MSG_WARN([Unable to find export file: Makefile.export.$1.  Using predefined library dependencies.])
  if test "only_makefile[]" != "yes"; then
    AC_CHECK_LIB(dak_lib, func, [], [AC_MSG_ERROR([Unable to find required library, dak_lib.])]) 
  fi
fi

m4_popdef([only_makefile],[func],[dak_lib]) dnl
])

