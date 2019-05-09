
# XYCE_LIBRARY_OPTION(LIBRARY[, INCLUDE-DIRECTORY[, FUNCTION[, ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]]]])
# Attempts to link with the library, and includes it in LIBS if its found.
# The default options mearly determine if the library can be linked against.
# If INCLUDE-DIRECTORY is specified, add it to the INCDIRS variable as "-I${ARCHDIR}/include/INCLUDE-DIRECTORY"
# Call FUNCTION in library test.  Defaults to "main".
# Run ACTION-IF-FOUND if successful.  Defaults to AC_CHECK_LIB's actions.
# Run ACTION-IF-NOT-FOUND if unsuccessful.  Defaults to calling AC_MSG_ERROR.

AC_DEFUN([XYCE_LIBRARY_OPTION],
[
dnl the library name, like "petra"
m4_pushdef([xyce_lib], [$1]) dnl

dnl the include directory name, like "epetra", if different from the package name
m4_ifval([$2],
[if test -d "${ARCHDIR}/include/$2"; then
  INCDIRS="${INCDIRS} -I${ARCHDIR}/include/$2"
fi])

m4_pushdef([func], [m4_default([$3], [main])]) dnl
AC_CHECK_LIB(xyce_lib, func, [$4],
  [m4_default([$5], [AC_MSG_ERROR([Unable to find required library, xyce_lib.])])]
)
m4_popdef([func]) dnl
m4_popdef([xyce_lib]) dnl
])

