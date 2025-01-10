AC_DEFUN([XYCE_LIBRARY_INCLUDE],
[

# the library name, like "petra"
xyce_lib_directory_name=$1
if test -d $ARCHDIR/include/${xyce_lib_directory_name}; then
   INCDIRS="$INCDIRS -I$ARCHDIR/include/${xyce_lib_directory_name}"
fi

])

