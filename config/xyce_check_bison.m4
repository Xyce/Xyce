#Check appropriate version of bison
AC_DEFUN([XYCE_CHECK_BISON],[

succeeded=no

AC_CHECK_PROG([BISON],[bison],[bison],[no])

if test "$BISON" != "no"
then
  BISON_VERSION=`$BISON --version | cut -d ' ' -f 4 | head -1`
  BISON_MAJOR_VERSION=`echo $BISON_VERSION | cut -d '.' -f 1`
  BISON_MINOR_VERSION=`echo $BISON_VERSION | cut -d '.' -f 2`

  m4_foreach([AC_xyce_bison_ver],[$1],
  [
  if test "$succeeded" != "yes"
  then
    DESIRED_MAJOR_VERSION=`echo AC_xyce_bison_ver|cut -d "." -f 1`
    DESIRED_MINOR_VERSION=`echo AC_xyce_bison_ver|cut -d "." -f 2`
    AC_MSG_CHECKING(for bison version AC_xyce_bison_ver)
  
    VERSION_CHECK=`expr $BISON_MAJOR_VERSION \= $DESIRED_MAJOR_VERSION \& $BISON_MINOR_VERSION \= $DESIRED_MINOR_VERSION`
    if test "$VERSION_CHECK" = "1" ; then
       AC_MSG_RESULT(yes)
       succeeded=yes
    else
       AC_MSG_RESULT(no)
    fi
   fi
   ]
  )
  fi

if test $succeeded = yes; then
     ifelse([$2], , :, [$2])
else
     BISON=no
     ifelse([$3], , AC_MSG_ERROR([Usable version of Bison not found.]), [$3])
fi

])
