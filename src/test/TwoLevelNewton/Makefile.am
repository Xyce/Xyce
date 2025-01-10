
AM_CPPFLAGS = @Xyce_INCS@


# standalone XyceLibTest executable
check_PROGRAMS = twoLevelNewtonLinearDC twoLevelNewtonLinearTran
#twoLevelNewton

#twoLevelNewton_SOURCES = twoLevelNewton.C
#twoLevelNewton_LDADD = $(top_builddir)/src/libxyce.la
#twoLevelNewton_LDFLAGS = -static $(AM_LDFLAGS)

twoLevelNewtonLinearDC_SOURCES = twoLevelNewtonLinearDC.C baseNewton.h twoLevelNewtonLinear.h
twoLevelNewtonLinearDC_LDADD = $(top_builddir)/src/libxyce.la
twoLevelNewtonLinearDC_LDFLAGS = -static $(AM_LDFLAGS)

twoLevelNewtonLinearTran_SOURCES = twoLevelNewtonLinearTran.C baseNewton.h twoLevelNewtonLinear.h
twoLevelNewtonLinearTran_LDADD = $(top_builddir)/src/libxyce.la
twoLevelNewtonLinearTran_LDFLAGS = -static $(AM_LDFLAGS)

