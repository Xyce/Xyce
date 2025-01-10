
AM_CPPFLAGS = @Xyce_INCS@

TEST_FFT_INTERFACE_SOURCES = \
  testFFTInterface.C

# standalone executable
check_PROGRAMS = testFFTInterface
testFFTInterface_SOURCES = $(TEST_FFT_INTERFACE_SOURCES)
testFFTInterface_LDADD = $(top_builddir)/src/libxyce.la
testFFTInterface_LDFLAGS = -static $(AM_LDFLAGS)
 
