
AM_CPPFLAGS = @Xyce_INCS@

# needed for Dakota 4.x not 5.0
if DAKOTA_OBJ_NEEDED 
  DAKOTA_OBJS = 
endif

TEST_FFT_INTERFACE_SOURCES = \
  testFFTInterface.C

# standalone executable
check_PROGRAMS = testFFTInterface
testFFTInterface_SOURCES = $(TEST_FFT_INTERFACE_SOURCES)
testFFTInterface_LDADD = $(top_builddir)/src/libxyce.la
testFFTInterface_LDFLAGS = -static $(AM_LDFLAGS) $(DAKOTA_OBJS)
 
