
AM_CPPFLAGS = @Xyce_INCS@


# standalone testHarness1685 executable
check_PROGRAMS = testHarness1685
 
BUG1685TESTSOURCES = \
  testHarness1685.C 

# standalone testHarness1685  executable
testHarness1685_SOURCES = $(BUG1685TESTSOURCES)
testHarness1685_LDADD = $(top_builddir)/src/libxyce.la
testHarness1685_LDFLAGS = -static $(AM_LDFLAGS)
 
