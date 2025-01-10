
AM_CPPFLAGS = @Xyce_INCS@

TEST_LINALG_SOURCES = \
  testBlockLinearSystems.C

# standalone executable
check_PROGRAMS = testBlockLinearSystems
testBlockLinearSystems_SOURCES = $(TEST_LINALG_SOURCES)
testBlockLinearSystems_LDADD = $(top_builddir)/src/libxyce.la
testBlockLinearSystems_LDFLAGS = -static $(AM_LDFLAGS)
 
