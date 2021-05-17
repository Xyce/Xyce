
AM_CPPFLAGS = @Xyce_INCS@

TEST_RANDOM_SOURCES = \
  testRandom.C

# standalone executable
check_PROGRAMS = testRandom
testRandom_SOURCES= $(TEST_RANDOM_SOURCES)
testRandom_LDADD = $(top_builddir)/src/libxyce.la
testRandom_LDFLAGS = -static $(AM_LDFLAGS)
 
