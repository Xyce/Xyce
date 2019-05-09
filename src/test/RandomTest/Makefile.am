
AM_CPPFLAGS = @Xyce_INCS@

# needed for Dakota 4.x not 5.0
if DAKOTA_OBJ_NEEDED 
  DAKOTA_OBJS = 
endif

TEST_RANDOM_SOURCES = \
  testRandom.C

# standalone executable
check_PROGRAMS = testRandom
testRandom_SOURCES= $(TEST_RANDOM_SOURCES)
testRandom_LDADD = $(top_builddir)/src/libxyce.la
testRandom_LDFLAGS = -static $(AM_LDFLAGS) $(DAKOTA_OBJS)
 
