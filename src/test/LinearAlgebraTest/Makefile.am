
AM_CPPFLAGS = @Xyce_INCS@

# needed for Dakota 4.x not 5.0
if DAKOTA_OBJ_NEEDED 
  DAKOTA_OBJS = 
endif

TEST_LINALG_SOURCES = \
  testBlockLinearSystems.C

# standalone executable
check_PROGRAMS = testBlockLinearSystems
testBlockLinearSystems_SOURCES = $(TEST_LINALG_SOURCES)
testBlockLinearSystems_LDADD = $(top_builddir)/src/libxyce.la
testBlockLinearSystems_LDFLAGS = -static $(AM_LDFLAGS) $(DAKOTA_OBJS)
 
