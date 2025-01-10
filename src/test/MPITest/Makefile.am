
AM_CPPFLAGS = @Xyce_INCS@

TEST_MPI_SOURCES = \
  testBUG967.c

# standalone executable
check_PROGRAMS = testMPI
testMPI_SOURCES = $(TEST_MPI_SOURCES)
testMPI_LDADD = 
testMPI_LDFLAGS = -static $(AM_LDFLAGS) 
 
