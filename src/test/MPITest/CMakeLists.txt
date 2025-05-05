
if( Xyce_PARALLEL_MPI )
  add_executable( testMPI testBUG967.c )
  target_link_libraries( testMPI XyceLib )
  if( BUILD_TESTING AND Xyce_PARALLEL_MPI )
    add_test( NAME testMPI COMMAND mpiexec -np 2 testMPI )
    set_tests_properties( testMPI PROPERTIES LABELS "parallel;unittest" FAIL_REGULAR_EXPRESSION "BAD BAD BAD")
  endif()
endif()

