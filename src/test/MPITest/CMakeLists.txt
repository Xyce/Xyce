
if( BUILD_TESTING AND Xyce_PARALLEL_MPI )
     add_executable( testMPI testBUG967.c )
     target_link_libraries( testMPI XyceLib )
     add_test( testMPI testMPI )
endif()

