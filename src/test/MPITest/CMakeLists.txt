
if( Xyce_PARALLEL_MPI )
     add_executable( testMPI testBUG967.c )
     target_link_libraries( testMPI XyceLib )
#     if( BUILD_TESTING )
#          add_test( testMPI testMPI )
#          set_tests_properties( testMPI PROPERTIES LABELS "parallel;unittest" )
#     endif()
endif()

