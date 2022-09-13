
add_executable( testBlockLinearSystems EXCLUDE_FROM_ALL testBlockLinearSystems.C )
target_link_libraries( testBlockLinearSystems XyceLib )
#if( BUILD_TESTING )
#     add_test( testBlockLinearSystems testBlockLinearSystems )
#endif()

