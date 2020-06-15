
add_executable( bug1574 EXCLUDE_FROM_ALL ${CMAKE_CURRENT_SOURCE_DIR}/bug1574.C )
target_link_libraries( bug1574 XyceStatic )

add_executable( testHarness1685 EXCLUDE_FROM_ALL ${CMAKE_CURRENT_SOURCE_DIR}/testHarness1685.C )
target_link_libraries( testHarness1685 XyceStatic )


#if( Xyce_ENABLED_SHARED )
#  target_link_libraries( bug1574 lib_xyce_shared )
#else( Xyce_ENABLED_SHARED )
#  target_link_libraries( bug1574 lib_xyce_static )
#endif( Xyce_ENABLED_SHARED )

