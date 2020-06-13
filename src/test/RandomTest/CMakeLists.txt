
add_executable( testRandom EXCLUDE_FROM_ALL ${CMAKE_CURRENT_SOURCE_DIR}/testRandom.C )
target_link_libraries( testRandom XyceStatic )

#if( Xyce_ENABLED_SHARED )
#  target_link_librarie( testRandom lib_xyce_shared )
#else( Xyce_ENABLED_SHARED )
#  target_link_libraries( testRandom lib_xyce_static )
#endif( Xyce_ENABLED_SHARED )

