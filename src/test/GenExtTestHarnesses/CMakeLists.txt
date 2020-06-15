
add_executable( testGenCoup EXCLUDE_FROM_ALL ${CMAKE_CURRENT_SOURCE_DIR}/testGenCoup.C )
target_link_libraries( testGenCoup XyceStatic )

#if( Xyce_ENABLED_SHARED )
#  target_link_libraries( testGenCoup lib_xyce_shared )
#else( Xyce_ENABLED_SHARED )
#  target_link_libraries( testGenCoup lib_xyce_static )
#endif( Xyce_ENABLED_SHARED )

