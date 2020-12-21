
add_executable( testGenCoup EXCLUDE_FROM_ALL testGenCoup.C )
target_link_libraries( testGenCoup XyceLib )

#if( Xyce_ENABLED_SHARED )
#  target_link_libraries( testGenCoup lib_xyce_shared )
#else( Xyce_ENABLED_SHARED )
#  target_link_libraries( testGenCoup lib_xyce_static )
#endif( Xyce_ENABLED_SHARED )

