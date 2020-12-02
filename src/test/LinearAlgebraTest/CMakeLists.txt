
add_executable( testBlockLinearSystems EXCLUDE_FROM_ALL testBlockLinearSystems.C )
target_link_libraries( testBlockLinearSystems XyceLib )

#if( Xyce_ENABLED_SHARED )
#  target_link_libraries( testBlockLinearSystems lib_xyce_shared )
#else( Xyce_ENABLED_SHARED )
#  target_link_libraries( testBlockLinearSystems lib_xyce_static )
#endif( Xyce_ENABLED_SHARED )

