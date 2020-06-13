
add_executable( Xygra EXCLUDE_FROM_ALL ${CMAKE_CURRENT_SOURCE_DIR}/Xygra.C )
target_link_libraries( Xygra XyceStatic )

add_executable( XygraFile EXCLUDE_FROM_ALL ${CMAKE_CURRENT_SOURCE_DIR}/XygraFile.C )
target_link_libraries( XygraFile XyceStatic )

add_executable( XygraSource EXCLUDE_FROM_ALL ${CMAKE_CURRENT_SOURCE_DIR}/XygraSource.C )
target_link_libraries( XygraSource XyceStatic )



#if( Xyce_ENABLED_SHARED )
#  target_link_librarie( Xygra lib_xyce_shared )
#else( Xyce_ENABLED_SHARED )
#  target_link_libraries( Xygra lib_xyce_static )
#endif( Xyce_ENABLED_SHARED )

