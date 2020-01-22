
add_executable( testFFTInterface EXCLUDE_FROM_ALL ${CMAKE_CURRENT_SOURCE_DIR}/testFFTInterface.C )
target_link_libraries( testFFTInterface XyceStatic )

#if( Xyce_ENABLED_SHARED )
#  target_link_libraries( testFFTInterface lib_xyce_shared )
#else( Xyce_ENABLED_SHARED )
#  target_link_libraries( testFFTInterface lib_xyce_static )
#endif( Xyce_ENABLED_SHARED )

