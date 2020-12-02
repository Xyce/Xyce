
add_executable( testFFTInterface EXCLUDE_FROM_ALL testFFTInterface.C )
target_link_libraries( testFFTInterface XyceLib )

#if( Xyce_ENABLED_SHARED )
#  target_link_libraries( testFFTInterface lib_xyce_shared )
#else( Xyce_ENABLED_SHARED )
#  target_link_libraries( testFFTInterface lib_xyce_static )
#endif( Xyce_ENABLED_SHARED )

