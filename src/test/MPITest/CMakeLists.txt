
add_executable( testMPI EXCLUDE_FROM_ALL ${CMAKE_CURRENT_SOURCE_DIR}/testBUG967.c )
target_link_libraries( testMPI XyceStatic )

#if( Xyce_ENABLED_SHARED )
#  target_link_libraries( testMPI lib_xyce_shared )
#else( Xyce_ENABLED_SHARED )
#  target_link_libraries( testMPI lib_xyce_static )
#endif( Xyce_ENABLED_SHARED )

