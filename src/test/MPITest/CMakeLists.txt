
# -- build targets -----------------------------------------------------------

# create binary
add_executable( testMPI
                EXCLUDE_FROM_ALL
                ${CMAKE_CURRENT_SOURCE_DIR}/testBUG967.c )

# link against available Xyce library 
if ( Xyce_ENABLE_SHARED )
  target_link_libraries( testMPI lib_xyce_shared )
else ( Xyce_ENABLE_SHARED )
  target_link_libraries( testMPI lib_xyce_static )
endif ( Xyce_ENABLE_SHARED )

