
# -- build targets -----------------------------------------------------------

# create binary
add_executable( testRandom 
                EXCLUDE_FROM_ALL
                ${CMAKE_CURRENT_SOURCE_DIR}/testRandom.C )

# link against available Xyce library 
if ( Xyce_ENABLE_SHARED )
  target_link_libraries( testRandom lib_xyce_shared ${DAKOTA_OBJS} )
else ( Xyce_ENABLE_SHARED )
  target_link_libraries( testRandom lib_xyce_static ${DAKOTA_OBJS} )
endif ( Xyce_ENABLE_SHARED )

