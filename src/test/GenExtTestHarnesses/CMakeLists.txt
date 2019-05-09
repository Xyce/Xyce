
# -- build targets -----------------------------------------------------------

# create binary
add_executable( testGenCoup
                EXCLUDE_FROM_ALL
                ${CMAKE_CURRENT_SOURCE_DIR}/testGenCoup.C )

# link against available Xyce library 
if ( Xyce_ENABLE_SHARED )
  target_link_libraries( testGenCoup lib_xyce_shared ${DAKOTA_OBJS} )
else ( Xyce_ENABLE_SHARED )
  target_link_libraries( testGenCoup lib_xyce_static ${DAKOTA_OBJS} )
endif ( Xyce_ENABLE_SHARED )

