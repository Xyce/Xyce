
# -- build targets -----------------------------------------------------------


# create binary 1 of 2
add_executable( testHarness1685 
                EXCLUDE_FROM_ALL
                ${CMAKE_CURRENT_SOURCE_DIR}/testHarness1685.C )

# link against available Xyce library 
if ( Xyce_ENABLE_SHARED )
  target_link_libraries( testHarness1685 lib_xyce_shared ${DAKOTA_OBJS} )
else ( Xyce_ENABLE_SHARED )
  target_link_libraries( testHarness1685 lib_xyce_static ${DAKOTA_OBJS} )
endif ( Xyce_ENABLE_SHARED )


