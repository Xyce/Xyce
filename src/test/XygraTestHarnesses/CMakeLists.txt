
# -- build targets -----------------------------------------------------------

# create binary
add_executable( XygraSource
                EXCLUDE_FROM_ALL
                ${CMAKE_CURRENT_SOURCE_DIR}/XygraSource.C )

# link against available Xyce library 
if ( Xyce_ENABLE_SHARED )
  target_link_libraries( XygraSource lib_xyce_shared ${DAKOTA_OBJS} )
else ( Xyce_ENABLE_SHARED )
  target_link_libraries( XygraSource lib_xyce_static ${DAKOTA_OBJS} )
endif ( Xyce_ENABLE_SHARED )

