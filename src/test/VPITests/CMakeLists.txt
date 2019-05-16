
# -- build targets -----------------------------------------------------------

# create binary
add_executable( runXyce 
                EXCLUDE_FROM_ALL
                ${CMAKE_CURRENT_SOURCE_DIR}/runXyce.c )
add_executable( runXyceInSteps 
                EXCLUDE_FROM_ALL
                ${CMAKE_CURRENT_SOURCE_DIR}/runXyceInSteps.c )
add_executable( runXyceWithDAC 
                EXCLUDE_FROM_ALL
                ${CMAKE_CURRENT_SOURCE_DIR}/runXyceWithDAC.c )

# link against available Xyce library 
if ( Xyce_ENABLE_SHARED )
  target_link_libraries( runXyce lib_xyce_shared ${DAKOTA_OBJS} )
  target_link_libraries( runXyceInSteps lib_xyce_shared ${DAKOTA_OBJS} )
  target_link_libraries( runXyceWithDAC lib_xyce_shared ${DAKOTA_OBJS} )
else ( Xyce_ENABLE_SHARED )
  target_link_libraries( runXyce lib_xyce_static ${DAKOTA_OBJS} )
  target_link_libraries( runXyceInSteps lib_xyce_static ${DAKOTA_OBJS} )
  target_link_libraries( runXyceWithDAC lib_xyce_static ${DAKOTA_OBJS} )
endif ( Xyce_ENABLE_SHARED )

