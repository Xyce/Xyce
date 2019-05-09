
# -- build targets -----------------------------------------------------------


# create binary
add_executable( DeviceInterfaceTest
                EXCLUDE_FROM_ALL
                ${CMAKE_CURRENT_SOURCE_DIR}/DeviceInterfaceTest.C )

# link against available Xyce library 
if ( Xyce_ENABLE_SHARED )
  target_link_libraries( DeviceInterfaceTest lib_xyce_shared ${DAKOTA_OBJS} )
else ( Xyce_ENABLE_SHARED )
  target_link_libraries( DeviceInterfaceTest lib_xyce_static ${DAKOTA_OBJS} )
endif ( Xyce_ENABLE_SHARED )
