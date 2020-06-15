
add_executable( DeviceInterfaceTest EXCLUDE_FROM_ALL ${CMAKE_CURRENT_SOURCE_DIR}/DeviceInterfaceTest.C )
target_link_libraries( DeviceInterfaceTest XyceStatic )

#if( Xyce_ENABLED_SHARED )
#  target_link_librarie( DeviceInterfaceTest lib_xyce_shared )
#else( Xyce_ENABLED_SHARED )
#  target_link_libraries( DeviceInterfaceTest lib_xyce_static )
#endif( Xyce_ENABLED_SHARED )

