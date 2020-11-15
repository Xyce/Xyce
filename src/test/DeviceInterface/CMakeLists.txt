
add_executable( DeviceInterfaceTest EXCLUDE_FROM_ALL DeviceInterfaceTest.C )
target_link_libraries( DeviceInterfaceTest XyceStatic )

#if( Xyce_ENABLED_SHARED )
#  target_link_librarie( DeviceInterfaceTest lib_xyce_shared )
#else( Xyce_ENABLED_SHARED )
#  target_link_libraries( DeviceInterfaceTest lib_xyce_static )
#endif( Xyce_ENABLED_SHARED )

