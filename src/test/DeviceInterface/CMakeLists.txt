
add_executable( DeviceInterfaceTest EXCLUDE_FROM_ALL DeviceInterfaceTest.C )
target_link_libraries( DeviceInterfaceTest XyceLib )

add_executable( ADC_DACRunTest EXCLUDE_FROM_ALL ADC_DACRunTest.C )
target_link_libraries( ADC_DACRunTest XyceLib )


#if( Xyce_ENABLED_SHARED )
#  target_link_librarie( DeviceInterfaceTest lib_xyce_shared )
#else( Xyce_ENABLED_SHARED )
#  target_link_libraries( DeviceInterfaceTest lib_xyce_static )
#endif( Xyce_ENABLED_SHARED )

