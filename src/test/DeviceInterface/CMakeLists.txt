
add_executable( DeviceInterfaceTest DeviceInterfaceTest.C )
target_link_libraries( DeviceInterfaceTest XyceLib )

add_executable( ADC_DACRunTest ADC_DACRunTest.C )
target_link_libraries( ADC_DACRunTest XyceLib )

if( BUILD_TESTING )
     add_test( DeviceInterfaceTest DeviceInterfaceTest ADC_DACtest.cir )
     set_tests_properties( DeviceInterfaceTest PROPERTIES REQUIRED_FILES "ADC_DACtest.cir" )
     file( COPY ADC_DACtest.cir DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

     add_test( ADC_DACRunTest ADC_DACRunTest ADC_DACRunTest.cir )
     set_tests_properties( ADC_DACRunTest PROPERTIES REQUIRED_FILES "ADC_DACRunTest.cir")
     file( COPY ADC_DACRunTest.cir DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
endif()

