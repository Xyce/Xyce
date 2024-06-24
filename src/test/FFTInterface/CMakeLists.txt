
add_executable( testFFTInterface testFFTInterface.C )
target_link_libraries( testFFTInterface XyceLib )

if( BUILD_TESTING )
     add_test( testFFTInterface testFFTInterface FFTInterface.cir )
     set_tests_properties( testFFTInterface PROPERTIES REQUIRED_FILES "FFTInterface.cir" )
     get_target_property(XyceLibDir XyceLib BINARY_DIR )
     set_tests_properties( testFFTInterface PROPERTIES ENVIRONMENT_MODIFICATION "PATH=path_list_prepend:${XyceLibDir}")
     set_tests_properties( testFFTInterface PROPERTIES LABELS "nightly;unittest;fft" )
     file( COPY FFTInterface.cir DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
endif()

