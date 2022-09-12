
if( BUILD_TESTING )
     add_executable( testFFTInterface testFFTInterface.C )
     target_link_libraries( testFFTInterface XyceLib )
     add_test( testFFTInterface testFFTInterface FFTInterface.cir )
     set_tests_properties( testFFTInterface PROPERTIES REQUIRED_FILES "FFTInterface.cir" )
     file( COPY FFTInterface.cir DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
endif()

