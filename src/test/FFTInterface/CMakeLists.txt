


if(BUILD_TESTING AND GTest_FOUND)
     #test executables
     add_executable(FFTUnitTests FFTUnitTests.C)
     target_link_libraries(FFTUnitTests PUBLIC XyceLib GTest::gtest)

     # copy input data files over so this test can be run from the build directory too
     file( COPY 
       ${CMAKE_CURRENT_SOURCE_DIR}/Sin1fEven.txt 
       ${CMAKE_CURRENT_SOURCE_DIR}/Sin1fOdd.txt 
       ${CMAKE_CURRENT_SOURCE_DIR}/Sin2fEven.txt 
       ${CMAKE_CURRENT_SOURCE_DIR}/Sin2fOdd.txt 
       ${CMAKE_CURRENT_SOURCE_DIR}/Sin3fEven.txt 
       ${CMAKE_CURRENT_SOURCE_DIR}/Sin3fOdd.txt 
       ${CMAKE_CURRENT_SOURCE_DIR}/Gauss1Even.txt 
       ${CMAKE_CURRENT_SOURCE_DIR}/Gauss1Odd.txt 
       ${CMAKE_CURRENT_SOURCE_DIR}/Gauss2Even.txt 
       ${CMAKE_CURRENT_SOURCE_DIR}/Gauss2Odd.txt 
       ${CMAKE_CURRENT_SOURCE_DIR}/Step1Even.txt 
       ${CMAKE_CURRENT_SOURCE_DIR}/Step1Odd.txt 
       DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
       
     # Wrap the GTest tests with CTest
     # Setting the test labels looks odd but is needed to overcome a bug in gtest-cmake see https://gitlab.kitware.com/cmake/cmake/-/issues/17812
     gtest_discover_tests( FFTUnitTests TEST_PREFIX FFTUnitTests: PROPERTIES "LABELS;nightly;LABELS;unittest;LABELS;fft;WORKING_DIRECTORY;${CMAKE_CURRENT_SOURCE_DIR}")
endif()

if( BUILD_TESTING )
     add_executable( testFFTInterface testFFTInterface.C )
     target_link_libraries( testFFTInterface XyceLib )
     add_test( testFFTInterface testFFTInterface FFTInterface.cir )
     set_tests_properties( testFFTInterface PROPERTIES REQUIRED_FILES "FFTInterface.cir" )
     get_target_property(XyceLibDir XyceLib BINARY_DIR )
     set_tests_properties( testFFTInterface PROPERTIES ENVIRONMENT_MODIFICATION "PATH=path_list_prepend:${XyceLibDir}")
     set_tests_properties( testFFTInterface PROPERTIES LABELS "nightly;unittest;fft" )
     file( COPY FFTInterface.cir DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
endif()

