if(GTest_FOUND)
     #test executables
     add_executable(XyceSimulatorUnitTests XyceSimulatorUnitTests.C)
     target_link_libraries(XyceSimulatorUnitTests PUBLIC XyceLib GTest::gtest)

     # Some tests require a data file; put them in the right place
     file(COPY
          ${CMAKE_CURRENT_SOURCE_DIR}/TestNetlist1.cir
          ${CMAKE_CURRENT_SOURCE_DIR}/TestNetlist2.cir
          ${CMAKE_CURRENT_SOURCE_DIR}/TestNetlist3.cir
          ${CMAKE_CURRENT_SOURCE_DIR}/MultiPort.cir
          ${CMAKE_CURRENT_SOURCE_DIR}/MultiPort2.cir
          DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

     # Wrap the GTest tests with CTest
     # Setting the test labels looks odd but is needed to overcome a bug in gtest-cmake see https://gitlab.kitware.com/cmake/cmake/-/issues/17812
     gtest_discover_tests( XyceSimulatorUnitTests TEST_PREFIX XyceSimulatorUnit: PROPERTIES "LABELS;nightly;LABELS;unittest;LABELS;required:mixedsignal")
endif()

get_target_property(XyceBuildDir Xyce BINARY_DIR)
add_test(NAME RunXyce COMMAND Xyce -capabilities WORKING_DIR ${XyceBuildDir} )
set_tests_properties( RunXyce PROPERTIES LABELS "nightly;unittest" )
