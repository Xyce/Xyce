if(Xyce_TEST_SUITE)
     include(GoogleTest)
     find_package(GTest REQUIRED)
     message(STATUS "Enabling test suite")

     #test executables
     add_executable(XyceSimulatorUnitTests XyceSimulatorUnitTests.C)
     target_link_libraries(XyceSimulatorUnitTests PUBLIC XyceLib GTest::GTest)

     # Some tests require a data file; put them in the right place
     file(COPY
          ${CMAKE_CURRENT_SOURCE_DIR}/TestNetlist1.cir
          ${CMAKE_CURRENT_SOURCE_DIR}/TestNetlist2.cir
          ${CMAKE_CURRENT_SOURCE_DIR}/TestNetlist3.cir
          DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

     # Wrap the GTest tests with CTest
     enable_testing()
     gtest_discover_tests( XyceSimulatorUnitTests TEST_PREFIX XyceSimulatorUnit:)

else()
     message(STATUS "Test suite is not enabled")
endif()


