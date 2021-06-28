if(Xyce_TEST_SUITE)
     include(GoogleTest)
     find_package(GTest REQUIRED)
     message(STATUS "Enabling test suite")

     #test executables
     add_executable(ParsingHelperTests ParsingHelperTests.C)
     target_link_libraries( ParsingHelperTests PUBLIC XyceLib GTest::GTest)

     #test executables
     add_executable(CmdParseTests CmdParseTests.C)
     target_link_libraries( CmdParseTests PUBLIC XyceLib GTest::GTest)
endif()
