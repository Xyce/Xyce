if(GTest_FOUND)
    #test executables
    add_executable(ParsingHelperTests ParsingHelperTests.C)
    target_link_libraries( ParsingHelperTests PUBLIC XyceLib GTest::gtest)

    #test executables
    add_executable(CmdParseTests CmdParseTests.C)
    target_link_libraries( CmdParseTests PUBLIC XyceLib GTest::gtest)

    gtest_discover_tests(ParsingHelperTests TEST_PREFIX ParsingHelper:)
    gtest_discover_tests(CmdParseTests TEST_PREFIX CmdParse:)
endif()