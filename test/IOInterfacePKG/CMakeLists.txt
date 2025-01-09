if(GTest_FOUND)
    #test executables
    add_executable(ParsingHelperTests ParsingHelperTests.C)
    target_link_libraries( ParsingHelperTests PUBLIC XyceLib GTest::gtest)

    #test executables
    add_executable(CmdParseTests CmdParseTests.C)
    target_link_libraries( CmdParseTests PUBLIC XyceLib GTest::gtest)

    # Wrap the GTest tests with CTest
    # Setting the test labels looks odd but is needed to overcome a bug in gtest-cmake see https://gitlab.kitware.com/cmake/cmake/-/issues/17812
    gtest_discover_tests(ParsingHelperTests TEST_PREFIX ParsingHelper: PROPERTIES "LABELS;nightly;LABELS;unittest;LABELS;parsing")
    gtest_discover_tests(CmdParseTests TEST_PREFIX CmdParse: PROPERTIES "LABELS;nightly;LABELS;unittest;LABELS;parsing")
    
endif()