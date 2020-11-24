if(TEST_SUITE)
     include(GoogleTest)
     find_package(GTest REQUIRED)
     message(STATUS "Enabling test suite")

     #test executables
     add_executable(astUnitTest astUnitTest.C)
     target_link_libraries(astUnitTest PUBLIC XyceStatic GTest::GTest)
     add_executable(astCopyConstUnitTest astCopyConstUnitTest.C)
     target_link_libraries(astCopyConstUnitTest PUBLIC XyceStatic GTest::GTest)

     # USE_TYPE_DOUBLE affects how the expression library is compiled. It
     # determines whether the type associated with the AST tree is double or
     # std::complex<double>. Right now it has to be defined at compile time.
     # Maybe, eventually, it can be specified at run time.
     #
     # Additionally, a version of the expression library based around double is
     # not currently compatible with the rest of Xyce. So, if USE_TYPE_DOUBLE
     # is defined, the parserUnitTest tests might run, but the rest of Xyce
     # won't work. This if/else is a placeholder to take that information into
     # account, but it will likely never be used in the current form. If it is
     # used, then USE_TYPE_DOUBLE will have to be added in a set() call to
     # cmake/feature_modes.cmake and in a line of src/Xyce_config.h.cmake.
     if (USE_TYPE_DOUBLE)
          #The parserUnitTest tests require double.
          add_executable(parserUnitTest parserUnitTest.C)
          target_link_libraries(parserUnitTest PUBLIC XyceStatic GTest::GTest)
     else()
          # The complexParserUnitTest tests require std::complex<double>.
          add_executable(complexParserUnitTest complexParserUnitTest.C)
          target_link_libraries(complexParserUnitTest PUBLIC XyceStatic GTest::GTest)
     endif()

     # The following is not particularly useful any more, so commenting out for now
     #add_executable(parserCopyUnitTest parserCopyUnitTest.C)
     #target_link_libraries(parserCopyUnitTest PUBLIC XyceStatic GTest::GTest )

     #convert gtest tests to ctest
     enable_testing()
     gtest_discover_tests(astUnitTest TEST_PREFIX astUnit:)
     gtest_discover_tests(astCopyConstUnitTest TEST_PREFIX astCopyConstUnit:)

     if (USE_TYPE_DOUBLE)
          gtest_discover_tests(parserUnitTest TEST_PREFIX parserUnit:)
     else()
          gtest_discover_tests(complexParserUnitTest TEST_PREFIX complexParserUnit:)
     endif()

     # The following is not particularly useful any more, so commenting out for now
     #gtest_discover_tests(parserCopyUnitTest TEST_PREFIX cparserCopyUnit:)
else()
     message(STATUS "Test suite is not enabled")
endif()
