if(TEST_SUITE)
  include(GoogleTest)
  find_package(GTest REQUIRED)
  message(STATUS "Enabling test suite")
  # Parser tests currently cause issue: std::complex<double> to double
  # when this is resolved, uncomment and parser tests should work as intended

  #flex-bison generation
  bison_target(ExpressionParser ExpressionParser.yxx ${CMAKE_CURRENT_BINARY_DIR}/ExpressionParser.cxx)
  flex_target(ExpressionLexer ExpressionLexer.l ${CMAKE_CURRENT_BINARY_DIR}/ExpressionLexer.cxx)

  add_flex_bison_dependency(ExpressionLexer ExpressionParser)

  #test executables
  add_executable(astUnitTest astUnitTest.C)
  target_link_libraries(astUnitTest PUBLIC XyceStatic GTest::GTest)
  add_executable(astCopyConstUnitTest astCopyConstUnitTest.C)
  target_link_libraries(astCopyConstUnitTest PUBLIC XyceStatic GTest::GTest)

  #add_executable(parserUnitTest parserUnitTest.C ExpressionParser ExpressionLexer)
  #target_link_libraries(parserUnitTest PUBLIC ExpressionLexer ExpressionParser XyceStatic GTest::GTest)
  #add_executable(parserCopyUnitTest parserCopyUnitTest.C ExpressionParser ExpressionLexer)
  #target_link_libraries(parserCopyUnitTest PUBLIC ExpressionLexer ExpressionParser XyceStatic GTest::GTest )

  #convert gtest tests to ctest
  enable_testing()
  gtest_discover_tests(astUnitTest TEST_PREFIX astUnit:)
  gtest_discover_tests(astCopyConstUnitTest TEST_PREFIX astCopyConstUnit:)

  #gtest_discover_tests(parserUnitTest TEST_PREFIX parserUnit:)
  #gtest_discover_tests(parserCopyUnitTest TEST_PREFIX parserCopyUnit:)
else()
  message(STATUS "Test suite is not enabled")
endif()
