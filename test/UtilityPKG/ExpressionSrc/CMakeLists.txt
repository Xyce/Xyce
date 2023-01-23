if(GTest_FOUND)
     set(TARGET_SOURCE_DIR "${CMAKE_SOURCE_DIR}/src/UtilityPKG/ExpressionSrc")

     #test executables
     add_executable(astUnitTest astUnitTest.C)
     target_link_libraries(astUnitTest PUBLIC XyceLib GTest::gtest)

     add_executable(astCopyConstUnitTest astCopyConstUnitTest.C)
     target_link_libraries(astCopyConstUnitTest PUBLIC XyceLib GTest::gtest)

     # The AST tree uses templates, so it needs a type specified in order to use it;
     # but, because of the inflexibility of flex/bison, the type needs to be specified
     # at compile time rather than run time. Currently this is done via the
     # `USE_TYPE_DOUBLE` variable. It is possible that using flex/bison variants will
     # allow templates to be used in a more complete way and then the type can be
     # specified at runtime. At this time, Xyce *requires* the expression
     # library to be in its "complex" form

     # The "complexParserUnitTest" program is written to test the expression
     # library using `std::complex<double>`.
     add_executable(complexParserUnitTest complexParserUnitTest.C)
     target_link_libraries(complexParserUnitTest PUBLIC XyceLib GTest::gtest)

     # The "parserUnitTest" program is written to test the expression library
     # using doubles.  Therefore, the tests will not pass when linking to
     # XyceLib, as is done above.  The following forces a
     # reconstruction/recompile of the parser with USE_TYPE_DOUBLE defined,
     # thus overriding the parser in XyceLib.  I think one should be able to
     # simply recompile the C++ source files created by `XyceExpressionParser"
     # and `XyceExpressionLexer` (as done in the BUILD_ONLY_UNIT_TESTS section,
     # below), but CMake seems to get confused by that.  That could be due to
     # the CMake2-style used in the Bison and flex modules.
     bison_target( ExpressionParserForDouble
          "${TARGET_SOURCE_DIR}/ExpressionParser.yxx"
          "${CMAKE_CURRENT_BINARY_DIR}/ExpressionParser.cxx")
     flex_target( ExpressionLexerForDouble
          "${TARGET_SOURCE_DIR}/ExpressionLexer.l"
          "${CMAKE_CURRENT_BINARY_DIR}/ExpressionLexer.cxx"
          COMPILE_FLAGS "-i -+ --full --prefix=exp")

     add_flex_bison_dependency( ExpressionLexerForDouble ExpressionParserForDouble )

     add_executable(parserUnitTest parserUnitTest.C
          ${BISON_ExpressionParserForDouble_OUTPUTS}
          ${FLEX_ExpressionLexerForDouble_OUTPUTS}
          ${TARGET_SOURCE_DIR}/newExpression.C
          ${TARGET_SOURCE_DIR}/expressionGroup.C)
     target_compile_definitions(parserUnitTest PUBLIC USE_TYPE_DOUBLE)
     target_link_libraries(parserUnitTest PUBLIC XyceLib GTest::gtest)

     # Some tests require a data file; put them in the right place
     file(COPY
          ${CMAKE_CURRENT_SOURCE_DIR}/SubDir1
          ${CMAKE_CURRENT_SOURCE_DIR}/Sub_Dir
          ${CMAKE_CURRENT_SOURCE_DIR}/SubDir1/test1.dat
          ${CMAKE_CURRENT_SOURCE_DIR}/test1-1.dat
          DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

     # Wrap the GTest tests with CTest
     gtest_discover_tests(astUnitTest TEST_PREFIX astUnit:)
     gtest_discover_tests(astCopyConstUnitTest TEST_PREFIX astCopyConstUnit:)
     gtest_discover_tests(complexParserUnitTest TEST_PREFIX complexParserUnit:)
     gtest_discover_tests(parserUnitTest TEST_PREFIX parserUnit:)

endif()

# The following should **NEVER** be enabled with a full Xyce build.
# See the README.md file in this directory.
if(BUILD_ONLY_UNIT_TESTS)

     cmake_minimum_required(VERSION 3.13 FATAL_ERROR)

     # Xyce should use the same compilers as Trilinos. Therefore, Trilinos must be
     # loaded before "project" is called.
     message(STATUS "Looking for required Trilinos package: Teuchos")
     find_package(Trilinos 13.5 REQUIRED  Teuchos)
     message(STATUS "Looking for required Trilinos package: Teuchoss - found")

     # Set the compiler, so "project" finds the right one.
     if( DEFINED CMAKE_CXX_COMPILER )
          message( "Using user defined CXX ${CMAKE_CXX_COMPILER}" )
     else()
          message( "Using Trilinos CXX compiler" )
          set( CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER} )
     endif()

     if( DEFINED CMAKE_C_COMPILER )
          message( "Using user defined CC ${CMAKE_C_COMPILER}" )
     else()
          message( "Using Trilinos C compiler" )
          set( CMAKE_C_COMPILER ${Trilinos_C_COMPILER} )
     endif()

     project(Xyce_Expression_Unit_Tests VERSION 1.0.0 LANGUAGES CXX C)

     set(CMAKE_CXX_STANDARD 11)
     set(CMAKE_CXX_STANDARD_REQUIRED ON)
     set(CMAKE_CXX_EXTENSIONS OFF)

     include(CheckIncludeFileCXX)
     include(GoogleTest)
     find_package(GTest REQUIRED)
     enable_testing()

     # Find flex and Bison
     message(STATUS "Looking for flex and Bison")
     find_package(FLEX REQUIRED)
     find_package(BISON 3.0.4 REQUIRED)

     set(TARGET_SOURCE_DIR "${CMAKE_SOURCE_DIR}/../../../src/UtilityPKG/ExpressionSrc")

     bison_target( XyceExpressionParser
          "${TARGET_SOURCE_DIR}/ExpressionParser.yxx"
          "${CMAKE_CURRENT_BINARY_DIR}/ExpressionParser.cxx")
     flex_target( XyceExpressionLexer
          "${TARGET_SOURCE_DIR}/ExpressionLexer.l"
          "${CMAKE_CURRENT_BINARY_DIR}/ExpressionLexer.cxx"
          COMPILE_FLAGS "-i -+ --full --prefix=exp")

     add_flex_bison_dependency( XyceExpressionLexer XyceExpressionParser  )

     # END Find flex and Bison

     # Several files include Xyce_config.h, though most don't seem to need it.
     # However, N_UTL_Marshal.h does leverage HAVE_SYS_STAT_H, so get that set up.
     check_include_file_cxx( "sys/stat.h" HAVE_SYS_STAT_H )
     configure_file( "${CMAKE_SOURCE_DIR}/Xyce_config.h.cmake"
     "${CMAKE_CURRENT_BINARY_DIR}/Xyce_config.h" )

     # Build a library using the common required source files.
     add_library(parserReqdComplex STATIC
          ${BISON_XyceExpressionParser_OUTPUTS}
          ${FLEX_XyceExpressionLexer_OUTPUTS}
          ${TARGET_SOURCE_DIR}/ast.C
          ${TARGET_SOURCE_DIR}/newExpression.C
          ${TARGET_SOURCE_DIR}/expressionGroup.C
          ${TARGET_SOURCE_DIR}/../N_UTL_CheckIfValidFile.C
          ${TARGET_SOURCE_DIR}/../N_UTL_LogStream.C
          ${TARGET_SOURCE_DIR}/../N_UTL_Marshal.C
          ${TARGET_SOURCE_DIR}/../N_UTL_NetlistLocation.C
          ${TARGET_SOURCE_DIR}/../N_UTL_ReportHandler.C
          ${TARGET_SOURCE_DIR}/../N_UTL_ExtendedString.C
          ${TARGET_SOURCE_DIR}/../N_UTL_NoCase.C
          ${TARGET_SOURCE_DIR}/../../ErrorHandlingPKG/N_ERH_Message.C
          ${TARGET_SOURCE_DIR}/../../ErrorHandlingPKG/N_ERH_Messenger.C
          )

     target_include_directories(parserReqdComplex PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}> 
          $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}> $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/../> 
          $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/../ExpressionSrc/ >
          $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/../../ErrorHandlingPKG/ >
          $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/../../IOInterfacePKG/ >
          $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/../../ParallelDistPKG/ >)

     target_link_libraries(parserReqdComplex PUBLIC Trilinos::all_selected_libs)

     target_link_libraries(parserReqdComplex PUBLIC GTest::gtest)

     # Build the test executables.
     add_executable(astUnitTest astUnitTest.C)
     target_link_libraries(astUnitTest PUBLIC parserReqdComplex)

     add_executable(astCopyConstUnitTest astCopyConstUnitTest.C)
     target_link_libraries(astCopyConstUnitTest PUBLIC parserReqdComplex)

     add_executable(complexParserUnitTest complexParserUnitTest.C)
     target_link_libraries(complexParserUnitTest PUBLIC parserReqdComplex)

     # The parserUnitTest requires USE_TYPE_DOUBLE to be true, so recompile the
     # parser with USE_TYPE_DOUBLE defined.  See the explanation in the
     # BUILD_TESTING section, above.
     add_executable(parserUnitTest parserUnitTest.C
          ${BISON_XyceExpressionParser_OUTPUTS}
          ${FLEX_XyceExpressionLexer_OUTPUTS}
          ${TARGET_SOURCE_DIR}/newExpression.C
          ${TARGET_SOURCE_DIR}/expressionGroup.C
          )
     target_compile_definitions(parserUnitTest PUBLIC USE_TYPE_DOUBLE)
     target_link_libraries(parserUnitTest PUBLIC parserReqdComplex)

     # Some tests require data files; put them in the right place
     file(COPY
          ${CMAKE_SOURCE_DIR}/SubDir1
          ${CMAKE_SOURCE_DIR}/Sub_Dir
          ${CMAKE_SOURCE_DIR}/SubDir1/test1.dat
          ${CMAKE_SOURCE_DIR}/test1-1.dat
          DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

     gtest_discover_tests(astUnitTest TEST_PREFIX astUnit:)
     gtest_discover_tests(astCopyConstUnitTest TEST_PREFIX astCopyConstUnit:)
     gtest_discover_tests(complexParserUnitTest TEST_PREFIX complexParserUnit:)
     gtest_discover_tests(parserUnitTest TEST_PREFIX parserUnit:)

endif()
