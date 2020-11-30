if(TEST_SUITE)
     include(GoogleTest)
     find_package(GTest REQUIRED)
     message(STATUS "Enabling test suite")

     #test executables
     add_executable(astUnitTest astUnitTest.C)
     target_link_libraries(astUnitTest PUBLIC XyceStatic GTest::GTest)

     add_executable(astCopyConstUnitTest astCopyConstUnitTest.C)
     target_link_libraries(astCopyConstUnitTest PUBLIC XyceStatic GTest::GTest)

     add_executable(complexParserUnitTest complexParserUnitTest.C)
     target_link_libraries(complexParserUnitTest PUBLIC XyceStatic GTest::GTest)
     # One complexParserUnitTest test requires a data file; put it in the right place
     file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/SubDir1/test1.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

     # The AST tree uses templates, so it needs a type specified in order to use it;
     # but, because of the inflexibility of flex/bison, the type needs to be specified
     # at compile time rather than run time. Currently this is done via the
     # `USE_TYPE_DOUBLE` variable. It is possible that using flex/bison variants will
     # allow templates to be used in a more complete way and then the type can be
     # specified at runtime. For the time being, though, Xyce requires the expression
     # library to be in its "complex" form, so "parserUnitTest" cannot currently be
     # run with the full Xyce build.
     #
     # The other parser unit test program, "complexParserUnitTest" (above) is
     # written to test out the expression library with `std::complex<double>`.
     #
     #    add_executable(parserUnitTest parserUnitTest.C)
     #    target_link_libraries(parserUnitTest PUBLIC XyceStatic GTest::GTest)
     #    # Some parserUnitTest tests require data files; put them in the right place
     #    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Sub_Dir ${CMAKE_CURRENT_BINARY_DIR}/SubDir1)
     #    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/SubDir1/test1.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/SubDir1)
     #    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/Sub_Dir/1test_5.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/Sub_Dir)

     # Wrap the GTest tests with CTest
     enable_testing()
     gtest_discover_tests(astUnitTest TEST_PREFIX astUnit:)
     gtest_discover_tests(astCopyConstUnitTest TEST_PREFIX astCopyConstUnit:)
     gtest_discover_tests(complexParserUnitTest TEST_PREFIX complexParserUnit:)
     #    gtest_discover_tests(parserUnitTest TEST_PREFIX parserUnit:)

else()
     message(STATUS "Test suite is not enabled")
endif()

if(BUILD_ONLY_UNIT_TESTS)
     # The following is ONLY to be run as an independent project.

     cmake_minimum_required(VERSION 3.13 FATAL_ERROR)

     # Xyce should use the same compilers as Trilinos. Therefore, Trilinos must be
     # loaded before "project" is called.
     message(STATUS "Looking for required Trilinos package: Teuchos")
     find_package(Trilinos 12.12 REQUIRED  Teuchos)
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

     # Find flex and Bison
     message(STATUS "Looking for flex and Bison")
     find_package(FLEX REQUIRED)
     find_package(BISON 3.0.4 REQUIRED)

     bison_target( XyceExpressionParser
          "${CMAKE_SOURCE_DIR}/../ExpressionSrc/ExpressionParser.yxx"
          "${CMAKE_CURRENT_BINARY_DIR}/ExpressionParser.cxx")
     flex_target( XyceExpressionLexer
          "${CMAKE_SOURCE_DIR}/../ExpressionSrc/ExpressionLexer.l"
          "${CMAKE_CURRENT_BINARY_DIR}/ExpressionLexer.cxx"
          COMPILE_FLAGS "-i -+ --full --prefix=exp")

     add_flex_bison_dependency( XyceExpressionLexer XyceExpressionParser  )

     # END Find flex and Bison

     # Set up a single ifdef
     check_include_file_cxx( "sys/stat.h" HAVE_SYS_STAT_H )
     configure_file( "${CMAKE_SOURCE_DIR}/Xyce_config.h.cmake"
          "${CMAKE_CURRENT_BINARY_DIR}/Xyce_config.h" )

     add_library(parserReqdComplex STATIC
          ${BISON_XyceExpressionParser_OUTPUTS}
          ${FLEX_XyceExpressionLexer_OUTPUTS}
          ${CMAKE_SOURCE_DIR}/../ExpressionSrc/ast.C
          ${CMAKE_SOURCE_DIR}/../ExpressionSrc/newExpression.C
          ${CMAKE_SOURCE_DIR}/../N_UTL_CheckIfValidFile.C
          ${CMAKE_SOURCE_DIR}/../N_UTL_LogStream.C
          ${CMAKE_SOURCE_DIR}/../N_UTL_Marshal.C
          ${CMAKE_SOURCE_DIR}/../N_UTL_NetlistLocation.C
          ${CMAKE_SOURCE_DIR}/../N_UTL_ReportHandler.C
          ${CMAKE_SOURCE_DIR}/../N_UTL_ExtendedString.C
          ${CMAKE_SOURCE_DIR}/../N_UTL_NoCase.C
          ${CMAKE_SOURCE_DIR}/../../ErrorHandlingPKG/N_ERH_Message.C
          ${CMAKE_SOURCE_DIR}/../../ErrorHandlingPKG/N_ERH_Messenger.C
          )

     target_include_directories(parserReqdComplex PUBLIC ${CMAKE_CURRENT_BINARY_DIR}
          ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/../
          ${CMAKE_SOURCE_DIR}/../ExpressionSrc/
          ${CMAKE_SOURCE_DIR}/../../ErrorHandlingPKG/
          ${CMAKE_SOURCE_DIR}/../../ParallelDistPKG/
          ${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS})

     target_link_libraries(parserReqdComplex PUBLIC
          ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES} GTest::GTest)

     add_executable(astUnitTest astUnitTest.C)
     target_include_directories(astUnitTest PUBLIC ${CMAKE_CURRENT_BINARY_DIR}
          ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/../
          ${CMAKE_SOURCE_DIR}/../ExpressionSrc/
          ${CMAKE_SOURCE_DIR}/../../ErrorHandlingPKG/
          ${CMAKE_SOURCE_DIR}/../../ParallelDistPKG/
          ${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS})
     target_link_libraries(astUnitTest PUBLIC parserReqdComplex
          ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES} GTest::GTest)

     add_executable(astCopyConstUnitTest astCopyConstUnitTest.C)
     target_include_directories(astCopyConstUnitTest PUBLIC ${CMAKE_CURRENT_BINARY_DIR}
          ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/../
          ${CMAKE_SOURCE_DIR}/../ExpressionSrc/
          ${CMAKE_SOURCE_DIR}/../../ErrorHandlingPKG/
          ${CMAKE_SOURCE_DIR}/../../ParallelDistPKG/
          ${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS})
     target_link_libraries(astCopyConstUnitTest PUBLIC parserReqdComplex
          ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES} GTest::GTest)

     add_executable(complexParserUnitTest complexParserUnitTest.C)
     target_include_directories(complexParserUnitTest PUBLIC ${CMAKE_CURRENT_BINARY_DIR}
          ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/../
          ${CMAKE_SOURCE_DIR}/../ExpressionSrc/
          ${CMAKE_SOURCE_DIR}/../../ErrorHandlingPKG/
          ${CMAKE_SOURCE_DIR}/../../ParallelDistPKG/
          ${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS})
     target_link_libraries(complexParserUnitTest PUBLIC parserReqdComplex
          ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES} GTest::GTest)

     # The parserUnitTest requires USE_TYPE_DOUBLE to be true, so the build is
     # incompatible with the above library. Therefore, all the sources need to
     # be rebuilt.
     add_executable(parserUnitTest parserUnitTest.C
          ${BISON_XyceExpressionParser_OUTPUTS}
          ${FLEX_XyceExpressionLexer_OUTPUTS}
          ${CMAKE_SOURCE_DIR}/../ExpressionSrc/newExpression.C
          )
     target_compile_options(parserUnitTest PUBLIC -DUSE_TYPE_DOUBLE)
     target_include_directories(parserUnitTest PUBLIC ${CMAKE_CURRENT_BINARY_DIR}
          ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/../
          ${CMAKE_SOURCE_DIR}/../ExpressionSrc/
          ${CMAKE_SOURCE_DIR}/../../ErrorHandlingPKG/
          ${CMAKE_SOURCE_DIR}/../../ParallelDistPKG/
          ${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS})
     target_link_libraries(parserUnitTest PUBLIC parserReqdComplex
          ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES} GTest::GTest)

     # Some tests require data files; put them in the right place
     file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Sub_Dir ${CMAKE_CURRENT_BINARY_DIR}/SubDir1)
     file(COPY ${CMAKE_SOURCE_DIR}/SubDir1/test1.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/SubDir1)
     file(COPY ${CMAKE_SOURCE_DIR}/Sub_Dir/1test_5.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/Sub_Dir)
     file(COPY ${CMAKE_SOURCE_DIR}/SubDir1/test1.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

     enable_testing()
     gtest_discover_tests(astUnitTest TEST_PREFIX astUnit:)
     gtest_discover_tests(astCopyConstUnitTest TEST_PREFIX astCopyConstUnit:)
     gtest_discover_tests(complexParserUnitTest TEST_PREFIX complexParserUnit:)
     gtest_discover_tests(parserUnitTest TEST_PREFIX parserUnit:)

endif()
