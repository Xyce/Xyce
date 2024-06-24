add_subdirectory(ExpressionSrc)

if(GTest_FOUND)
    # Get system information to compare with with test functions
    cmake_host_system_information(RESULT systemInfo QUERY HOSTNAME FQDN OS_NAME OS_VERSION OS_PLATFORM PROCESSOR_NAME PROCESSOR_DESCRIPTION)
    list(GET systemInfo 0 VAR)
    add_compile_definitions(UTL_Platform_Tests PUBLIC EXPECTED_HOSTNAME="${VAR}")

    # Test Executables with XyceLib
    add_executable(UTL_Platform_Tests N_UTL_Platform_UnitTests.C)
    target_link_libraries(UTL_Platform_Tests PUBLIC XyceLib GTest::gtest)

    # Test Executables with mocks
    add_executable(UTL_Platform_System_Error_Tests N_UTL_Platform_System_Errors_UnitTests.C N_UTL_Platform_System_Error_Mocks.C ${CMAKE_SOURCE_DIR}/src/UtilityPKG/N_UTL_Platform.C)
    target_include_directories(UTL_Platform_System_Error_Tests PRIVATE ${CMAKE_SOURCE_DIR}/src/UtilityPKG)
    target_include_directories(UTL_Platform_System_Error_Tests PUBLIC $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>/src)
    target_link_libraries(UTL_Platform_System_Error_Tests PUBLIC GTest::gtest)

    add_executable(UTL_Platform_Config_Error_Tests ${CMAKE_SOURCE_DIR}/src/UtilityPKG/N_UTL_Platform.C N_UTL_Platform_Config_Errors_UnitTests.C)
    target_compile_definitions(UTL_Platform_Config_Error_Tests PUBLIC Xyce_CONFIG_H)
    target_include_directories(UTL_Platform_Config_Error_Tests PUBLIC $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>/src)
    target_include_directories(UTL_Platform_Config_Error_Tests PRIVATE ${CMAKE_SOURCE_DIR}/src/UtilityPKG)
    target_link_libraries(UTL_Platform_Config_Error_Tests PUBLIC GTest::gtest)
    
    # Wrap the GTest tests with CTest
    # Setting the test labels looks odd but is needed to overcome a bug in gtest-cmake see https://gitlab.kitware.com/cmake/cmake/-/issues/17812
    gtest_discover_tests(UTL_Platform_Tests TEST_PREFIX UTL_Platform: PROPERTIES "LABELS;nightly;LABELS;unittest;LABELS;platform")
    gtest_discover_tests(UTL_Platform_System_Error_Tests TEST_PREFIX UTL_Platform_Errors: PROPERTIES "LABELS;nightly;LABELS;unittest;LABELS;platform")
    gtest_discover_tests(UTL_Platform_Config_Error_Tests TEST_PREFIX UTL_Platform_Errors: PROPERTIES "LABELS;nightly;LABELS;unittest;LABELS;platform")
endif()