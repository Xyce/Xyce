
add_library( xycecinterface EXCLUDE_FROM_ALL )
target_sources( xycecinterface PRIVATE N_CIR_XyceCInterface.C PUBLIC N_CIR_XyceCInterface.h )
target_include_directories( xycecinterface PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> )
target_link_libraries( xycecinterface XyceLib )

if(GTest_FOUND)
     #test executables
     add_executable(XyceCInterfaceUnitTests XyceCInterfaceUnitTests.C)
     target_link_libraries(XyceCInterfaceUnitTests PUBLIC xycecinterface GTest::gtest)

     # Some tests require a data file; put them in the right place
     file(COPY
          ${CMAKE_CURRENT_SOURCE_DIR}/TestNetlist1.cir
          ${CMAKE_CURRENT_SOURCE_DIR}/TestNetlist2.cir
          ${CMAKE_CURRENT_SOURCE_DIR}/TestNetlist3.cir
          DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

     # Wrap the GTest tests with CTest
     gtest_discover_tests( XyceCInterfaceUnitTests TEST_PREFIX XyceCInterfaceUnit: PROPERTIES "LABELS;nightly;LABELS;required:mixedsignal")
endif()

if(BUILD_SHARED_LIBS AND CMAKE_HOST_WIN32 )
  install( TARGETS xycecinterface DESTINATION bin OPTIONAL)
else()
  install( TARGETS xycecinterface DESTINATION lib OPTIONAL)
endif()

configure_file(xyce_interface.py.cmake xyce_interface.py @ONLY )
configure_file(SetupPythonEnv.m.cmake SetupPythonEnv.m @ONLY )
install( FILES N_CIR_XyceCInterface.h DESTINATION include OPTIONAL)
install( FILES ${CMAKE_CURRENT_BINARY_DIR}/xyce_interface.py DESTINATION share OPTIONAL)
install( FILES ${CMAKE_CURRENT_BINARY_DIR}/SetupPythonEnv.m DESTINATION share OPTIONAL)
install( FILES XyceRest.py DESTINATION share OPTIONAL)
install( DIRECTORY PythonExamples DESTINATION share OPTIONAL)
install( DIRECTORY MatlabExamples DESTINATION share OPTIONAL)



