if (Xyce_REGRESSION AND DEFINED Xyce_REGRESSION_DIR AND EXISTS "${Xyce_REGRESSION_DIR}/CMakeLists.txt")
  add_subdirectory (${Xyce_REGRESSION_DIR} ${CMAKE_CURRENT_BINARY_DIR}/Xyce_Regression)
endif ()

# unit tests
if(Xyce_BUILD_UNIT_TESTS)
  if(NOT Xyce_PARALLEL_MPI)
    add_subdirectory(CircuitPKG)
  endif()

  add_subdirectory(IOInterfacePKG)
  add_subdirectory(UtilityPKG)
endif()
