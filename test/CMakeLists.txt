if(NOT Xyce_PARALLEL_MPI)
  add_subdirectory(CircuitPKG)
endif()

if (Xyce_REGRESSION AND DEFINED Xyce_REGRESSION_DIR AND EXISTS "${Xyce_REGRESSION_DIR}/CMakeLists.txt")
  add_subdirectory (${Xyce_REGRESSION_DIR} ${CMAKE_CURRENT_BINARY_DIR}/Xyce_Regression)
endif ()
add_subdirectory(IOInterfacePKG)
add_subdirectory(UtilityPKG)
