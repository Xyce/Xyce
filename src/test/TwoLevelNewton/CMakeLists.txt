
add_executable(twoLevelNewtonLinearDC EXCLUDE_FROM_ALL twoLevelNewtonLinearDC.C )
target_include_directories(twoLevelNewtonLinearDC PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(twoLevelNewtonLinearDC XyceLib )

add_executable(twoLevelNewtonLinearTran twoLevelNewtonLinearTran.C )
target_include_directories(twoLevelNewtonLinearTran PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(twoLevelNewtonLinearTran XyceLib )


add_executable(TwoLevelTest1 TwoLevelTest1.C )
target_include_directories(TwoLevelTest1 PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(TwoLevelTest1 XyceLib )
if( BUILD_TESTING AND (NOT Xyce_PARALLEL_MPI) )
  add_test( TwoLevelTest1 TwoLevelTest1  )
  set_tests_properties( TwoLevelTest1 PROPERTIES REQUIRED_FILES "compInner.cir" )
  get_target_property(XyceLibDir XyceLib BINARY_DIR )
  set_tests_properties( TwoLevelTest1 PROPERTIES ENVIRONMENT_MODIFICATION "PATH=path_list_prepend:${XyceLibDir}")
  set_tests_properties( TwoLevelTest1 PROPERTIES LABELS "nightly;unittest;twolevel" )
  file( COPY compInner.cir DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
endif()
