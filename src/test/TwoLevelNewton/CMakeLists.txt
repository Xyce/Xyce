
if( BUILD_TESTING )

     add_executable(twoLevelNewtonLinearDC EXCLUDE_FROM_ALL twoLevelNewtonLinearDC.C )
     target_include_directories(twoLevelNewtonLinearDC PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
     target_link_libraries(twoLevelNewtonLinearDC XyceLib )

     add_executable(twoLevelNewtonLinearTran twoLevelNewtonLinearTran.C )
     target_include_directories(twoLevelNewtonLinearTran PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
     target_link_libraries(twoLevelNewtonLinearTran XyceLib )

endif()

