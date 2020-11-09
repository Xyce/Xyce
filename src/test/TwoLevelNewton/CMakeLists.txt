
add_executable(twoLevelNewtonLinearDC EXCLUDE_FROM_ALL twoLevelNewtonLinearDC.C )
target_include_directories(twoLevelNewtonLinearDC PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(twoLevelNewtonLinearDC XyceStatic )


add_executable(twoLevelNewtonLinearTran EXCLUDE_FROM_ALL twoLevelNewtonLinearTran.C )
target_include_directories(twoLevelNewtonLinearTran PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(twoLevelNewtonLinearTran XyceStatic )


#if( Xyce_ENABLED_SHARED )
#  target_link_librarie(twoLevelNewtonLinearDC lib_xyce_shared )
#else( Xyce_ENABLED_SHARED )
#  target_link_libraries(twoLevelNewtonLinearDC lib_xyce_static )
#endif( Xyce_ENABLED_SHARED )

