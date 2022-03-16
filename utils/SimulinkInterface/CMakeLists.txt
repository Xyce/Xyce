

if( Matlab_FOUND )
  matlab_add_mex( NAME xyce_sfunction MODULE  SRC xyce_sfunction.cpp LINK_TO XyceLib )
  target_include_directories(xyce_sfunction PUBLIC
      ${CMAKE_CURRENT_SOURCE_DIR}
      ${Matlab_INCLUDE_DIRS} "${Matlab_ROOT_DIR}/simulink/include")
  install(TARGETS xyce_sfunction DESTINATION lib)
  install(FILES XyceBlocks.slx DESTINATION share/simulink)
  
  matlab_add_mex( NAME xyce2nLv_sfunction MODULE  SRC xyce2nLv_sfunction.cpp LINK_TO XyceLib )
  target_include_directories(xyce2nLv_sfunction PUBLIC
      ${CMAKE_CURRENT_SOURCE_DIR}
      ${Matlab_INCLUDE_DIRS} "${Matlab_ROOT_DIR}/simulink/include")
  
  #  EXCLUDE_FROM_ALL
  #target_link_libraries( xyce_sfunction XyceStatic )

endif()

