
# -- process files for this library ------------------------------------------

# class source list
set ( TimeIntegration_SOURCES   
      ${CMAKE_CURRENT_SOURCE_DIR}/N_TIA_DataStore.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_TIA_NoTimeIntegration.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_TIA_Gear12.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_TIA_OneStep.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_TIA_StepErrorControl.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_TIA_TIAParams.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_TIA_TimeIntegrationMethods.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_TIA_RegisterTimeIntegrationMethods.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_TIA_WorkingIntegrationMethod.C 
      CACHE INTERNAL "X_S" )
  

# add header files to project
set ( XYCE_INCLUDE_DIRECTORIES ${XYCE_INCLUDE_DIRECTORIES} 
      ${CMAKE_CURRENT_SOURCE_DIR}
      CACHE INTERNAL "X_I" )




# -- build requested libraries -----------------------------------------------
set ( XYCE_PACKAGE_LIBRARIES ${XYCE_PACKAGE_LIBRARIES} 
      TimeIntegration 
      CACHE INTERNAL "X_L" )

