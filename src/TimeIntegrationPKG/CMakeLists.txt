
# -- process files for this library ------------------------------------------

# class source list
target_sources(XyceStatic PUBLIC
      N_TIA_DataStore.C
      N_TIA_NoTimeIntegration.C
      N_TIA_Gear12.C
      N_TIA_OneStep.C
      N_TIA_StepErrorControl.C
      N_TIA_TIAParams.C
      N_TIA_TimeIntegrationMethods.C
      N_TIA_RegisterTimeIntegrationMethods.C
      N_TIA_WorkingIntegrationMethod.C )

# add header files to project
target_include_directories( XyceStatic PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# -- build requested libraries -----------------------------------------------
set ( XYCE_PACKAGE_LIBRARIES ${XYCE_PACKAGE_LIBRARIES}
      TimeIntegration
      CACHE INTERNAL "X_L" )

