#
# CMakeLists.txt for TimeIntegrationPKG
#

set ( TimeIntegration_SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/N_TIA_DataStore.C
	${CMAKE_CURRENT_SOURCE_DIR}/N_TIA_NoTimeIntegration.C
	${CMAKE_CURRENT_SOURCE_DIR}/N_TIA_Gear12.C
	${CMAKE_CURRENT_SOURCE_DIR}/N_TIA_OneStep.C
	${CMAKE_CURRENT_SOURCE_DIR}/N_TIA_StepErrorControl.C
	${CMAKE_CURRENT_SOURCE_DIR}/N_TIA_TIAParams.C
	${CMAKE_CURRENT_SOURCE_DIR}/N_TIA_TimeIntegrationMethods.C
	${CMAKE_CURRENT_SOURCE_DIR}/N_TIA_TwoLevelError.C
	${CMAKE_CURRENT_SOURCE_DIR}/N_TIA_RegisterTimeIntegrationMethods.C
	${CMAKE_CURRENT_SOURCE_DIR}/N_TIA_WorkingIntegrationMethod.C
)

add_to_xyce_library( NAME TimeIntegration
	FILES ${TimeIntegration_SOURCES}
)

# Add TimeIntegration Headers
set( TimeIntegration_Headers
	N_TIA_TwoLevelError.h
	N_TIA_fwd.h
)
install( FILES ${TimeIntegration_Headers}
	DESTINATION include
)


