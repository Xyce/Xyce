#
# CMakeLists.txt for Core Device Model stuff.
#
if( Xyce_REACTION_PARSER )
    BISON_TARGET( XyceReactionParser "${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_ReactionParser.yxx"
        "${CMAKE_CURRENT_BINARY_DIR}/N_DEV_ReactionParser.cxx" )
    FLEX_TARGET( XyceReactionScanner "${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_ReactionLexer.l"
        "${CMAKE_CURRENT_BINARY_DIR}/N_DEV_ReactionLexer.cxx")
    ADD_FLEX_BISON_DEPENDENCY( XyceReactionScanner XyceReactionParser )

    set( ReactionParserSrc
        ${BISON_XyceReactionParser_OUTPUTS}
        ${FLEX_XyceReactionScanner_OUTPUTS} )

    include_directories( ${CMAKE_CURRENT_BINARY_DIR} )
endif()

add_to_xyce_library( NAME Core
    FILES
        ${ReactionParserSrc}
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_RegisterDevices.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_DeviceInstance.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_DeviceBlock.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_DeviceState.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_DeviceModel.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_DeviceMgr.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_DeviceSensitivities.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_DeviceMaster.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_Message.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_MatrixLoadData.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_NumericalJacobian.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_DeviceEntity.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_CompositeParam.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_Source.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_SourceData.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_RxnSet.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_Region.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_RegionData.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_SpecieSource.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_Specie.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_DeviceSupport.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_Param.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_DeviceOptions.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_ScalingVars.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_SolverState.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_Reaction.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_ReactionNetwork.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_RateConstantCalculators.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_Pars.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_LaTexDoc.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_Op.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_OpBuilders.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_Dump.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_InstanceName.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_Units.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_Configuration.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_Algorithm.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_Print.C )

set( Core_Headers
    N_DEV_Configuration.h
    N_DEV_Const.h
    N_DEV_Device.h
    N_DEV_DeviceBlock.h
    N_DEV_DeviceEntity.h
    N_DEV_DeviceInstance.h
    N_DEV_DeviceMaster.h
    N_DEV_DeviceModel.h
    N_DEV_DeviceOptions.h
    N_DEV_DeviceSupport.h
    N_DEV_InstanceName.h
    N_DEV_ExternalSimulationData.h
    N_DEV_ExternData.h
    N_DEV_fwd.h
    N_DEV_Message.h
    N_DEV_Param.h
    N_DEV_Pars.h
    N_DEV_SolverState.h
    N_DEV_Units.h
    N_DEV_VectorComputeInterface.h )

install( FILES ${Core_Headers}
    DESTINATION include )

