# -- process generated files for this library ---------------------------------
# this code block had to be moved into ../src for include_directories reasons. Vile Cmake2 reasons. 
#if( Xyce_REACTION_PARSER )
#    BISON_TARGET( XyceReactionParser N_DEV_ReactionParser.yxx ${CMAKE_CURRENT_BINARY_DIR}/N_DEV_ReactionParser.cxx )
#    FLEX_TARGET( XyceReactionScanner "${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_ReactionLexer.l"
#        "${CMAKE_CURRENT_BINARY_DIR}/N_DEV_ReactionLexer.cxx")
#    ADD_FLEX_BISON_DEPENDENCY( XyceReactionScanner XyceReactionParser )
#
#    set( REACTION_PARSER_SRC
#        ${BISON_XyceReactionParser_OUTPUT_SOURCE}
#        ${FLEX_XyceReactionScanner_OUTPUT_SOURCE} )
#
#    include_directories( ${CMAKE_CURRENT_BINARY_DIR} )
#endif()

# -- process files for this library ------------------------------------------

# class source list
set ( Core_SOURCES
#      ${BISON_XyceReactionParser_OUTPUTS}
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
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_Print.C
      )
set (XYCE_SOURCES ${XYCE_SOURCES} ${Core_SOURCES} PARENT_SCOPE)
# add header files to project
set ( XYCE_INCLUDE_DIRECTORIES ${XYCE_INCLUDE_DIRECTORIES}
      ${CMAKE_CURRENT_SOURCE_DIR}
      ${REACTION_PARSER_INC}
      PARENT_SCOPE )

# -- build requested libraries -----------------------------------------------
set ( XYCE_PACKAGE_LIBRARIES ${XYCE_PACKAGE_LIBRARIES}
      Core
      CACHE INTERNAL "X_L" )



