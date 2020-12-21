# -- process generated files for this library ---------------------------------
# this code block had to be moved into ../src for include_directories reasons. Vile Cmake2 reasons.
#if( Xyce_REACTION_PARSER )
#    BISON_TARGET( XyceReactionParser N_DEV_ReactionParser.yxx ${CMAKE_CURRENT_BINARY_DIR}/N_DEV_ReactionParser.cxx )
#    FLEX_TARGET( XyceReactionScanner "N_DEV_ReactionLexer.l"
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
target_sources(XyceLib PRIVATE
#      ${BISON_XyceReactionParser_OUTPUTS}
      N_DEV_RegisterDevices.C
      N_DEV_DeviceInstance.C
      N_DEV_DeviceBlock.C
      N_DEV_DeviceState.C
      N_DEV_DeviceModel.C
      N_DEV_DeviceMgr.C
      N_DEV_DeviceSensitivities.C
      N_DEV_DeviceMaster.C
      N_DEV_Message.C
      N_DEV_MatrixLoadData.C
      N_DEV_NumericalJacobian.C
      N_DEV_DeviceEntity.C
      N_DEV_CompositeParam.C
      N_DEV_Source.C
      N_DEV_SourceData.C
      N_DEV_RxnSet.C
      N_DEV_Region.C
      N_DEV_RegionData.C
      N_DEV_SpecieSource.C
      N_DEV_Specie.C
      N_DEV_DeviceSupport.C
      N_DEV_Param.C
      N_DEV_DeviceOptions.C
      N_DEV_ScalingVars.C
      N_DEV_SolverState.C
      N_DEV_Reaction.C
      N_DEV_ReactionNetwork.C
      N_DEV_RateConstantCalculators.C
      N_DEV_Pars.C
      N_DEV_LaTexDoc.C
      N_DEV_Op.C
      N_DEV_OpBuilders.C
      N_DEV_Dump.C
      N_DEV_InstanceName.C
      N_DEV_Units.C
      N_DEV_Configuration.C
      N_DEV_Algorithm.C
      N_DEV_Print.C)
# add header files to project
target_include_directories(XyceLib PUBLIC
      ${CMAKE_CURRENT_SOURCE_DIR}
      ${REACTION_PARSER_INC})

# -- build requested libraries -----------------------------------------------
set ( XYCE_PACKAGE_LIBRARIES ${XYCE_PACKAGE_LIBRARIES}
      Core
      CACHE INTERNAL "X_L" )



