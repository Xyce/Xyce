# -- process generated files for this library ---------------------------------

if ( Xyce_REACTION_PARSER )

  # TODO: Use BISON_TARGET and FLEX_TARGET here.

  # create source and header files
  add_custom_command ( OUTPUT ${Xyce_BINARY_DIR}/src/DeviceModelPKG/Core/N_DEV_ReactionParser.cxx
                              ${Xyce_BINARY_DIR}/src/DeviceModelPKG/Core/N_DEV_ReactionParser.hxx
                              ${Xyce_BINARY_DIR}/src/DeviceModelPKG/Core/location.hh
                              ${Xyce_BINARY_DIR}/src/DeviceModelPKG/Core/position.hh
                              ${Xyce_BINARY_DIR}/src/DeviceModelPKG/Core/stack.hh
                              ${Xyce_BINARY_DIR}/src/DeviceModelPKG/Core/N_DEV_ReactionLexer.cxx

                      COMMAND ${BISON_EXECUTABLE}  -o ${Xyce_BINARY_DIR}/src/DeviceModelPKG/Core/N_DEV_ReactionParser.cxx
                              ${Xyce_SOURCE_DIR}/src/DeviceModelPKG/Core/N_DEV_ReactionParser.yxx

                      COMMAND ${FLEX_EXECUTABLE} -+ -o${Xyce_BINARY_DIR}/src/DeviceModelPKG/Core/N_DEV_ReactionLexer.cxx
                              ${Xyce_SOURCE_DIR}/src/DeviceModelPKG/Core/N_DEV_ReactionLexer.l

                      DEPENDS ${Xyce_SOURCE_DIR}/src/DeviceModelPKG/Core/N_DEV_ReactionParser.yxx
                              ${Xyce_SOURCE_DIR}/src/DeviceModelPKG/Core/N_DEV_ReactionLexer.l )

  # force (re-)generation at build time
  add_custom_target ( reaction_parser_sources
                      DEPENDS ${Xyce_BINARY_DIR}/src/DeviceModelPKG/Core/N_DEV_ReactionParser.cxx
                              ${Xyce_BINARY_DIR}/src/DeviceModelPKG/Core/N_DEV_ReactionParser.hxx
                              ${Xyce_BINARY_DIR}/src/DeviceModelPKG/Core/location.hh
                              ${Xyce_BINARY_DIR}/src/DeviceModelPKG/Core/position.hh
                              ${Xyce_BINARY_DIR}/src/DeviceModelPKG/Core/stack.hh
                              ${Xyce_BINARY_DIR}/src/DeviceModelPKG/Core/N_DEV_ReactionLexer.cxx )

  # add to list of subpackage target dependencies
  set ( XYCE_TARGET_DEPENDENCIES ${XYCE_TARGET_DEPENDENCIES}
        "reaction_parser_sources"
        CACHE INTERNAL "X_T_D" )

  # save source location and header paths
  set ( REACTION_PARSER_SRC ${Xyce_BINARY_DIR}/src/DeviceModelPKG/Core/N_DEV_ReactionParser.cxx
                            ${Xyce_BINARY_DIR}/src/DeviceModelPKG/Core/N_DEV_ReactionLexer.cxx )

  set ( REACTION_PARSER_INC ${Xyce_BINARY_DIR}/src/DeviceModelPKG/Core )

  # pass list of generated files into global scope
  set ( XYCE_GENERATED_SOURCES ${XYCE_GENERATED_SOURCES}
        ${REACTION_PARSER_SRC}
        CACHE INTERNAL "X_G_S" )

endif ( Xyce_REACTION_PARSER )

# -- process files for this library ------------------------------------------

# class source list
set ( Core_SOURCES
      ${REACTION_PARSER_SRC}
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
      CACHE INTERNAL "X_S" )

# add header files to project
set ( XYCE_INCLUDE_DIRECTORIES ${XYCE_INCLUDE_DIRECTORIES}
      ${CMAKE_CURRENT_SOURCE_DIR}
      ${REACTION_PARSER_INC}
      CACHE INTERNAL "X_I" )

# -- build requested libraries -----------------------------------------------
set ( XYCE_PACKAGE_LIBRARIES ${XYCE_PACKAGE_LIBRARIES}
      Core
      CACHE INTERNAL "X_L" )



