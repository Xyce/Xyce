set( XYCE_CONFIG_FILEPATH "${CMAKE_CURRENT_BINARY_DIR}/Xyce_config.h" )
configure_file( "Xyce_config.h.cmake" "${XYCE_CONFIG_FILEPATH}" )

####
# Bison and Flex modules don't work in just any directory. Under the covers,
# they are CMake2 style, where they use include_directories.
if( Xyce_REACTION_PARSER )
     bison_target( XyceReactionParser DeviceModelPKG/Core/N_DEV_ReactionParser.yxx
          ${CMAKE_CURRENT_BINARY_DIR}/DeviceModelPKG/Core/N_DEV_ReactionParser.cxx )
     flex_target( XyceReactionLexer "DeviceModelPKG/Core/N_DEV_ReactionLexer.l"
          "${CMAKE_CURRENT_BINARY_DIR}/DeviceModelPKG/Core/N_DEV_ReactionLexer.cxx"
          COMPILE_FLAGS "-+")

     add_flex_bison_dependency( XyceReactionLexer XyceReactionParser )

     include_directories( ${CMAKE_CURRENT_BINARY_DIR}/DeviceModelPKG/Core )
endif()

bison_target( XyceExpressionParser
     "UtilityPKG/ExpressionSrc/ExpressionParser.yxx"
     "${CMAKE_CURRENT_BINARY_DIR}/UtilityPKG/ExpressionSrc/ExpressionParser.cxx")
flex_target( XyceExpressionLexer "UtilityPKG/ExpressionSrc/ExpressionLexer.l"
     "${CMAKE_CURRENT_BINARY_DIR}/UtilityPKG/ExpressionSrc/ExpressionLexer.cxx"
     COMPILE_FLAGS "-i -+ --full --prefix=exp")

add_flex_bison_dependency( XyceExpressionLexer XyceExpressionParser  )

include_directories( ${CMAKE_CURRENT_BINARY_DIR}/UtilityPKG/ExpressionSrc )

#End bison and Flex non-sense.

# ---- This is the declaration and linking of the desired artifacts ----------------

add_library(XyceLib
 ${BISON_XyceExpressionParser_OUTPUTS}
 ${FLEX_XyceExpressionLexer_OUTPUTS}
 ${BISON_XyceReactionParser_OUTPUTS}
 ${FLEX_XyceReactionLexer_OUTPUTS})

###
# add_subdirectories here

add_subdirectory(CircuitPKG)
add_subdirectory(AnalysisPKG)
add_subdirectory(DeviceModelPKG)
add_subdirectory(UtilityPKG)
add_subdirectory(ErrorHandlingPKG)
add_subdirectory(IOInterfacePKG)
add_subdirectory(LinearAlgebraServicesPKG)
add_subdirectory(LoaderServicesPKG)
add_subdirectory(ParallelDistPKG)
add_subdirectory(TopoManagerPKG)
add_subdirectory(NonlinearSolverPKG)
add_subdirectory(TimeIntegrationPKG)
add_subdirectory(DakotaLinkPKG)
add_subdirectory(MultiTimePDEPKG)
add_subdirectory(test)

# For debugging:

# get_target_property(XYCE_SOURCES XyceLib SOURCES)
# message(STATUS "Xyce Sources ${XYCE_SOURCES}")
# get_target_property(XYCE_INCLUDE_DIRECTORIES XyceLib INCLUDE_DIRECTORIES)
# message(STATUS "Xyce Inc ${XYCE_INCLUDE_DIRECTORIES}")


target_include_directories(XyceLib PUBLIC ${CMAKE_CURRENT_BINARY_DIR}
     ${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS} ${FLEX_INCLUDE_DIRS}
     ${Boost_INCLUDE_DIRS} ${CURL_INCLUDE_DIRS})
target_link_libraries(XyceLib PUBLIC ${FFT} ${CURL_LIBRARIES}
     ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES} ${CMAKE_DL_LIBS})

# In the case where CURL has been built with CMake, then CURL_LIBRARIES and
# CURL_INCLUDE_DIRS are NOT set; instead, there is an import library defined.
# Check for that here:
if (Xyce_USE_CURL AND  NOT (DEFINED CURL_LIBRARIES))
  target_link_libraries(XyceLib PUBLIC CURL::libcurl)
endif()

add_executable(Xyce Xyce.C)
target_include_directories (Xyce PRIVATE ${XYCE_INCLUDE_DIRECTORIES})
target_link_libraries(Xyce XyceLib)

if(WIN32)
  if(DEFINED ENV{MKLROOT})

    if( CMAKE_SIZEOF_VOID_P EQUAL 8 )
      set( POINTER_SIZE  "64" )
    else()
      set( POINTER_SIZE "32" )
    endif()
    add_custom_command(
      TARGET Xyce POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different
              $ENV{MKLROOT}/../redist/intel${POINTER_SIZE}/mkl/mkl_sequential.dll
              ${CMAKE_CURRENT_BINARY_DIR})
  endif()
endif()


# ---- And the installation of them. ----------------
     # A word to the wise about Windows:
     # A Static library named "xyce" will produce a file named xyce.lib.
     # A Shared library named "xyce" will produce two files, one named xyce.dll
     # and another named xyce.lib, which is not the same as the static library.
     # It is an interface library.  Interface libraries may be quite foriegn to
     # Unix programmers. For right now the only way I can tell xyce.lib the
     # static library apart from xyce.lib the interface library is that the
     # interface library is smaller in size.  Which isn't a very clear
     # distinction.  MORAL: if you would like to build both a static and shared
     # library, use names such as xyce-static and xyce-shared to make things
     # explicit.
set_target_properties ( XyceLib PROPERTIES
			OUTPUT_NAME xyce
			CLEAN_DIRECT_OUTPUT ON)

install(TARGETS XyceLib DESTINATION lib)
install(TARGETS Xyce DESTINATION bin)

# Add header files needed by other codes to link against libxyce.*
#
# We need to do the following because Xyce was not designed to have a clean set
# of interface header files. It is possible could include *all* .h files using:
#    install(DIRECTORY src/ DESTINATION include/myproj
#            FILES_MATCHING PATTERN "*.h")
# However, that would result in something close to 600 header files being
# installed (over 500 of which are not needed).
install(FILES
    AnalysisPKG/N_ANP_fwd.h
    AnalysisPKG/N_ANP_NoiseData.h
    CircuitPKG/N_CIR_Xyce.h
    CircuitPKG/N_CIR_GenCouplingSimulator.h
    DeviceModelPKG/Core/N_DEV_Configuration.h
    DeviceModelPKG/Core/N_DEV_Const.h
    DeviceModelPKG/Core/N_DEV_Device.h
    DeviceModelPKG/Core/N_DEV_DeviceBlock.h
    DeviceModelPKG/Core/N_DEV_DeviceEntity.h
    DeviceModelPKG/Core/N_DEV_DeviceInstance.h
    DeviceModelPKG/Core/N_DEV_DeviceMaster.h
    DeviceModelPKG/Core/N_DEV_DeviceModel.h
    DeviceModelPKG/Core/N_DEV_DeviceOptions.h
    DeviceModelPKG/Core/N_DEV_DeviceSupport.h
    DeviceModelPKG/Core/N_DEV_InstanceName.h
    DeviceModelPKG/Core/N_DEV_ExternalSimulationData.h
    DeviceModelPKG/Core/N_DEV_ExternData.h
    DeviceModelPKG/Core/N_DEV_fwd.h
    DeviceModelPKG/Core/N_DEV_Message.h
    DeviceModelPKG/Core/N_DEV_Param.h
    DeviceModelPKG/Core/N_DEV_Pars.h
    DeviceModelPKG/Core/N_DEV_SolverState.h
    DeviceModelPKG/Core/N_DEV_Units.h
    DeviceModelPKG/Core/N_DEV_VectorComputeInterface.h
    DeviceModelPKG/OpenModels/N_DEV_ADC.h
    DeviceModelPKG/OpenModels/N_DEV_DAC.h
    DeviceModelPKG/OpenModels/N_DEV_Resistor.h
    DeviceModelPKG/OpenModels/N_DEV_Capacitor.h
    DeviceModelPKG/OpenModels/N_DEV_BJT.h
    DeviceModelPKG/OpenModels/N_DEV_Diode.h
    DeviceModelPKG/OpenModels/N_DEV_JFET.h
    DeviceModelPKG/OpenModels/N_DEV_MOSFET1.h
    ErrorHandlingPKG/N_ERH_fwd.h
    ErrorHandlingPKG/N_ERH_Message.h
    IOInterfacePKG/N_IO_CmdParse.h
    IOInterfacePKG/N_IO_ExtOutInterface.h
    IOInterfacePKG/N_IO_HangingResistor.h
    IOInterfacePKG/N_IO_OutputTypes.h
    IOInterfacePKG/N_IO_fwd.h
    LinearAlgebraServicesPKG/N_LAS_fwd.h
    LinearAlgebraServicesPKG/N_LAS_Vector.h
    LinearAlgebraServicesPKG/N_LAS_Matrix.h
    LinearAlgebraServicesPKG/N_LAS_MultiVector.h
    LoaderServicesPKG/N_LOA_fwd.h
    LoaderServicesPKG/N_LOA_Loader.h
    NonlinearSolverPKG/N_NLS_fwd.h
    NonlinearSolverPKG/N_NLS_TwoLevelEnum.h
    ParallelDistPKG/N_PDS_Manager.h
    ParallelDistPKG/N_PDS_fwd.h
    ParallelDistPKG/N_PDS_MPI.h
    ParallelDistPKG/N_PDS_ParallelMachine.h
    ParallelDistPKG/N_PDS_Serial.h
    TimeIntegrationPKG/N_TIA_TwoLevelError.h
    TimeIntegrationPKG/N_TIA_fwd.h
    TopoManagerPKG/N_TOP_fwd.h
    UtilityPKG/N_UTL_BreakPoint.h
    UtilityPKG/N_UTL_Diagnostic.h
    UtilityPKG/N_UTL_Expression.h
    UtilityPKG/N_UTL_ExpressionSymbolTable.h
    UtilityPKG/N_UTL_ExpressionData.h
    UtilityPKG/N_UTL_fwd.h
    UtilityPKG/N_UTL_FeatureTest.h
    UtilityPKG/N_UTL_FormatTime.h
    UtilityPKG/N_UTL_Interface_Enum_Types.h
    UtilityPKG/N_UTL_JSON.h
    UtilityPKG/N_UTL_MachDepParams.h
    UtilityPKG/N_UTL_Marshal.h
    UtilityPKG/N_UTL_Misc.h
    UtilityPKG/N_UTL_Math.h
    UtilityPKG/N_UTL_NameLevelKey.h
    UtilityPKG/N_UTL_NetlistLocation.h
    UtilityPKG/N_UTL_NoCase.h
    UtilityPKG/N_UTL_NodeSymbols.h
    UtilityPKG/N_UTL_Op.h
    UtilityPKG/N_UTL_Pack.h
    UtilityPKG/N_UTL_Param.h
    UtilityPKG/N_UTL_ReportHandler.h
    UtilityPKG/N_UTL_Stats.h
    UtilityPKG/N_UTL_StatMetricTraits.h
    UtilityPKG/N_UTL_TypeIndex.h
    UtilityPKG/N_UTL_Xyce.h
    ${CMAKE_BINARY_DIR}/src/Xyce_config.h
  DESTINATION include)
