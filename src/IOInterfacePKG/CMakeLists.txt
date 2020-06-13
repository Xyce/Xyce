
# required source subdirectories
add_subdirectory( Output )
# -- process files for this library ------------------------------------------

# class source list
set ( IOInterface_SOURCES
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_CmdParse.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_NetlistImportTool.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_DistToolFactory.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_DistToolDefault.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_DistToolDevBalanced.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_DistToolFlatRoundRobin.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_DistToolBase.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_ExtOutWrapper.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_MORAnalysisTool.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputMgr.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_ParsingMgr.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_LoadManager.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputResponse.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputMacroResults.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputMOR.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputROM.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputResults.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputFileBase.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputPrn.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputCsd.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputHDF5.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_PrintDeviceCount.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_InitialConditions.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterExternal.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterLocal.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterFrequencyCSV.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterFrequencyPrn.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterFrequencyProbe.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterFrequencyRaw.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterFrequencyRawASCII.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterFrequencyTecplot.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterHBTimeCSV.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterHBFreqCSV.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterHBICCSV.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterHBTimePrn.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterHBFreqPrn.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterHBICPrn.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterHBTimeTecplot.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterHBFreqTecplot.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterHBICTecplot.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterHDF5.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterHomotopy.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterHomotopyCSV.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterHomotopyPrn.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterHomotopyTecplot.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterMPDEPrn.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterMPDECSV.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterMPDETecplot.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterOverrideRaw.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterOverrideRawASCII.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterNoisePrn.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterNoiseCSV.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterNoiseTecplot.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterPCECSV.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterPCEPrn.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterPCETecplot.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterSensitivityPrn.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterSensitivityCSV.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterSensitivityTecplot.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterSensitivityDakota.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterSensitivityACPrn.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterSensitivityACCSV.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterSensitivityACTecplot.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterSParamTS1.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterSParamTS2.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterTimeCSV.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterTimePrn.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterTimeProbe.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterTimeRaw.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterTimeRawASCII.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterTimeTecplot.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_Probe.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_Tecplot.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_Outputter.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterAC.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterDC.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterHB.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterMPDE.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterRawOverride.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterNoise.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterPCE.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterSensitivity.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterSParam.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterTransient.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_RestartMgr.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_SpiceSeparatedFieldTool.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_CircuitMetadata.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_CircuitBlock.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_CircuitContext.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_DeviceBlock.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OptionBlock.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_ParameterBlock.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_Report.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_FunctionBlock.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_RestartNode.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_Op.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OpBuilders.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_ActiveOutput.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_mmio.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_PkgOptionsMgr.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_FourierMgr.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_ParsingHelpers.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterEmbeddedSamplingTecplot.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterEmbeddedSamplingPrn.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterEmbeddedSamplingCSV.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_OutputterEmbeddedSampling.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_IO_WildcardSupport.C
      CACHE INTERNAL "X_S" )

set ( XYCE_SOURCES ${XYCE_SOURCES} ${IOInterface_SOURCES} PARENT_SCOPE )
# add header files to project
set ( XYCE_INCLUDE_DIRECTORIES ${XYCE_INCLUDE_DIRECTORIES}
      ${CMAKE_CURRENT_SOURCE_DIR}
      PARENT_SCOPE )




# -- build requested libraries -----------------------------------------------
set ( XYCE_PACKAGE_LIBRARIES ${XYCE_PACKAGE_LIBRARIES}
      IOInterface
      CACHE INTERNAL "X_L" )

