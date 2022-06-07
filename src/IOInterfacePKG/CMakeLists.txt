
# required source subdirectories
add_subdirectory( Output )
# -- process files for this library ------------------------------------------

# class source list
target_sources(XyceLib PRIVATE
      N_IO_CmdParse.C
      N_IO_NetlistImportTool.C
      N_IO_DistToolFactory.C
      N_IO_DistToolDefault.C
      N_IO_DistToolDevBalanced.C
      N_IO_DistToolFlatRoundRobin.C
      N_IO_DistToolBase.C
      N_IO_ExtOutWrapper.C
      N_IO_MORAnalysisTool.C
      N_IO_OutputMgr.C
      N_IO_ParsingMgr.C
      N_IO_LoadManager.C
      N_IO_OutputResponse.C
      N_IO_OutputMacroResults.C
      N_IO_OutputMOR.C
      N_IO_OutputROM.C
      N_IO_OutputResults.C
      N_IO_OutputFileBase.C
      N_IO_OutputPrn.C
      N_IO_OutputCsd.C
      N_IO_OutputHDF5.C
      N_IO_PrintDeviceCount.C
      N_IO_InitialConditions.C
      N_IO_OutputterExternal.C
      N_IO_OutputterLocal.C
      N_IO_OutputterFrequencyCSV.C
      N_IO_OutputterFrequencyPrn.C
      N_IO_OutputterFrequencyProbe.C
      N_IO_OutputterFrequencyRaw.C
      N_IO_OutputterFrequencyRawASCII.C
      N_IO_OutputterFrequencyTecplot.C
      N_IO_OutputterHBTimeCSV.C
      N_IO_OutputterHBFreqCSV.C
      N_IO_OutputterHBICCSV.C
      N_IO_OutputterHBTimePrn.C
      N_IO_OutputterHBFreqPrn.C
      N_IO_OutputterHBICPrn.C
      N_IO_OutputterHBTimeTecplot.C
      N_IO_OutputterHBFreqTecplot.C
      N_IO_OutputterHBICTecplot.C
      N_IO_OutputterHDF5.C
      N_IO_OutputterHomotopy.C
      N_IO_OutputterHomotopyCSV.C
      N_IO_OutputterHomotopyPrn.C
      N_IO_OutputterHomotopyTecplot.C
      N_IO_OutputterMPDEPrn.C
      N_IO_OutputterMPDECSV.C
      N_IO_OutputterMPDETecplot.C
      N_IO_OutputterOverrideRaw.C
      N_IO_OutputterOverrideRawASCII.C
      N_IO_OutputterNoisePrn.C
      N_IO_OutputterNoiseCSV.C
      N_IO_OutputterNoiseTecplot.C
      N_IO_OutputterPCECSV.C
      N_IO_OutputterPCEPrn.C
      N_IO_OutputterPCETecplot.C
      N_IO_OutputterSensitivityPrn.C
      N_IO_OutputterSensitivityCSV.C
      N_IO_OutputterSensitivityTecplot.C
      N_IO_OutputterSensitivityDakota.C
      N_IO_OutputterSensitivityACPrn.C
      N_IO_OutputterSensitivityACCSV.C
      N_IO_OutputterSensitivityACTecplot.C
      N_IO_OutputterSParamTS1.C
      N_IO_OutputterSParamTS2.C
      N_IO_OutputterTimeCSV.C
      N_IO_OutputterTimePrn.C
      N_IO_OutputterTimeProbe.C
      N_IO_OutputterTimeRaw.C
      N_IO_OutputterTimeRawASCII.C
      N_IO_OutputterTimeTecplot.C
      N_IO_Probe.C
      N_IO_Tecplot.C
      N_IO_Outputter.C
      N_IO_OutputterAC.C
      N_IO_OutputterDC.C
      N_IO_OutputterHB.C
      N_IO_OutputterMPDE.C
      N_IO_OutputterRawOverride.C
      N_IO_OutputterNoise.C
      N_IO_OutputterPCE.C
      N_IO_OutputterSensitivity.C
      N_IO_OutputterSParam.C
      N_IO_OutputterTransient.C
      N_IO_RestartMgr.C
      N_IO_SpiceSeparatedFieldTool.C
      N_IO_CircuitMetadata.C
      N_IO_CircuitBlock.C
      N_IO_CircuitContext.C
      N_IO_DeviceBlock.C
      N_IO_OptionBlock.C
      N_IO_ParameterBlock.C
      N_IO_Report.C
      N_IO_FunctionBlock.C
      N_IO_RestartNode.C
      N_IO_Op.C
      N_IO_OpBuilders.C
      N_IO_ActiveOutput.C
      N_IO_mmio.C
      N_IO_PkgOptionsMgr.C
      N_IO_FourierMgr.C
      N_IO_FFTMgr.C
      N_IO_FFTAnalysis.C
      N_IO_ParsingHelpers.C
      N_IO_OutputterEmbeddedSamplingTecplot.C
      N_IO_OutputterEmbeddedSamplingPrn.C
      N_IO_OutputterEmbeddedSamplingCSV.C
      N_IO_OutputterEmbeddedSampling.C
      N_IO_WildcardSupport.C )

# add header files to project
target_include_directories(XyceLib PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)




# -- build requested libraries -----------------------------------------------
set ( XYCE_PACKAGE_LIBRARIES ${XYCE_PACKAGE_LIBRARIES}
      IOInterface
      CACHE INTERNAL "X_L" )

