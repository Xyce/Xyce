
# -- process files for this library ------------------------------------------

# class source list
target_sources(XyceStatic PUBLIC
      N_ANP_AC.C
      N_ANP_AnalysisBase.C
      N_ANP_AnalysisEvent.C
      N_ANP_AnalysisManager.C
      N_ANP_DCSweep.C
      N_ANP_Dakota.C
      N_ANP_EmbeddedSampling.C
      N_ANP_PCE.C
      N_ANP_HB.C
      N_ANP_InnerLevelSolve.C
      N_ANP_MOR.C
      N_ANP_MPDE.C
      N_ANP_MixedSignalManager.C
      N_ANP_NOISE.C
      N_ANP_Op.C
      N_ANP_OpBuilders.C
      N_ANP_OutputMgrAdapter.C
      N_ANP_ROL.C
      N_ANP_Report.C
      N_ANP_RegisterAnalysis.C
      N_ANP_Sampling.C
      N_ANP_SecondLevelManager.C
      N_ANP_Step.C
      N_ANP_SweepParam.C
      N_ANP_Transient.C
      N_ANP_UQSupport.C)

# add header files to project
target_include_directories (XyceStatic PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})




# -- build requested libraries ------------------------------------------------
set ( XYCE_PACKAGE_LIBRARIES ${XYCE_PACKAGE_LIBRARIES}
      Analysis
      CACHE INTERNAL "X_L" )

