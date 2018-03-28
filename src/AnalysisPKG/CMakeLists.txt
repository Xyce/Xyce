#
# CMakeLists.txt for AnalysisPKG
#

add_to_xyce_library( NAME Analysis
	FILES
      ${CMAKE_CURRENT_SOURCE_DIR}/N_ANP_AC.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_ANP_AnalysisBase.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_ANP_AnalysisEvent.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_ANP_AnalysisManager.C  
      ${CMAKE_CURRENT_SOURCE_DIR}/N_ANP_DCSweep.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_ANP_Dakota.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_ANP_EmbeddedSampling.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_ANP_PCE.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_ANP_HB.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_ANP_InnerLevelSolve.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_ANP_MOR.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_ANP_MPDE.C    
      ${CMAKE_CURRENT_SOURCE_DIR}/N_ANP_MixedSignalManager.C  
      ${CMAKE_CURRENT_SOURCE_DIR}/N_ANP_NOISE.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_ANP_Op.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_ANP_OpBuilders.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_ANP_OutputMgrAdapter.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_ANP_ROL.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_ANP_Report.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_ANP_RegisterAnalysis.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_ANP_Sampling.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_ANP_SecondLevelManager.C  
      ${CMAKE_CURRENT_SOURCE_DIR}/N_ANP_Step.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_ANP_SweepParam.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_ANP_Transient.C 
)

set( AnalysisHeaders
	N_ANP_fwd.h
	N_ANP_NoiseData.h
)
install( FILES ${AnalysisHeaders}
	DESTINATION include
)

