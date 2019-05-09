
# -- process files for this library ------------------------------------------

# add class files that should be included in this directory's library
set ( NonlinearSolver_SOURCES  
      ${CMAKE_CURRENT_SOURCE_DIR}/N_NLS_ConductanceExtractor.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_NLS_ConstraintBT.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_NLS_DampedNewton.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_NLS_LOCA_Group.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_NLS_Manager.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_NLS_MatrixFreeEpetraOperator.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_NLS_NLParams.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_NLS_NonLinearSolver.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_NLS_NOX_AugmentLinSys_GStepping.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_NLS_NOX_AugmentLinSys_IC.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_NLS_NOX_AugmentLinSys_IC_Gmin.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_NLS_NOX_AugmentLinSys_OPStart.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_NLS_NOX_AugmentLinSys_PseudoTransient.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_NLS_NOX.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_NLS_NOX_Group.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_NLS_NOX_Interface.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_NLS_NOX_ParameterSet.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_NLS_NOX_PseudoTransientSolver.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_NLS_NOX_PseudoTransientTest.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_NLS_NOX_SharedSystem.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_NLS_NOX_Vector.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_NLS_NOX_XyceTests.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_NLS_ParamMgr.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_NLS_Sensitivity.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_NLS_SensitivityResiduals.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_NLS_TwoLevelNewton.C 
      CACHE INTERNAL "X_S" )

  
# add header files to project
set ( XYCE_INCLUDE_DIRECTORIES ${XYCE_INCLUDE_DIRECTORIES} 
      ${CMAKE_CURRENT_SOURCE_DIR}
      CACHE INTERNAL "X_I" )




# -- build requested libraries -----------------------------------------------
set ( XYCE_PACKAGE_LIBRARIES ${XYCE_PACKAGE_LIBRARIES} 
      NonlinearSolver 
      CACHE INTERNAL "X_L" )
