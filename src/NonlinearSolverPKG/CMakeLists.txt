#
# CMakeLists.txt for NonlinearSolverPKG
#


# -- process files for this library ------------------------------------------

# add class files that should be included in this directory's library
target_sources(XyceLib PRIVATE
      N_NLS_ConductanceExtractor.C
      N_NLS_ConstraintBT.C
      N_NLS_DampedNewton.C
      N_NLS_LOCA_Group.C
      N_NLS_Manager.C
      N_NLS_NLParams.C
      N_NLS_NonLinearSolver.C
      N_NLS_NOX_AugmentLinSys_GStepping.C
      N_NLS_NOX_AugmentLinSys_IC.C
      N_NLS_NOX_AugmentLinSys_IC_Gmin.C
      N_NLS_NOX_AugmentLinSys_PseudoTransient.C
      N_NLS_NOX.C
      N_NLS_NOX_Group.C
      N_NLS_NOX_Interface.C
      N_NLS_NOX_ParameterSet.C
      N_NLS_NOX_PseudoTransientSolver.C
      N_NLS_NOX_PseudoTransientTest.C
      N_NLS_NOX_SharedSystem.C
      N_NLS_NOX_Vector.C
      N_NLS_NOX_XyceTests.C
      N_NLS_ParamMgr.C
      N_NLS_ReturnCodes.C
      N_NLS_Sensitivity.C
      N_NLS_SensitivityResiduals.C
      N_NLS_TwoLevelNewton.C )
# add header files to project
target_include_directories(XyceLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# -- build requested libraries -----------------------------------------------
set ( XYCE_PACKAGE_LIBRARIES ${XYCE_PACKAGE_LIBRARIES}
      NonlinearSolver
      CACHE INTERNAL "X_L" )
