#
# CMakeLists.txt for NonlinearSolverPKG
#

# This instructs N_NLS_NOX_Vector to use NOX::SizeType, which is in turn
# defined through Trilinos. I'm not sure why it is still an option, but
# let's define it here for now.
add_definitions( -DXyce_NOX_SIZETYPE )

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
     ${CMAKE_CURRENT_SOURCE_DIR}/N_NLS_TwoLevelNewton.C )

add_to_xyce_library( NAME NonlinearSolver
    FILES ${NonlinearSolver_SOURCES} )

# Add headers for NonlinearSolver
set( NonlinearSolver_Headers
    N_NLS_fwd.h
    N_NLS_TwoLevelEnum.h )

install( FILES ${NonlinearSolver_Headers}
    DESTINATION include )

