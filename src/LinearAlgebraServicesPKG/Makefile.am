
AM_CPPFLAGS = @Xyce_INCS@

noinst_LTLIBRARIES = libLinearAlgebraServices.la

if BUILD_BELOS_SUPPORT
belos_SOURCES = \
  N_LAS_BelosSolver.C \
  N_LAS_BelosSolver.h
endif

if BUILD_SHYLU_SUPPORT
shylu_SOURCES = \
  N_LAS_ShyLUSolver.C \
  N_LAS_ShyLUSolver.h
endif

if BUILD_AMESOS2_SUPPORT
amesos2_SOURCES = \
  N_LAS_Amesos2Solver.C \
  N_LAS_Amesos2Solver.h
endif

ksparse_SOURCES = \
  N_LAS_KSparseSolver.C \
  N_LAS_KSparseSolver.h 

SUBDIRS=ksparse

libLinearAlgebraServices_la_SOURCES = \
  $(belos_SOURCES) \
  $(ml_SOURCES) \
  $(shylu_SOURCES) \
  $(amesos2_SOURCES) \
  $(ksparse_SOURCES) \
  N_LAS_SimpleSolver.C \
  N_LAS_AmesosSolver.C \
  N_LAS_AztecOOSolver.C \
  N_LAS_Problem.C \
  N_LAS_EpetraProblem.C \
  N_LAS_MultiVector.C \
  N_LAS_Matrix.C \
  N_LAS_FilteredMatrix.C \
  N_LAS_FilteredMultiVector.C \
  N_LAS_MOROperators.C \
  N_LAS_QueryUtil.C \
  N_LAS_Graph.C \
  N_LAS_System.C \
  N_LAS_SystemHelpers.C \
  N_LAS_Builder.C \
  N_LAS_ESBuilder.C \
  N_LAS_PCEBuilder.C \
  N_LAS_HBBuilder.C \
  N_LAS_Vector.C \
  N_LAS_TransformTool.C \
  N_LAS_TranSolverFactory.C \
  N_LAS_BlockVector.C \
  N_LAS_BlockMultiVector.C \
  N_LAS_BlockMatrix.C \
  N_LAS_BlockSystemHelpers.C \
  N_LAS_EpetraHelpers.C \
  N_LAS_MatrixFreeEpetraOperator.C \
  N_LAS_IfpackPrecond.C \
  N_LAS_ESDirectSolver.C \
  N_LAS_ESSolverFactory.C \
  N_LAS_PCEDirectSolver.C \
  N_LAS_PCESolverFactory.C \
  N_LAS_HBDirectSolver.C \
  N_LAS_HBBlockJacobiPrecond.C \
  N_LAS_HBBlockJacobiEpetraOperator.C \
  N_LAS_HBPrecondFactory.C \
  N_LAS_HBSolverFactory.C \
  N_LAS_TrilinosPrecondFactory.C \
  N_LAS_fwd.h \
  N_LAS_SimpleSolver.h \
  N_LAS_AmesosSolver.h \
  N_LAS_AztecOOSolver.h \
  N_LAS_Builder.h \
  N_LAS_HBBuilder.h \
  N_LAS_ESBuilder.h \
  N_LAS_PCEBuilder.h \
  N_LAS_Solver.h \
  N_LAS_Problem.h \
  N_LAS_EpetraProblem.h \
  N_LAS_Matrix.h \
  N_LAS_FilteredMatrix.h \
  N_LAS_FilteredMultiVector.h \
  N_LAS_MOROperators.h \
  N_LAS_MultiVector.h \
  N_LAS_QueryUtil.h \
  N_LAS_Graph.h \
  N_LAS_System.h \
  N_LAS_SystemHelpers.h \
  N_LAS_Vector.h \
  N_LAS_TransformTool.h \
  N_LAS_SolverFactory.h \
  N_LAS_TranSolverFactory.h \
  N_LAS_BlockVector.h \
  N_LAS_BlockMultiVector.h \
  N_LAS_BlockMatrix.h \
  N_LAS_BlockSystemHelpers.h \
  N_LAS_EpetraHelpers.h \
  N_LAS_MatrixFreeEpetraOperator.h \
  N_LAS_Operator.h \
  N_LAS_Preconditioner.h \
  N_LAS_NoPrecond.h \
  N_LAS_IfpackPrecond.h \
  N_LAS_ESDirectSolver.h \
  N_LAS_ESSolverFactory.h \
  N_LAS_PCEDirectSolver.h \
  N_LAS_PCESolverFactory.h \
  N_LAS_HBDirectSolver.h \
  N_LAS_HBBlockMatrixEntry.h \
  N_LAS_HBBlockJacobiPrecond.h \
  N_LAS_HBBlockJacobiEpetraOperator.h \
  N_LAS_PrecondFactory.h \
  N_LAS_HBPrecondFactory.h \
  N_LAS_HBSolverFactory.h \
  N_LAS_TrilinosPrecondFactory.h \
  N_LAS_Util.h

libLinearAlgebraServices_la_LIBADD=ksparse/libksparse.la
