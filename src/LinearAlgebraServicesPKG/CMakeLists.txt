
# -- process files for this library ------------------------------------------

# add option class sources
if ( Xyce_SHYLU )
  set ( SHYLU_SRC N_LAS_ShyLUSolver.C )
endif()

if ( Xyce_AMESOS2 )
  set ( AMESOS2_SRC N_LAS_Amesos2Solver.C )
endif()


set ( KSPARSE_SRC N_LAS_KSparseSolver.C
                  ksparse/Epetra_CrsKundertSparse.C
                  ksparse/alloc.c
#                  ksparse/driver.c  #this has main {} in it. Not suitable for inclusion into Xyce's shared library.
                  ksparse/spalloc.c
                  ksparse/spbuild.c
                  ksparse/spfactor.c
                  ksparse/spoutput.c
                  ksparse/spsolve.c
                  ksparse/sputils.c )
set ( KSPARSE_INC ksparse/ )


# class source list
target_sources(XyceLib PRIVATE
      ${SHYLU_SRC}
      ${AMESOS2_SRC}
      ${KSPARSE_SRC}
      N_LAS_BelosSolver.C
      N_LAS_SimpleSolver.C
      N_LAS_AmesosSolver.C
      N_LAS_AztecOOSolver.C
      N_LAS_BlockVector.C
      N_LAS_BlockSystemHelpers.C
      N_LAS_EpetraBlockMatrix.C
      N_LAS_EpetraBlockMultiVector.C
      N_LAS_EpetraGraph.C
      N_LAS_EpetraHelpers.C
      N_LAS_EpetraImporter.C
      N_LAS_EpetraMatrix.C
      N_LAS_EpetraProblem.C
      N_LAS_MatrixFreeEpetraOperator.C
      N_LAS_Builder.C
      N_LAS_HBDirectSolver.C
      N_LAS_ESDirectSolver.C
      N_LAS_PCEDirectSolver.C
      N_LAS_HBBlockJacobiPrecond.C
      N_LAS_HBBlockJacobiEpetraOperator.C
      N_LAS_ESBuilder.C
      N_LAS_PCEBuilder.C
      N_LAS_HBBuilder.C
      N_LAS_HBPrecondFactory.C
      N_LAS_HBSolverFactory.C
      N_LAS_ESSolverFactory.C
      N_LAS_PCESolverFactory.C
      N_LAS_IfpackPrecond.C
      N_LAS_FilteredMatrix.C
      N_LAS_FilteredMultiVector.C
      N_LAS_MultiVector.C
      N_LAS_MOROperators.C
      N_LAS_Problem.C
      N_LAS_QueryUtil.C
      N_LAS_TranSolverFactory.C
      N_LAS_System.C
      N_LAS_SystemHelpers.C
      N_LAS_TransformTool.C
      N_LAS_TrilinosPrecondFactory.C
      N_LAS_Vector.C )

# add header files to project
target_include_directories(XyceLib PUBLIC
      ${CMAKE_CURRENT_SOURCE_DIR}
      ${KSPARSE_INC} )


# -- build requested libraries -----------------------------------------------
set ( XYCE_PACKAGE_LIBRARIES ${XYCE_PACKAGE_LIBRARIES}
      LinearAlgebraServices
      CACHE INTERNAL "X_L" )
