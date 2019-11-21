
# -- process files for this library ------------------------------------------

# add option class sources
if ( Xyce_BELOS )
  set ( BELOS_SRC ${CMAKE_CURRENT_SOURCE_DIR}/N_LAS_BelosSolver.C )
endif ( Xyce_BELOS )

if ( Xyce_SHYLU )
  set ( SHYLU_SRC ${CMAKE_CURRENT_SOURCE_DIR}/N_LAS_ShyLUSolver.C )
endif ( Xyce_SHYLU )

if ( Xyce_AMESOS2 )
  set ( AMESOS2_SRC ${CMAKE_CURRENT_SOURCE_DIR}/N_LAS_Amesos2Solver.C )
endif ( Xyce_AMESOS2 )


set ( KSPARSE_SRC ${CMAKE_CURRENT_SOURCE_DIR}/N_LAS_KSparseSolver.C
                  ${CMAKE_CURRENT_SOURCE_DIR}/ksparse/Epetra_CrsKundertSparse.C
                  ${CMAKE_CURRENT_SOURCE_DIR}/ksparse/alloc.c
                  ${CMAKE_CURRENT_SOURCE_DIR}/ksparse/driver.c
                  ${CMAKE_CURRENT_SOURCE_DIR}/ksparse/spalloc.c
                  ${CMAKE_CURRENT_SOURCE_DIR}/ksparse/spbuild.c
                  ${CMAKE_CURRENT_SOURCE_DIR}/ksparse/spfactor.c
                  ${CMAKE_CURRENT_SOURCE_DIR}/ksparse/spoutput.c
                  ${CMAKE_CURRENT_SOURCE_DIR}/ksparse/spsolve.c
                  ${CMAKE_CURRENT_SOURCE_DIR}/ksparse/sputils.c )
set ( KSPARSE_INC ${CMAKE_CURRENT_SOURCE_DIR}/ksparse/ )


# class source list
set ( LinearAlgebraServices_SOURCES   
      ${BELOS_SRC}
      ${SHYLU_SRC}
      ${AMESOS2_SRC}
      ${KSPARSE_SRC}
      ${CMAKE_CURRENT_SOURCE_DIR}/N_LAS_SimpleSolver.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_LAS_AmesosSolver.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_LAS_AztecOOSolver.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_LAS_BlockMatrix.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_LAS_BlockVector.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_LAS_BlockMultiVector.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_LAS_BlockSystemHelpers.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_LAS_EpetraHelpers.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_LAS_Builder.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_LAS_HBDirectSolver.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_LAS_ESDirectSolver.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_LAS_PCEDirectSolver.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_LAS_HBBlockJacobiPrecond.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_LAS_HBBlockJacobiEpetraOperator.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_LAS_ESBuilder.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_LAS_PCEBuilder.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_LAS_HBBuilder.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_LAS_HBPrecondFactory.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_LAS_HBSolverFactory.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_LAS_ESSolverFactory.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_LAS_PCESolverFactory.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_LAS_IfpackPrecond.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_LAS_LAFactory.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_LAS_Matrix.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_LAS_FilteredMatrix.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_LAS_FilteredMultiVector.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_LAS_MultiVector.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_LAS_MOROperators.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_LAS_Problem.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_LAS_QueryUtil.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_LAS_TranSolverFactory.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_LAS_System.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_LAS_SystemHelpers.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_LAS_TransformTool.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_LAS_TrilinosPrecondFactory.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_LAS_Vector.C
      CACHE INTERNAL "X_S" )


# add header files to project
set ( XYCE_INCLUDE_DIRECTORIES ${XYCE_INCLUDE_DIRECTORIES} 
      ${CMAKE_CURRENT_SOURCE_DIR}
      ${KSPARSE_INC}
      CACHE INTERNAL "X_I" )


# -- build requested libraries -----------------------------------------------
set ( XYCE_PACKAGE_LIBRARIES ${XYCE_PACKAGE_LIBRARIES} 
      LinearAlgebraServices 
      CACHE INTERNAL "X_L" )
