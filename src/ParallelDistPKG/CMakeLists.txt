
# -- process files for this library ------------------------------------------

# class source list
set ( ParallelDist_SOURCES
      ${CMAKE_CURRENT_SOURCE_DIR}/N_PDS_Comm.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_PDS_GlobalAccessor.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_PDS_MPI.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_PDS_MPIComm.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_PDS_Manager.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_PDS_ParMap.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_PDS_ParallelMachine.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_PDS_SerialComm.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_PDS_Node.C 
      CACHE INTERNAL "X_S" )

# add header files to project
set ( XYCE_INCLUDE_DIRECTORIES ${XYCE_INCLUDE_DIRECTORIES} 
      ${CMAKE_CURRENT_SOURCE_DIR}
      CACHE INTERNAL "X_I" )




# -- build requested libraries -----------------------------------------------
set ( XYCE_PACKAGE_LIBRARIES ${XYCE_PACKAGE_LIBRARIES} 
      ParallelDist 
      CACHE INTERNAL "X_L" )
