
# -- process files for this library ------------------------------------------

# class source list
target_sources(XyceLib PRIVATE
      N_PDS_GlobalAccessor.C
      N_PDS_MPI.C
      N_PDS_Manager.C
      N_PDS_ParallelMachine.C
      N_PDS_EpetraMPIComm.C 
      N_PDS_EpetraSerialComm.C 
      N_PDS_EpetraHelpers.C 
      N_PDS_EpetraParMap.C
      N_PDS_Node.C  )

# add header files to project
target_include_directories(XyceLib PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)




# -- build requested libraries -----------------------------------------------
set ( XYCE_PACKAGE_LIBRARIES ${XYCE_PACKAGE_LIBRARIES}
      ParallelDist
      CACHE INTERNAL "X_L" )
