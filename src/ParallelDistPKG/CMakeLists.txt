#
# CMakeLists.txt for ParallelDistPKG
#

set ( ParallelDist_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/N_PDS_Comm.C 
    ${CMAKE_CURRENT_SOURCE_DIR}/N_PDS_GlobalAccessor.C 
    ${CMAKE_CURRENT_SOURCE_DIR}/N_PDS_MPI.C 
    ${CMAKE_CURRENT_SOURCE_DIR}/N_PDS_MPIComm.C 
    ${CMAKE_CURRENT_SOURCE_DIR}/N_PDS_Manager.C 
#    ${CMAKE_CURRENT_SOURCE_DIR}/N_PDS_ParComm.C 
    ${CMAKE_CURRENT_SOURCE_DIR}/N_PDS_ParMap.C 
    ${CMAKE_CURRENT_SOURCE_DIR}/N_PDS_ParallelMachine.C 
    ${CMAKE_CURRENT_SOURCE_DIR}/N_PDS_SerialComm.C 
#    ${CMAKE_CURRENT_SOURCE_DIR}/N_PDS_SerialParComm.C 
    ${CMAKE_CURRENT_SOURCE_DIR}/N_PDS_Node.C 
)
add_to_xyce_library( NAME ParallelDist
    FILES ${ParallelDist_SOURCES}
)

# Add ParallelDist headers
set( ParallelDist_Headers
    N_PDS_Manager.h
    N_PDS_fwd.h
    N_PDS_MPI.h
    N_PDS_ParallelMachine.h
    N_PDS_Serial.h
)
install( FILES ${ParallelDist_Headers}
    DESTINATION include
)

