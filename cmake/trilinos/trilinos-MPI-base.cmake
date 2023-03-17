include(${CMAKE_CURRENT_LIST_DIR}/trilinos-base.cmake)
set( TPL_ENABLE_MPI                          ON  CACHE BOOL "" )
set( Trilinos_ENABLE_Zoltan                  ON  CACHE BOOL "" )
set( Trilinos_ENABLE_Isorropia               ON  CACHE BOOL "" )

# The following must be ON to enable the ParMETIS matrix partitioner, which
# must already be installed on the system
#set(     Zoltan_ENABLE_ParMETIS              ON  CACHE BOOL "" )
#set( TPL_ENABLE_ParMETIS                     ON  CACHE BOOL "" )

# The following must be ON to enable the SuperLUDist solver, which must already
# be installed on the system
#set(     Amesos_ENABLE_SuperLUDist           ON  CACHE BOOL "" )
#set( TPL_ENABLE_SuperLUDist                  ON  CACHE BOOL "" )
