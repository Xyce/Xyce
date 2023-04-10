include(${CMAKE_CURRENT_LIST_DIR}/trilinos-MPI-base.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/trilinos-beta.cmake)
# The following enables the ShyLU hybrid/hybrid solver
set( Trilinos_ENABLE_ShyLU_DDCore            ON  CACHE BOOL "" )

# NOTES
# ShyLU requires OpenMP, which does not work on MacOS
# The OMP_NUM_THREADS environment variable must be set to a reasonable number
