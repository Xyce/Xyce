include(${CMAKE_CURRENT_LIST_DIR}/trilinos-base.cmake)
# The following enables ShyLU-Basker (assumes METIS is already installed)
set( Trilinos_ENABLE_OpenMP                  ON  CACHE BOOL "" )
set( Trilinos_ENABLE_ShyLU_NodeBasker        ON  CACHE BOOL "" )
set( TPL_ENABLE_METIS                        ON  CACHE BOOL "" )

# NOTES
# ShyLU-Basker requires OpenMP
# The OMP_NUM_THREADS environment variable must be set to a reasonable number
