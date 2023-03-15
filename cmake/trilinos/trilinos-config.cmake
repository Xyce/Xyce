set( Trilinos_ENABLE_NOX                     ON  CACHE BOOL "" )
set(     NOX_ENABLE_LOCA                     ON  CACHE BOOL "" )
set( Trilinos_ENABLE_EpetraExt               ON  CACHE BOOL "" )
set(     EpetraExt_BUILD_BTF                 ON  CACHE BOOL "" )
set(     EpetraExt_BUILD_EXPERIMENTAL        ON  CACHE BOOL "" )
set(     EpetraExt_BUILD_GRAPH_REORDERINGS   ON  CACHE BOOL "" )
set( Trilinos_ENABLE_TrilinosCouplings       ON  CACHE BOOL "" )
set( Trilinos_ENABLE_Ifpack                  ON  CACHE BOOL "" )
set( Trilinos_ENABLE_AztecOO                 ON  CACHE BOOL "" )
set( Trilinos_ENABLE_Belos                   ON  CACHE BOOL "" )
set( Trilinos_ENABLE_Teuchos                 ON  CACHE BOOL "" )
set( Trilinos_ENABLE_Amesos                  ON  CACHE BOOL "" )
set(     Amesos_ENABLE_KLU                   ON  CACHE BOOL "" )
set( Trilinos_ENABLE_Sacado                  ON  CACHE BOOL "" )
set( Trilinos_ENABLE_Stokhos                 ON  CACHE BOOL "" )
set( Trilinos_ENABLE_ROL                     ON  CACHE BOOL "" )
set( Trilinos_ENABLE_Amesos2                 ON  CACHE BOOL "" )
set(     Amesos2_ENABLE_Basker               ON  CACHE BOOL "" )
set( Trilinos_ENABLE_COMPLEX_DOUBLE          ON  CACHE BOOL "" )
set( Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES   OFF CACHE BOOL "" )
set( TPL_ENABLE_AMD                          ON  CACHE BOOL "" )
set( TPL_ENABLE_BLAS                         ON  CACHE BOOL "" )
set( TPL_ENABLE_LAPACK                       ON  CACHE BOOL "" )
set( CMAKE_POSITION_INDEPENDENT_CODE        TRUE CACHE BOOL "" )

# Enabling this with the Trilinos 14 branch resulted in a failed compilation:
# set( TPL_ENABLE_SuperLU                      ON  CACHE BOOL "" )

# I don't know how to enable SuperLUDist (which requires MPI?)

# ShyLU requires MPI
# Is this the correct invocation?
#set( Trilinos_ENABLE_ShyLU                   ON  CACHE BOOL "" )

# ShyLU-Basker requires OpenMP, which won't work on MacOS
# Is this the correct invocation?
#set( Trilinos_ENABLE_OpenMP                  ON  CACHE BOOL "" )
#set( Trilinos_ENABLE_ShyLU_NodeBasker        ON  CACHE BOOL "" )
#set(     Amesos2_ENABLE_ShyLU_NodeBasker     ON  CACHE BOOL "" )

# How to enable PARDISO_MKL?

