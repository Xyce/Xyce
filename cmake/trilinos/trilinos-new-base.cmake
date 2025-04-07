# create a list of variables to turn ON for trilinos cmake configuration
set(Trilinos_ENABLE_Amesos ON CACHE BOOL "")
set(Trilinos_ENABLE_Amesos2 ON CACHE BOOL "")
set(Trilinos_ENABLE_AztecOO ON CACHE BOOL "")
set(Trilinos_ENABLE_Belos ON CACHE BOOL "")
set(Trilinos_ENABLE_Epetra ON CACHE BOOL "")
set(Trilinos_ENABLE_EpetraExt ON CACHE BOOL "")
set(Trilinos_ENABLE_Ifpack ON CACHE BOOL "")
set(Trilinos_ENABLE_Isorropia ON CACHE BOOL "")
set(Trilinos_ENABLE_Kokkos ON CACHE BOOL "")
set(Trilinos_ENABLE_NOX ON CACHE BOOL "")
set(Trilinos_ENABLE_Sacado ON CACHE BOOL "")
set(Trilinos_ENABLE_Stokhos ON CACHE BOOL "")
set(Trilinos_ENABLE_Teuchos ON CACHE BOOL "")
set(Trilinos_ENABLE_TrilinosCouplings ON CACHE BOOL "")
set(Trilinos_ENABLE_Zoltan ON CACHE BOOL "")

set(Trilinos_ENABLE_CXX11 ON CACHE BOOL "")
set(Trilinos_ENABLE_COMPLEX_DOUBLE ON CACHE BOOL "")

set(TPL_ENABLE_AMD ON CACHE BOOL "")
set(TPL_ENABLE_BLAS ON CACHE BOOL "")
set(TPL_ENABLE_LAPACK ON CACHE BOOL "")

set(NOX_ENABLE_AztecOO ON CACHE BOOL "")
set(NOX_ENABLE_Epetra ON CACHE BOOL "")
set(NOX_ENABLE_Ifpack ON CACHE BOOL "")

set(EpetraExt_BUILD_BTF ON CACHE BOOL "")
set(EpetraExt_BUILD_EXPERIMENTAL ON CACHE BOOL "")
set(EpetraExt_BUILD_GRAPH_REORDERINGS ON CACHE BOOL "")

set(TPL_ENABLE_AMD ON CACHE BOOL "")
set(TPL_ENABLE_BLAS ON CACHE BOOL "")
set(TPL_ENABLE_LAPACK ON CACHE BOOL "")

set(Amesos_ENABLE_KLU ON CACHE BOOL "")
set(Amesos2_ENABLE_Basker ON CACHE BOOL "")
set(Amesos2_ENABLE_KLU2 ON CACHE BOOL "")

# if the NKL library is loaded set things appropriately
if(DEFINED ENV{MKLROOT})
  set(TPL_ENABLE_MKL ON CACHE BOOL "")

  set(MKL_LIBRARY_DIRS $ENV{MKLROOT}/lib CACHE STRING "")
  set(MKL_INCLUDE_DIRS $ENV{MKLROOT}/include CACHE STRING "")

  set(BLAS_LIBRARY_NAMES mkl_rt CACHE STRING "")
  set(BLAS_LIBRARY_DIRS $ENV{MKLROOT}/lib CACHE STRING "")

  set(LAPACK_LIBRARY_NAMES mkl_rt CACHE STRING "")
  set(LAPACK_LIBRARY_DIRS $ENV{MKLROOT}/lib CACHE STRING "")
endif()

# create a list of variables to turn OFF for trilinos cmake configuration
set(TPL_ENABLE_DLlib OFF CACHE BOOL "")
set(Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES OFF CACHE BOOL "")
