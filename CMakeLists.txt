cmake_minimum_required(VERSION 3.10 FATAL_ERROR) 

set(CMAKE_PREFIX_PATH ${Trilinos_DIRECTORY} ${CMAKE_PREFIX_PATH})

# Trilinos recommends it be found BEFORE project() is called, so the compilers
# and flags are properly set.
message("-- Looking for required Trilinos packages:\n"
        "        Amesos Epetra EpetraExt Ifpack NOX Teuchos Sacado\n"
        "        Triutils AztecOO Belos TrilinosCouplings")
find_package(Trilinos 12.12
     REQUIRED Amesos Epetra EpetraExt Ifpack NOX Teuchos Sacado Triutils
              AztecOO Belos TrilinosCouplings
     OPTIONAL_COMPONENTS Isorropia Zoltan ShyLU Basker Amesos2 Stokhos ROL)
message("-- Looking for required Trilinos packages - found")

# Make sure to use same compilers and flags as Trilinos
set( CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER} )
set( CMAKE_C_COMPILER ${Trilinos_C_COMPILER} )
set( CMAKE_Fortran_COMPILER ${Trilinos_Fortran_COMPILER} )

set( CMAKE_CXX_FLAGS  "${Trilinos_CXX_COMPILER_FLAGS} ${CMAKE_CXX_FLAGS}")
set( CMAKE_C_FLAGS  "${Trilinos_C_COMPILER_FLAGS} ${CMAKE_C_FLAGS}")
set( CMAKE_Fortran_FLAGS  "${Trilinos_Fortran_COMPILER_FLAGS} ${CMAKE_Fortran_FLAGS}")

message("CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")

project(Xyce VERSION 7.0.0 LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
include(cmake/config.cmake)

add_subdirectory(src) 
