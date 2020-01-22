cmake_minimum_required(VERSION 3.10 FATAL_ERROR) 

# Xyce should use the same compilers as Trilinos. Therefore, Trilinos must be
# loaded before "project" is called.
message(STATUS "Looking for required Trilinos packages:\n"
     "        Amesos Epetra EpetraExt Ifpack NOX Teuchos Sacado\n"
     "        Triutils AztecOO Belos TrilinosCouplings\n"
     "   Also looking for the optional packages:\n"
     "        Isorropia Zoltan ShyLU Basker Amesos2 Stokhos ROL")
find_package(Trilinos 12.12
     REQUIRED Amesos Epetra EpetraExt Ifpack NOX Teuchos Sacado Triutils
              AztecOO Belos TrilinosCouplings
     OPTIONAL_COMPONENTS Isorropia Zoltan ShyLU Basker Amesos2 Stokhos ROL)
message(STATUS "Looking for required Trilinos packages - found")

# Set the compilers, so "project" finds the right ones
set( CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER} )
set( CMAKE_C_COMPILER ${Trilinos_C_COMPILER} )

project(Xyce VERSION 8.0.0 LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
include(cmake/config.cmake)

add_subdirectory(src) 

# add in unit tests subdirectory if it exists
if ( EXISTS ${Xyce_SOURCE_DIR}/Xyce_UnitTests)
  add_subdirectory (Xyce_UnitTests)
endif ( EXISTS ${Xyce_SOURCE_DIR}/Xyce_UnitTests)

include(CPack)
#See the CPackConfig.cmake in Xyce/cmake for settings
