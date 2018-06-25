#
# FindTrilinos.cmake
#
# This CMake script will search for Trilinos and the various libraries and options.
#

set( Trilinos_INCLUDE_DIR "/Users/asgibso/Documents/Xyce/Trilinos_build/install/include" )

file( GLOB Trilinos_TEUCHOS_STATIC_LIB "/Users/asgibso/Documents/Xyce/Trilinos_build/install/lib/*.a" )
file( GLOB Trilinos_TEUCHOS_SHARED_LIB "/Users/asgibso/Documents/Xyce/Trilinos_build/install/lib/*.dylib" )

file( GLOB SuiteSparse_STATIC_LIB "/Users/asgibso/Documents/Xyce/SuiteSparse/lib/*.a" )
file( GLOB SuiteSparse_SHARED_LIB "/Users/asgibso/Documents/Xyce/SuiteSparse/lib/*.dylib" )

set( Trilinos_LIBRARIES
    ${Trilinos_TEUCHOS_STATIC_LIB}
    ${Trilinos_TEUCHOS_SHARED_LIB}
    ${SuiteSparse_STATIC_LIB}
    ${SuiteSparse_SHARED_LIB} )

