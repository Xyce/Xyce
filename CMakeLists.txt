cmake_minimum_required(VERSION 3.13 FATAL_ERROR)
cmake_policy(SET CMP0076 NEW)

option(Xyce_USE_SUPERBUILD "Use superbuild to install Trilinos?"  )
if(Xyce_USE_SUPERBUILD)
  include("XyceSuperBuild.cmake")
  return()
endif()

# Xyce should use the same compilers as Trilinos. Therefore, Trilinos must be
# loaded before "project" is called.
message(STATUS "Looking for required Trilinos packages:\n"
  "        Amesos Epetra EpetraExt Ifpack NOX Teuchos Sacado\n"
  "        Triutils AztecOO Belos TrilinosCouplings\n"
  "   Also looking for the optional packages:\n"
  "        Isorropia Zoltan ShyLU Basker Amesos2 Stokhos ROL")
find_package(Trilinos
  REQUIRED Amesos Epetra EpetraExt Ifpack NOX Teuchos Sacado Triutils
       AztecOO Belos TrilinosCouplings
  OPTIONAL_COMPONENTS Isorropia Zoltan ShyLU Basker Amesos2 Stokhos ROL)
message(STATUS "Looking for required Trilinos packages - found")
if(Trilinos_VERSION VERSION_LESS "12.12")
  message(FATAL_ERROR
    "ERROR: Trilinos version ${Trilinos_VERSION} is less than the required minimum of 12.12.1. Install a version of Trilinos of 12.12.1 or greater.")
endif()

# Set the compilers, so "project" finds the right ones
# if the trilinos ones don't exist fall back on
# environment defaults
if( DEFINED CMAKE_CXX_COMPILER )
  message( "Using user-defined CXX ${CMAKE_CXX_COMPILER}" )
else()
  message( STATUS "Using Trilinos CXX compiler" )
  set( CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER} )
endif()

if( DEFINED CMAKE_C_COMPILER )
  message( "Using user-defined CC ${CMAKE_C_COMPILER}" )
else()
  message( STATUS "Using Trilinos C compiler" )
  set( CMAKE_C_COMPILER ${Trilinos_C_COMPILER} )
endif()

project(Xyce VERSION 8.0.0 LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib")
option(BUILD_SHARED_LIBS "Build shared libraries?" OFF)

# Don't force a build type in mutli-configuration platforms
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build.")
  message("Setting the build type to \"Release\", since it was not explicitly set.")
endif()

include(cmake/config.cmake)

if(Xyce_PLUGIN_SUPPORT)
  if(NOT BUILD_SHARED_LIBS)
    message(WARNING "Forcing BUILD_SHARED_LIBS to ON for plugin support.")
    set(BUILD_SHARED_LIBS ON CACHE BOOL "Build shared libraries?" FORCE)
  endif()
endif()

add_subdirectory(src)

# add in unit tests subdirectory if it exists
if ( EXISTS ${Xyce_SOURCE_DIR}/Xyce_UnitTests)
  add_subdirectory (Xyce_UnitTests)
endif ( EXISTS ${Xyce_SOURCE_DIR}/Xyce_UnitTests)

# Enable Xyce plugin capability
if(Xyce_PLUGIN_SUPPORT)
  find_program(ADMS_XML admsXml)

  if(ADMS_XML)
    configure_file(
      utils/buildxyceplugin.cmake.in
      buildxyceplugin.sh
      @ONLY)
    configure_file(
      utils/XycePluginProject.cmake.in
      XycePluginProject.cmake
      @ONLY)

    install(
      PROGRAMS
        ${ADMS_XML}
        ${CMAKE_BINARY_DIR}/buildxyceplugin.sh
      DESTINATION bin)

    install(
      FILES
        utils/ADMS/adms.implicit.xml
        utils/ADMS/xyceVersion_nosac.xml
        utils/ADMS/xyceBasicTemplates_nosac.xml
        utils/ADMS/xyceAnalogFunction_nosac.xml
        utils/ADMS/xyceHeaderFile_nosac.xml
        utils/ADMS/xyceImplementationFile_nosac.xml
        ${CMAKE_BINARY_DIR}/XycePluginProject.cmake
      DESTINATION share)

    install(
      DIRECTORY utils/ADMS/examples/toys
      DESTINATION share/examples)

    message(STATUS "Plugin compatibility enabled")
  else()
    message(WARNING "ADMS not found.  Disabling the plugin capability.")
    set(Xyce_PLUGIN_SUPPORT OFF CACHE BOOL "Install Xyce with plugin compatibility" FORCE)
  endif()
else()
  message(STATUS "Plugin compatibility not enabled")
endif()

# if no generator type is specified, use RPM
set(GEN_TYPE "RPM" CACHE STRING "What generator to use.")

include(CPack)
#See the CPackConfig.cmake in Xyce/cmake for settings
