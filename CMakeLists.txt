cmake_minimum_required(VERSION 3.13 FATAL_ERROR)
cmake_policy(SET CMP0075 NEW)
cmake_policy(SET CMP0076 NEW)

option(Xyce_USE_SUPERBUILD "Use superbuild to install Trilinos?")
if(Xyce_USE_SUPERBUILD)
  include("XyceSuperBuild.cmake")
  return()
endif()

# Xyce should use the same compilers as Trilinos.  To ensure this, Trilinos
# should be loaded before "project" is called.  However, "project" probes the
# system and sets up some critical variables.  Therefore, make an initial call
# to "project" with "NONE" as the compiler set, then find Trilinos, then set up
# the compilers and call project again.
project(Xyce NONE)

message(STATUS "Looking for Trilinos\n"
  "   Required packages:\n"
  "        Amesos Epetra EpetraExt Ifpack NOX Teuchos Sacado\n"
  "        Triutils AztecOO Belos TrilinosCouplings\n"
  "   Optional packages:\n"
  "        Isorropia Zoltan ShyLU ShyLU_NodeBasker Amesos2 Stokhos ROL")
find_package(Trilinos CONFIG
  REQUIRED Amesos Epetra EpetraExt Ifpack NOX Teuchos Sacado Triutils
       AztecOO Belos TrilinosCouplings
  OPTIONAL_COMPONENTS Isorropia Zoltan ShyLU ShyLU_NodeBasker Amesos2 Stokhos ROL)
message(STATUS "Looking for Trilinos - found")
if(Trilinos_VERSION VERSION_LESS "12.12")
  message(FATAL_ERROR
    "ERROR: Trilinos version ${Trilinos_VERSION} is less than the required minimum of 12.12. Install a version of Trilinos of 12.12 or greater.")
endif()

# If the user explicitly set the compilers, use those; otherwise use the same
# compilers as Trilinos.
if( DEFINED CMAKE_C_COMPILER )
  message( "Using the user-defined C compiler: ${CMAKE_C_COMPILER}" )
else()
  message( STATUS "Using the Trilinos C compiler: ${Trilinos_C_COMPILER}" )
  set( CMAKE_C_COMPILER ${Trilinos_C_COMPILER} )
endif()

if( DEFINED CMAKE_CXX_COMPILER )
  message( "Using the user-defined C++ compiler: ${CMAKE_CXX_COMPILER}" )
else()
  message( STATUS "Using the Trilinos C++ compiler: ${Trilinos_CXX_COMPILER}" )
  set( CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER} )
endif()

enable_language(C)
enable_language(CXX)

# Set a default build type if none was specified (using the fact that
# CMAKE_CONFIGURATION_TYPES is non-zero for mutli-configuration generators).
#
# This has to be before "project" is called, because, for compilers targeting
# the MSVC ABI on Windows, "project" sets the build type.  This behavior was
# changed in CMake version 3.15, but was not made the *default* behavior.
# Since we're compatible with CMake 3.13, we can't use the new behavior. When
# the minimum CMake version requirement is increased, see CMake policy CMP0091:
#    <https://cmake.org/cmake/help/latest/policy/CMP0091.html>
# (Originally, this block was placed right before the
# "include(cmake/config.cmake)" call, below.)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  message("Setting the build type to \"Release\", since it was not explicitly set.")
  # Set the possible values of build type for ccmake and cmake-gui.
  # If a new build type is ever added to Xyce, then this list should be updated.
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "None" "Debug" "Release"
     "RelWithDebInfo" "MinSizeRel")
endif()

project(Xyce VERSION 7.3.0 LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib")
option(BUILD_SHARED_LIBS "Build shared libraries?" OFF)

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
