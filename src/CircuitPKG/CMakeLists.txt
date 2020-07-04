
# -- process files for this library ------------------------------------------

# class source list
set ( Circuit_SOURCES
      # DO NOT list Xyce.C here.  It is part of the Xyce binary target. 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_CIR_Xyce.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_CIR_SecondLevelSimulator.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_CIR_MixedSignalSimulator.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_CIR_GenCouplingSimulator.C
      CACHE INTERNAL "X_S" )

set ( XYCE_SOURCES ${XYCE_SOURCES} ${Circuit_SOURCES} PARENT_SCOPE)
# add header files to project
set ( XYCE_INCLUDE_DIRECTORIES ${XYCE_INCLUDE_DIRECTORIES} 
      ${CMAKE_CURRENT_SOURCE_DIR}
      PARENT_SCOPE )

# -- build requested libraries -----------------------------------------------
# libxyce is created at top level.  See ${Xyce_SOURCE_DIR}/src/CMakeLists.txt 
