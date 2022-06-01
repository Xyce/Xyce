
# -- process files for this library ------------------------------------------

# class source list
target_sources(XyceLib PRIVATE
      # DO NOT list Xyce.C here.  It is part of the Xyce binary target.
      N_CIR_Xyce.C
      N_CIR_SecondLevelSimulator.C
      N_CIR_MixedSignalSimulator.C
      N_CIR_GenCouplingSimulator.C)
# add header files to project
target_include_directories(XyceLib PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

add_subdirectory(test)
# -- build requested libraries -----------------------------------------------
# libxyce is created at top level.  See ${Xyce_SOURCE_DIR}/src/CMakeLists.txt
