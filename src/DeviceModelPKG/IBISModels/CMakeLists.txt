
# -- process files for this library ------------------------------------------

# class source list
target_sources(XyceLib PRIVATE
      N_DEV_RegisterIBISDevices.C
      N_DEV_IBIS.C)
# add header files to project
target_include_directories(XyceLib PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> )

# -- build requested libraries -----------------------------------------------
set ( XYCE_PACKAGE_LIBRARIES ${XYCE_PACKAGE_LIBRARIES}
      IBISModels
      CACHE INTERNAL "X_L" )
