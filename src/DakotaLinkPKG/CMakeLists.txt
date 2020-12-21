
# -- process files for this library ------------------------------------------

# class source list
target_sources(XyceLib PRIVATE
      N_DAK_DakotaController.C
      N_DAK_DakotaInterface.C)

# add header files to project
target_include_directories(XyceLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )

# -- build requested libraries -----------------------------------------------
set ( XYCE_PACKAGE_LIBRARIES ${XYCE_PACKAGE_LIBRARIES}
      DakotaLink
      CACHE INTERNAL "X_L" )

