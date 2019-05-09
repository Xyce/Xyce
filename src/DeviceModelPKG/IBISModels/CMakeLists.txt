
# -- process files for this library ------------------------------------------

# class source list
set ( IBISModels_SOURCES  
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_RegisterIBISDevices.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_DEV_IBIS.C
      CACHE INTERNAL "X_S" )
  

# add header files to project
set ( XYCE_INCLUDE_DIRECTORIES ${XYCE_INCLUDE_DIRECTORIES} 
      ${CMAKE_CURRENT_SOURCE_DIR}
      CACHE INTERNAL "X_I" )

# -- build requested libraries -----------------------------------------------
set ( XYCE_PACKAGE_LIBRARIES ${XYCE_PACKAGE_LIBRARIES} 
      IBISModels
      CACHE INTERNAL "X_L" )
