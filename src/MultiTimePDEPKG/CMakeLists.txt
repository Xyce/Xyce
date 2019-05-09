
# -- process files for this library ------------------------------------------

# class source list
set ( MPDE_SOURCES  
      ${CMAKE_CURRENT_SOURCE_DIR}/N_MPDE_Builder.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_MPDE_Discretization.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_MPDE_Loader.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_MPDE_Manager.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_MPDE_SawtoothLoader.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_MPDE_WarpedPhaseCondition.C 
      CACHE INTERNAL "X_S" )

# add header files to project
set ( XYCE_INCLUDE_DIRECTORIES ${XYCE_INCLUDE_DIRECTORIES} 
      ${CMAKE_CURRENT_SOURCE_DIR}
      CACHE INTERNAL "X_I" )




# -- build requested libraries -----------------------------------------------
set ( XYCE_PACKAGE_LIBRARIES ${XYCE_PACKAGE_LIBRARIES} 
      MPDE 
      CACHE INTERNAL "X_L" )

## add_dependencies ( MPDE LinearAlgebraServices )
