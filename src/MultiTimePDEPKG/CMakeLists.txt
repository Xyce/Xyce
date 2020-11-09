
# -- process files for this library ------------------------------------------

# class source list
target_sources(XyceStatic PUBLIC
      N_MPDE_Builder.C
      N_MPDE_Discretization.C
      N_MPDE_Loader.C
      N_MPDE_Manager.C
      N_MPDE_SawtoothLoader.C
      N_MPDE_WarpedPhaseCondition.C )
# add header files to project
target_include_directories(XyceStatic PUBLIC
      ${CMAKE_CURRENT_SOURCE_DIR})




# -- build requested libraries -----------------------------------------------
set ( XYCE_PACKAGE_LIBRARIES ${XYCE_PACKAGE_LIBRARIES}
      MPDE
      CACHE INTERNAL "X_L" )

## add_dependencies ( MPDE LinearAlgebraServices )
