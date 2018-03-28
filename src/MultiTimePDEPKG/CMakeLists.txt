#
# CMakeLists.txt for MultiTimePDE
#

set( MultiTimePDE_Src
	${CMAKE_CURRENT_SOURCE_DIR}/N_MPDE_Builder.C 
	${CMAKE_CURRENT_SOURCE_DIR}/N_MPDE_Discretization.C 
	${CMAKE_CURRENT_SOURCE_DIR}/N_MPDE_Loader.C 
	${CMAKE_CURRENT_SOURCE_DIR}/N_MPDE_Manager.C 
	${CMAKE_CURRENT_SOURCE_DIR}/N_MPDE_SawtoothLoader.C 
	${CMAKE_CURRENT_SOURCE_DIR}/N_MPDE_WarpedPhaseCondition.C 
)
add_to_xyce_library( NAME MultiTimePDE
	FILES ${MultiTimePDE_Src}
)
