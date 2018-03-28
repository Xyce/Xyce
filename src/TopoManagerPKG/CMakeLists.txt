#
# CMakeLists.txt for TopoManagerPKG
#

set ( TopoManager_SOURCES 
	${CMAKE_CURRENT_SOURCE_DIR}/N_TOP_CktGraph.C 
	${CMAKE_CURRENT_SOURCE_DIR}/N_TOP_CktGraphBasic.C 
	${CMAKE_CURRENT_SOURCE_DIR}/N_TOP_CktNode.C 
	${CMAKE_CURRENT_SOURCE_DIR}/N_TOP_CktNode_Dev.C 
	${CMAKE_CURRENT_SOURCE_DIR}/N_TOP_CktNode_V.C 
	${CMAKE_CURRENT_SOURCE_DIR}/N_TOP_Directory.C 
	${CMAKE_CURRENT_SOURCE_DIR}/N_TOP_Indexor.C 
	${CMAKE_CURRENT_SOURCE_DIR}/N_TOP_Node.C 
 	${CMAKE_CURRENT_SOURCE_DIR}/N_TOP_NodeDevBlock.C 
 	${CMAKE_CURRENT_SOURCE_DIR}/N_TOP_ParNode.C 
 	${CMAKE_CURRENT_SOURCE_DIR}/N_TOP_Topology.C 
 	${CMAKE_CURRENT_SOURCE_DIR}/N_TOP_SerialLSUtil.C 
 	${CMAKE_CURRENT_SOURCE_DIR}/N_TOP_ParLSUtil.C 
 	${CMAKE_CURRENT_SOURCE_DIR}/N_TOP_LSUtilFactory.C 
)
add_to_xyce_library( NAME TopoManager
	FILES ${TopoManager_SOURCES}
)

# Add TopoManager headers
set( TopoManager_Headers
	N_TOP_fwd.h
)
install( FILES ${TopoManager_Headers}
	DESTINATION include
)


