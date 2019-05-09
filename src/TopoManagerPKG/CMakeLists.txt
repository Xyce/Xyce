
# -- process files for this library ------------------------------------------

# class source list
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
      CACHE INTERNAL "X_S" )
  

# add header files to project
set ( XYCE_INCLUDE_DIRECTORIES ${XYCE_INCLUDE_DIRECTORIES} 
      ${CMAKE_CURRENT_SOURCE_DIR}
      CACHE INTERNAL "X_I" )




# -- build requested libraries -----------------------------------------------
set ( XYCE_PACKAGE_LIBRARIES ${XYCE_PACKAGE_LIBRARIES} 
      TopoManager 
      CACHE INTERNAL "X_L" )

