
# -- process files for this library ------------------------------------------

# class source list
target_sources(XyceLib PRIVATE
      N_TOP_CktGraph.C
      N_TOP_CktGraphBasic.C
      N_TOP_CktNode.C
      N_TOP_CktNode_Dev.C
      N_TOP_CktNode_V.C
      N_TOP_Directory.C
      N_TOP_Indexor.C
      N_TOP_Node.C
      N_TOP_NodeDevBlock.C
      N_TOP_ParNode.C
      N_TOP_Topology.C
      N_TOP_SerialLSUtil.C
      N_TOP_ParLSUtil.C
      N_TOP_LSUtilFactory.C )
# add header files to project
target_include_directories(XyceLib PUBLIC
      ${CMAKE_CURRENT_SOURCE_DIR})




# -- build requested libraries -----------------------------------------------
set ( XYCE_PACKAGE_LIBRARIES ${XYCE_PACKAGE_LIBRARIES}
      TopoManager
      CACHE INTERNAL "X_L" )

