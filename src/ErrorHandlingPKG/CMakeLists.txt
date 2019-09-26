
# -- process files for this library ------------------------------------------

# class source list
set ( ErrorHandling_SOURCES 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_ERH_ErrorMgr.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_ERH_Message.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_ERH_Messenger.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_ERH_Progress.C 
      CACHE INTERNAL "X_S" )

set (XYCE_SOURCES ${XYCE_SOURCES} ${ErrorHandling_SOURCES} PARENT_SCOPE)
# add header files to project
set ( XYCE_INCLUDE_DIRECTORIES ${XYCE_INCLUDE_DIRECTORIES} 
      ${CMAKE_CURRENT_SOURCE_DIR}
      PARENT_SCOPE )




# -- build requested libraries -----------------------------------------------
set ( XYCE_PACKAGE_LIBRARIES ${XYCE_PACKAGE_LIBRARIES} 
      ErrorHandling 
      CACHE INTERNAL "X_L" )

