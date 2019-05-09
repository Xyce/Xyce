
# -- process files for this library ------------------------------------------

# class source list
set ( ErrorHandling_SOURCES 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_ERH_ErrorMgr.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_ERH_Message.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_ERH_Messenger.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/N_ERH_Progress.C 
      CACHE INTERNAL "X_S" )


# add header files to project
set ( XYCE_INCLUDE_DIRECTORIES ${XYCE_INCLUDE_DIRECTORIES} 
      ${CMAKE_CURRENT_SOURCE_DIR}
      CACHE INTERNAL "X_I" )




# -- build requested libraries -----------------------------------------------
set ( XYCE_PACKAGE_LIBRARIES ${XYCE_PACKAGE_LIBRARIES} 
      ErrorHandling 
      CACHE INTERNAL "X_L" )

