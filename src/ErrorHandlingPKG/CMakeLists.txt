
# -- process files for this library ------------------------------------------

# class source list
target_sources(XyceLib PRIVATE
      N_ERH_ErrorMgr.C
      N_ERH_Message.C
      N_ERH_Messenger.C
      N_ERH_Progress.C )

# add header files to project
target_include_directories(XyceLib PUBLIC
      ${CMAKE_CURRENT_SOURCE_DIR})




# -- build requested libraries -----------------------------------------------
set ( XYCE_PACKAGE_LIBRARIES ${XYCE_PACKAGE_LIBRARIES}
      ErrorHandling
      CACHE INTERNAL "X_L" )

