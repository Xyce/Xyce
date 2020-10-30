
# -- process files for this library ------------------------------------------

# class source list
target_sources(XyceStatic PUBLIC
      N_LOA_CktLoader.C
      N_LOA_HBLoader.C
      N_LOA_Loader.C
      N_LOA_NonlinearEquationLoader.C
      N_LOA_ESLoader.C
      N_LOA_PCELoader.C )
# add header files to project
target_include_directories(XyceStatic PUBLIC
      ${CMAKE_CURRENT_SOURCE_DIR})




# -- build requested libraries -----------------------------------------------
set ( XYCE_PACKAGE_LIBRARIES ${XYCE_PACKAGE_LIBRARIES}
      LoaderServices
      CACHE INTERNAL "X_L" )

