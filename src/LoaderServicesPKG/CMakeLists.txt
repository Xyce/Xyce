#
# CMakeLists.txt for LoaderServicesPKG
#

add_to_xyce_library( NAME LoaderServices
	FILES
		${CMAKE_CURRENT_SOURCE_DIR}/N_LOA_CktLoader.C 
		${CMAKE_CURRENT_SOURCE_DIR}/N_LOA_HBLoader.C 
		${CMAKE_CURRENT_SOURCE_DIR}/N_LOA_Loader.C 
		${CMAKE_CURRENT_SOURCE_DIR}/N_LOA_NonlinearEquationLoader.C
)

# Add headers from LoaderServicesPKG
set( LoaderServices_Headers
	N_LOA_fwd.h
	N_LOA_Loader.h
)
install(FILES ${LoaderServices_Headers}
	DESTINATION include
)

# -- process files for this library ------------------------------------------

# class source list
#set ( LoaderServices_SOURCES  
#      ${CMAKE_CURRENT_SOURCE_DIR}/N_LOA_CktLoader.C 
#      ${CMAKE_CURRENT_SOURCE_DIR}/N_LOA_HBLoader.C 
#      ${CMAKE_CURRENT_SOURCE_DIR}/N_LOA_Loader.C 
#      ${CMAKE_CURRENT_SOURCE_DIR}/N_LOA_NonlinearEquationLoader.C
#      CACHE INTERNAL "X_S" )


# add header files to project
#set ( XYCE_INCLUDE_DIRECTORIES ${XYCE_INCLUDE_DIRECTORIES} 
#      ${CMAKE_CURRENT_SOURCE_DIR}
#      CACHE INTERNAL "X_I" )




# -- build requested libraries -----------------------------------------------
#set ( XYCE_PACKAGE_LIBRARIES ${XYCE_PACKAGE_LIBRARIES} 
#      LoaderServices
#      CACHE INTERNAL "X_L" )

