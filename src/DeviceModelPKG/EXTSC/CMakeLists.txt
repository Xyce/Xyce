
# -- process files for this library ------------------------------------------

# class source list
target_sources(XyceStatic PUBLIC
      N_DEV_1D_Interface.C
      N_DEV_CharonInterface.C
      N_DEV_ExternCodeInterface.C
      N_DEV_ExternDevice.C
      N_DEV_RegisterExternalDevices.C
      N_DEV_XyceInterface.C )

# add header files to project
target_include_directories(XyceStatic PUBLIC
      ${CMAKE_CURRENT_SOURCE_DIR})




# -- build requested libraries -----------------------------------------------
set ( XYCE_PACKAGE_LIBRARIES ${XYCE_PACKAGE_LIBRARIES}
      EXT
      CACHE INTERNAL "X_L" )
