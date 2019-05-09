
# ----------------------------------------------------------------------------
# CMake Setup          
# ----------------------------------------------------------------------------


# DEBUG - OSX for Dakota 4 explicity link and path required for this object
# if ( CMAKE_HOST_APPLE AND Xyce_Dakota AND NOT Xyce_Dakota50 )
#  message ( WARNING "forcing CParam.o object to link..." )
#  set ( DAKOTA_OBJS ${Xyce_ARCHDIR}/dakota/lib/Cparam.o )
#endif ( CMAKE_HOST_APPLE AND Xyce_Dakota AND NOT Xyce_Dakota50 )


# -- process subdirectories --------------------------------------------------

add_subdirectory ( FFTInterface )
add_subdirectory ( LinearAlgebraTest )
add_subdirectory ( XyceAsLibrary )
add_subdirectory ( DeviceInterface )
add_subdirectory ( RandomTest )
add_subdirectory ( XygraTestHarnesses )
add_subdirectory ( GenExtTestHarnesses )
add_subdirectory ( MPITest )
add_subdirectory ( TwoLevelNewton )

