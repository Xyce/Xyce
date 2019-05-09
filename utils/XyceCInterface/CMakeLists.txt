

# add this directory to include directories
#set ( XYCE_INCLUDE_DIRECTORIES ${XYCE_INCLUDE_DIRECTORIES} 
#      ${CMAKE_CURRENT_SOURCE_DIR}
#      CACHE INTERNAL "X_I" )


#
# only build the C interface library if 
# this is a shared library build
#
if ( Xyce_ENABLE_SHARED )

# add header directories to -Ilist 
include_directories ( ${XYCE_INCLUDE_DIRECTORIES} )

# add the Xyce C Interface library as a build target
add_library( lib_xyce_cinterface  SHARED N_CIR_XyceCInterface.C ) 

target_include_directories( lib_xyce_cinterface PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# make this target installable
install ( TARGETS lib_xyce_cinterface LIBRARY DESTINATION lib   RUNTIME DESTINATION bin)


# set the target name, use RPATH on MacOSX 
set_target_properties( lib_xyce_cinterface PROPERTIES 
                          MACOSX_RPATH True
                          OUTPUT_NAME xycecinterface 
                          CLEAN_DIRECT_OUTPUT ON )

# set dependencies for this target 
target_link_libraries ( lib_xyce_cinterface lib_xyce_shared ${ARCH_LIBRARIES} ${EXTERNAL_LIBRARIES} )

endif ( Xyce_ENABLE_SHARED )
