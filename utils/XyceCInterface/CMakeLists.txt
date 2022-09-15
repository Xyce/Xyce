
add_library( xycecinterface EXCLUDE_FROM_ALL )
target_sources( xycecinterface PRIVATE N_CIR_XyceCInterface.C PUBLIC N_CIR_XyceCInterface.h )
target_include_directories( xycecinterface PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> )
target_link_libraries( xycecinterface XyceLib )

# with xycecinterface no longer built by default it shouldn't be
# part of the installation.
#install( TARGETS xycecinterface DESTINATION lib )
#install( FILES N_CIR_XyceCInterface.h DESTINATION include )

