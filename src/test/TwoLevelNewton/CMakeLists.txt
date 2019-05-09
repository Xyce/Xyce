
# -- build targets -----------------------------------------------------------

# create binary
add_executable( twoLevelNewtonLinearDC
                EXCLUDE_FROM_ALL
                ${CMAKE_CURRENT_SOURCE_DIR}/twoLevelNewtonLinearDC.C )
add_executable( twoLevelNewtonLinearTran
                EXCLUDE_FROM_ALL
                ${CMAKE_CURRENT_SOURCE_DIR}/twoLevelNewtonLinearTran.C )

target_include_directories( twoLevelNewtonLinearDC PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
)

target_include_directories( twoLevelNewtonLinearTran PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
)

# link against available Xyce library 
if ( Xyce_ENABLE_SHARED )
  target_link_libraries( twoLevelNewtonLinearDC lib_xyce_shared ${DAKOTA_OBJS} )
else ( Xyce_ENABLE_SHARED )
  target_link_libraries( twoLevelNewtonLinearDC lib_xyce_static ${DAKOTA_OBJS} )
endif ( Xyce_ENABLE_SHARED )
if ( Xyce_ENABLE_SHARED )
  target_link_libraries( twoLevelNewtonLinearTran lib_xyce_shared ${DAKOTA_OBJS} )
else ( Xyce_ENABLE_SHARED )
  target_link_libraries( twoLevelNewtonLinearTran lib_xyce_static ${DAKOTA_OBJS} )
endif ( Xyce_ENABLE_SHARED )

