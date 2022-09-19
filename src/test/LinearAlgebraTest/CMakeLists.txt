add_executable( testBlockLinearSystems testBlockLinearSystems.C )
target_link_libraries( testBlockLinearSystems XyceLib )

if( BUILD_TESTING )
     add_test( testBlockLinearSystems testBlockLinearSystems testBlockLinearSystems.C )
endif()
