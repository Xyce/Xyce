add_executable( testBlockLinearSystems testBlockLinearSystems.C )
target_link_libraries( testBlockLinearSystems XyceLib )

if( BUILD_TESTING )
     add_test( testBlockLinearSystems testBlockLinearSystems )
     set_tests_properties( testBlockLinearSystems PROPERTIES
          REQUIRED_FILES "hb_matrix1.mm"
          REQUIRED_FILES "hb_rhs1.mm" )
     get_target_property(XyceLibDir XyceLib BINARY_DIR )
     set_tests_properties( testBlockLinearSystems PROPERTIES ENVIRONMENT_MODIFICATION "PATH=path_list_prepend:${XyceLibDir}")
     set_tests_properties( testBlockLinearSystems PROPERTIES LABELS "nightly;unittest" )
     file( COPY hb_matrix1.mm hb_rhs1.mm DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
endif()
