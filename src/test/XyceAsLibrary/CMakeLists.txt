
if( BUILD_TESTING )

     add_executable( bug1574 EXCLUDE_FROM_ALL bug1574.C )
     target_link_libraries( bug1574 XyceLib )

     add_executable( testHarness1685 testHarness1685.C )
     target_link_libraries( testHarness1685 XyceLib )

endif()

