set( XYCE_CONFIG_FILEPATH "${CMAKE_CURRENT_BINARY_DIR}/Xyce_config.h" )
configure_file( "${CMAKE_CURRENT_SOURCE_DIR}/Xyce_config.h.cmake" "${XYCE_CONFIG_FILEPATH}" )
install( FILES "${XYCE_CONFIG_FILEPATH}"
        DESTINATION include )


set( XYCE_SOURCES "Xyce.C")
set( XYCE_INCLUDE_DIRECTORIES ${CMAKE_CURRENT_BINARY_DIR})
###
# add_subdirectories here
#add_subdirectory()

add_subdirectory(CircuitPKG)
add_subdirectory(AnalysisPKG) 
add_subdirectory(DeviceModelPKG)
add_subdirectory(UtilityPKG)
add_subdirectory(ErrorHandlingPKG)
add_subdirectory(IOInterfacePKG)
add_subdirectory(LinearAlgebraServicesPKG)
add_subdirectory(LoaderServicesPKG)
add_subdirectory(ParallelDistPKG)
add_subdirectory(TopoManagerPKG)
add_subdirectory(NonlinearSolverPKG)
add_subdirectory(TimeIntegrationPKG)
add_subdirectory(DakotaLinkPKG)
add_subdirectory(MultiTimePDEPKG)
# For debugging:
#    message(STATUS "Xyce Source ${XYCE_SOURCES}")
#    message(STATUS "Xyce Inc ${XYCE_INCLUDE_DIRECTORIES}")

####
# Bison and Flex modules don't work in any directory. Under the covers, they
# are CMake2 style, where they use include_directories :(
if( Xyce_REACTION_PARSER )
    BISON_TARGET( XyceReactionParser ${CMAKE_CURRENT_SOURCE_DIR}/DeviceModelPKG/Core/N_DEV_ReactionParser.yxx
         ${CMAKE_CURRENT_BINARY_DIR}/DeviceModelPKG/Core/N_DEV_ReactionParser.cxx )
    FLEX_TARGET( XyceReactionScanner "${CMAKE_CURRENT_SOURCE_DIR}/DeviceModelPKG/Core/N_DEV_ReactionLexer.l"
        "${CMAKE_CURRENT_BINARY_DIR}/DeviceModelPKG/Core/N_DEV_ReactionLexer.cxx"
	COMPILE_FLAGS "-+"
		)
    ADD_FLEX_BISON_DEPENDENCY( XyceReactionScanner XyceReactionParser )

    set( REACTION_PARSER_SRC
        ${BISON_XyceReactionParser_OUTPUT_SOURCE}
        ${FLEX_XyceReactionScanner_OUTPUT_SOURCE} )

    include_directories( ${CMAKE_CURRENT_BINARY_DIR}/DeviceModelPKG/Core )
endif()

add_library(XyceStatic STATIC ${XYCE_SOURCES} ${BISON_XyceReactionParser_OUTPUTS} ${FLEX_XyceReactionScanner_OUTPUTS}) 
target_include_directories(XyceStatic PUBLIC ${Boost_INCLUDE_DIRS}
${XYCE_INCLUDE_DIRECTORIES} ${Trilinos_INCLUDE_DIRS}
${Trilinos_TPL_INCLUDE_DIRS} ${FLEX_INCLUDE_DIR})
target_link_libraries(XyceStatic PUBLIC FFTW::FFTW ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})

add_executable(Xyce Xyce.C)
target_include_directories (Xyce PRIVATE ${XYCE_INCLUDE_DIRECTORIES})
target_link_libraries(Xyce XyceStatic) 
