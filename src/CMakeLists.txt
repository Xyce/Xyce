set( XYCE_CONFIG_FILEPATH "${CMAKE_CURRENT_BINARY_DIR}/Xyce_config.h" )
configure_file( "${CMAKE_CURRENT_SOURCE_DIR}/Xyce_config.h.cmake" "${XYCE_CONFIG_FILEPATH}" )


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
add_subdirectory(test)

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
#End bison and Flex non-sense. 

# ---- This is the declaration and linking of the desired
# artifacts ----------------

add_library(XyceStatic STATIC ${XYCE_SOURCES} ${BISON_XyceReactionParser_OUTPUTS} ${FLEX_XyceReactionScanner_OUTPUTS}) 
target_include_directories(XyceStatic PUBLIC ${XYCE_INCLUDE_DIRECTORIES}
     ${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS} ${FLEX_INCLUDE_DIRS}
     ${Boost_INCLUDE_DIRS} ${CURL_INCLUDE_DIRS})
target_link_libraries(XyceStatic PUBLIC ${FFT} ${CURL_LIBRARIES} ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES} ${CMAKE_DL_LIBS})

add_executable(Xyce Xyce.C)
target_include_directories (Xyce PRIVATE ${XYCE_INCLUDE_DIRECTORIES})
target_link_libraries(Xyce XyceStatic) 


# ---- And the installation of them. ----------------
# rename the library		# a word to the wise about Windows
				# A Static library named "xyce" will produce 
				# a file named xyce.lib.
				# A Shared library named "xyce" will produce
				# two files, one named xyce.dll and another 
				# named xyce.lib, which is not the same as the
				# static library. It is an interface library.
				# Interface libraries may be quite foriegn to
				# Unix programmers. For right now the only way
				# I can tell xyce.lib the static library apart
				# from xyce.lib the interface library is that
				# the interface library is smaller in size.
				# Which isn't a very clear distinction. 
				# MORAL: if you would like to build both a
				# static and shared library use names such as
				# xyce-static and xyce-shared to make things
				# explicit. 
set_target_properties ( XyceStatic PROPERTIES
			OUTPUT_NAME xyce
			CLEAN_DIRECT_OUTPUT ON)
					
install(TARGETS XyceStatic DESTINATION lib)
install(TARGETS Xyce DESTINATION bin)


