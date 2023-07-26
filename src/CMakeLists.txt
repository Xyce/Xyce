set( XYCE_CONFIG_FILEPATH "${CMAKE_CURRENT_BINARY_DIR}/Xyce_config.h" )
configure_file( "Xyce_config.h.cmake" "${XYCE_CONFIG_FILEPATH}" )

# check for the updated flex/bison installation includ directory and
# add it to the include list to avoid the system default older
# version(s)
find_path(XYCE_FLEX_INCLUDE_PATH NAMES FlexLexer.h
  HINTS /projects/xyce/flexbison/include
  HINTS ${FLEX_DIR}/include)

####
# Bison and Flex modules don't work in just any directory. Under the covers,
# they are CMake2 style, where they use include_directories.
if( Xyce_REACTION_PARSER )
     bison_target( XyceReactionParser DeviceModelPKG/Core/N_DEV_ReactionParser.yxx
          ${CMAKE_CURRENT_BINARY_DIR}/DeviceModelPKG/Core/N_DEV_ReactionParser.cxx )
     flex_target( XyceReactionLexer "DeviceModelPKG/Core/N_DEV_ReactionLexer.l"
          "${CMAKE_CURRENT_BINARY_DIR}/DeviceModelPKG/Core/N_DEV_ReactionLexer.cxx"
          COMPILE_FLAGS "-+")

     add_flex_bison_dependency( XyceReactionLexer XyceReactionParser )

     include_directories( ${CMAKE_CURRENT_BINARY_DIR}/DeviceModelPKG/Core )
endif()

bison_target( XyceExpressionParser
     "UtilityPKG/ExpressionSrc/ExpressionParser.yxx"
     "${CMAKE_CURRENT_BINARY_DIR}/UtilityPKG/ExpressionSrc/ExpressionParser.cxx")
flex_target( XyceExpressionLexer "UtilityPKG/ExpressionSrc/ExpressionLexer.l"
     "${CMAKE_CURRENT_BINARY_DIR}/UtilityPKG/ExpressionSrc/ExpressionLexer.cxx"
     COMPILE_FLAGS "-i -+ --full --prefix=exp")

add_flex_bison_dependency( XyceExpressionLexer XyceExpressionParser  )

if(XYCE_FLEX_INCLUDE_PATH)
  include_directories( ${XYCE_FLEX_INCLUDE_PATH} ${CMAKE_CURRENT_BINARY_DIR}/UtilityPKG/ExpressionSrc )
else()
  include_directories( ${CMAKE_CURRENT_BINARY_DIR}/UtilityPKG/ExpressionSrc )
endif()

#End bison and Flex non-sense.

# ---- This is the declaration and linking of the desired artifacts ----------------

add_library(XyceLib
 ${BISON_XyceExpressionParser_OUTPUTS}
 ${FLEX_XyceExpressionLexer_OUTPUTS}
 ${BISON_XyceReactionParser_OUTPUTS}
 ${FLEX_XyceReactionLexer_OUTPUTS})

###
# add_subdirectories here

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

# get_target_property(XYCE_SOURCES XyceLib SOURCES)
# message(STATUS "Xyce Sources ${XYCE_SOURCES}")
# get_target_property(XYCE_INCLUDE_DIRECTORIES XyceLib INCLUDE_DIRECTORIES)
# message(STATUS "Xyce Inc ${XYCE_INCLUDE_DIRECTORIES}")


target_include_directories(XyceLib PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
     ${FLEX_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS} ${CURL_INCLUDE_DIRS})

#
# Xyce-Dakota linkage
#
if( Xyce_Dakota )
  target_include_directories(XyceLib PUBLIC ${Dakota_INCLUDE_DIRS} )
 
  target_link_libraries(XyceLib PUBLIC ${Dakota_LIBRARIES} )
  target_link_libraries(XyceLib PUBLIC ${Dakota_TPL_LIBRARIES} )
endif()

# 
# if we are using FFTW we need to add it's components rather than it as a 
# target so that tools that link to Xyce externally will have the right
# dependencies
#
if( Xyce_USE_FFTW )
  get_property(FFTWLibraries TARGET ${FFT} PROPERTY INTERFACE_LINK_LIBRARIES )
  target_link_libraries(XyceLib PUBLIC ${FFTWLibraries} )
  
  get_property(FFTWIncludeDir TARGET ${FFT} PROPERTY INTERFACE_INCLUDE_DIRECTORIES )
  target_include_directories(XyceLib PUBLIC ${FFTWIncludeDir} )
endif()

target_link_libraries(XyceLib PUBLIC ${CURL_LIBRARIES} ${CMAKE_DL_LIBS})
target_link_libraries(XyceLib PUBLIC ${TriBITS_prefix}::all_selected_libs)

if(OpenMP_FOUND)
     target_link_libraries(XyceLib PUBLIC OpenMP::OpenMP_CXX)
endif()

# See https://cmake.org/cmake/help/latest/module/FindCURL.html
# The "FindCURL" provided by CMake defines CURL_LIBRARIES and CURL_INCLUDE_DIRS
# (used above).  If curl was built using its CMake buildsystem, though, then it
# provides its own CURLConfig.cmake file, which defines an import library.
# Handle that case here:
if (Xyce_USE_CURL AND  NOT (DEFINED CURL_LIBRARIES))
  target_link_libraries(XyceLib PUBLIC CURL::libcurl)
endif()

add_executable(Xyce Xyce.C)
#target_include_directories (Xyce PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}> ${XYCE_INCLUDE_DIRECTORIES})
target_include_directories (XyceLib PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}> $<INSTALL_INTERFACE:include>)
target_link_libraries(Xyce XyceLib)

if(WIN32)
  if(DEFINED ENV{MKLROOT})

    file(TO_CMAKE_PATH "$ENV{MKLROOT}" MKLROOT)

    add_custom_command(
      TARGET Xyce POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different
              ${MKLROOT}/redist/intel64/mkl_sequential.2.dll
              ${CMAKE_CURRENT_BINARY_DIR})
  endif()
endif()


# ---- And the installation of them. ----------------
# Keep this in mind if Xyce needs to support building both static and shared libraries
# in the same build tree.
#
# According to <https://gitlab.kitware.com/cmake/community/-/wikis/FAQ>, the
# above is not needed: "One more detail is needed with CMake 2.6.x and lower
# (but not CMake 2.8 or higher). If you are building your shared and static
# libraries in the same directory, you will also need the following to keep
# your shared and static libraries from clobbering each other during the
# build."
#
# From <https://cmake.org/cmake/help/cmake2.6docs.html>:
# "When a library is built, CMake by default generates code to remove any
# existing library using all possible names. This is needed to support
# libraries that switch between STATIC and SHARED by a user option. However
# when using OUTPUT_NAME to build a static and shared library of the same name
# using different logical target names the two targets will remove each other's
# files. This can be prevented by setting the CLEAN_DIRECT_OUTPUT property to 1."

if(BUILD_SHARED_LIBS AND CMAKE_HOST_WIN32 )
  # if we are building shared libraries under windows then a byproduct will
  # be an interface library with a suffix ".lib" that is used in linking 
  # the executable to the dynamic library ".dll"  This interface library is 
  # needed and should not be confused with a static library which also ends
  # in ".lib" under Windwos.  So to avoid confusion, make the shared library
  # name xyce-share.dll and xyce-share.lib for the interface library under 
  # Windows. 
  #
  # A word to the wise about Windows:
  # A Static library named "xyce" will produce a file named xyce.lib.
  # A Shared library named "xyce" will produce two files, one named xyce.dll
  # and another named xyce.lib, which is not the same as the static library.
  # It is an interface library.  Interface libraries may be quite foriegn to
  # Unix programmers. For right now the only way I can tell xyce.lib the
  # static library apart from xyce.lib the interface library is that the
  # interface library is smaller in size.  Which isn't a very clear
  # distinction.  MORAL: if you would like to build both a static and shared
  # library, use names such as xyce-static and xyce-shared to make things
  # explicit.
  
  set_target_properties ( XyceLib PROPERTIES
	 		OUTPUT_NAME xyce-share
			CLEAN_DIRECT_OUTPUT ON)
  install(TARGETS XyceLib DESTINATION bin EXPORT XyceLibTarget)
else()
  set_target_properties ( XyceLib PROPERTIES
			OUTPUT_NAME xyce
			CLEAN_DIRECT_OUTPUT ON)
  install(TARGETS XyceLib DESTINATION lib EXPORT XyceLibTarget)
endif()

install(TARGETS Xyce DESTINATION bin EXPORT XyceTarget)

# Add header files needed by other codes to link against libxyce.*
#
# We need to do the following because Xyce was not designed to have a clean set
# of interface header files. It is possible could include *all* .h files using:
#    install(DIRECTORY src/ DESTINATION include/myproj
#            FILES_MATCHING PATTERN "*.h")
# However, that would result in something close to 600 header files being
# installed (over 500 of which are not needed).
install(FILES
    AnalysisPKG/N_ANP_fwd.h
    AnalysisPKG/N_ANP_NoiseData.h
    CircuitPKG/N_CIR_Xyce.h
    CircuitPKG/N_CIR_GenCouplingSimulator.h
    DeviceModelPKG/Core/N_DEV_Configuration.h
    DeviceModelPKG/Core/N_DEV_Const.h
    DeviceModelPKG/Core/N_DEV_Depend.h
    DeviceModelPKG/Core/N_DEV_Device.h
    DeviceModelPKG/Core/N_DEV_DeviceBlock.h
    DeviceModelPKG/Core/N_DEV_DeviceEntity.h
    DeviceModelPKG/Core/N_DEV_DeviceInstance.h
    DeviceModelPKG/Core/N_DEV_DeviceMaster.h
    DeviceModelPKG/Core/N_DEV_DeviceModel.h
    DeviceModelPKG/Core/N_DEV_DeviceOptions.h
    DeviceModelPKG/Core/N_DEV_DeviceSupport.h
    DeviceModelPKG/Core/N_DEV_InstanceName.h
    DeviceModelPKG/Core/N_DEV_ExternalSimulationData.h
    DeviceModelPKG/Core/N_DEV_ExternData.h
    DeviceModelPKG/Core/N_DEV_fwd.h
    DeviceModelPKG/Core/N_DEV_Message.h
    DeviceModelPKG/Core/N_DEV_Param.h
    DeviceModelPKG/Core/N_DEV_Pars.h
    DeviceModelPKG/Core/N_DEV_SolverState.h
    DeviceModelPKG/Core/N_DEV_Units.h
    DeviceModelPKG/Core/N_DEV_VectorComputeInterface.h
    DeviceModelPKG/OpenModels/N_DEV_ADC.h
    DeviceModelPKG/OpenModels/N_DEV_DAC.h
    DeviceModelPKG/OpenModels/N_DEV_Resistor.h
    DeviceModelPKG/OpenModels/N_DEV_Capacitor.h
    DeviceModelPKG/OpenModels/N_DEV_BJT.h
    DeviceModelPKG/OpenModels/N_DEV_Diode.h
    DeviceModelPKG/OpenModels/N_DEV_JFET.h
    DeviceModelPKG/OpenModels/N_DEV_MOSFET1.h
    ErrorHandlingPKG/N_ERH_fwd.h
    ErrorHandlingPKG/N_ERH_Message.h
    IOInterfacePKG/N_IO_CmdParse.h
    IOInterfacePKG/N_IO_ExtOutInterface.h
    IOInterfacePKG/N_IO_HangingResistor.h
    IOInterfacePKG/N_IO_OutputTypes.h
    IOInterfacePKG/N_IO_fwd.h
    LinearAlgebraServicesPKG/N_LAS_fwd.h
    LinearAlgebraServicesPKG/N_LAS_Vector.h
    LinearAlgebraServicesPKG/N_LAS_Matrix.h
    LinearAlgebraServicesPKG/N_LAS_MultiVector.h
    LoaderServicesPKG/N_LOA_fwd.h
    LoaderServicesPKG/N_LOA_Loader.h
    NonlinearSolverPKG/N_NLS_fwd.h
    NonlinearSolverPKG/N_NLS_TwoLevelEnum.h
    ParallelDistPKG/N_PDS_Manager.h
    ParallelDistPKG/N_PDS_fwd.h
    ParallelDistPKG/N_PDS_MPI.h
    ParallelDistPKG/N_PDS_ParallelMachine.h
    ParallelDistPKG/N_PDS_Serial.h
    TimeIntegrationPKG/N_TIA_TwoLevelError.h
    TimeIntegrationPKG/N_TIA_fwd.h
    TopoManagerPKG/N_TOP_fwd.h
    UtilityPKG/N_UTL_BreakPoint.h
    UtilityPKG/N_UTL_Diagnostic.h
    UtilityPKG/N_UTL_Expression.h
    UtilityPKG/N_UTL_ExpressionSymbolTable.h
    UtilityPKG/N_UTL_ExpressionData.h
    UtilityPKG/N_UTL_fwd.h
    UtilityPKG/N_UTL_FeatureTest.h
    UtilityPKG/N_UTL_FormatTime.h
    UtilityPKG/N_UTL_Interface_Enum_Types.h
    UtilityPKG/N_UTL_JSON.h
    UtilityPKG/N_UTL_MachDepParams.h
    UtilityPKG/N_UTL_Marshal.h
    UtilityPKG/N_UTL_Misc.h
    UtilityPKG/N_UTL_Math.h
    UtilityPKG/N_UTL_NameLevelKey.h
    UtilityPKG/N_UTL_NetlistLocation.h
    UtilityPKG/N_UTL_NoCase.h
    UtilityPKG/N_UTL_NodeSymbols.h
    UtilityPKG/N_UTL_Op.h
    UtilityPKG/N_UTL_Pack.h
    UtilityPKG/N_UTL_Param.h
    UtilityPKG/N_UTL_ReportHandler.h
    UtilityPKG/N_UTL_Stats.h
    UtilityPKG/N_UTL_StatMetricTraits.h
    UtilityPKG/N_UTL_TypeIndex.h
    UtilityPKG/N_UTL_Xyce.h
    ${CMAKE_BINARY_DIR}/src/Xyce_config.h
DESTINATION include)

# When Charon uses a specific part Xyce as a TPL, it links to code deep in
# Xyce.  Therefore, it requires a lot more headers to be installed for proper
# linking.  We're installing the headers for all Charon builds, because it is
# not worth dealing with yet another special use case.
if(${TriBITS_prefix} STREQUAL "tcad-charon")
  install(FILES 
    AnalysisPKG/N_ANP_AnalysisBase.h
    AnalysisPKG/N_ANP_AnalysisEvent.h
    AnalysisPKG/N_ANP_AnalysisManager.h
    AnalysisPKG/N_ANP_RegisterAnalysis.h
    AnalysisPKG/N_ANP_SecondLevelManager.h
    AnalysisPKG/N_ANP_StepEvent.h
    AnalysisPKG/N_ANP_SweepParam.h
    CircuitPKG/N_CIR_SecondLevelSimulator.h
    DeviceModelPKG/Core/N_DEV_CompositeParam.h
    DeviceModelPKG/Core/N_DEV_ScalingVars.h
    DeviceModelPKG/Core/N_DEV_SpecieSource.h
    DeviceModelPKG/Core/N_DEV_TransportHelper.h
    DeviceModelPKG/SandiaModels/N_DEV_1D_QASPR.h
    DeviceModelPKG/TCADModels/N_DEV_BernouliSupport.h
    DeviceModelPKG/TCADModels/N_DEV_DevicePDEInstance.h
    DeviceModelPKG/TCADModels/N_DEV_DevicePDEModel.h
    DeviceModelPKG/TCADModels/N_DEV_DiodePDE.h
    DeviceModelPKG/TCADModels/N_DEV_DopeInfo.h
    DeviceModelPKG/TCADModels/N_DEV_FermiIntegrals.h
    DeviceModelPKG/TCADModels/N_DEV_MaterialLayer.h
    DeviceModelPKG/TCADModels/N_DEV_MaterialSupport.h
    DeviceModelPKG/TCADModels/N_DEV_PDE_Electrode.h
    DeviceModelPKG/TCADModels/N_DEV_bcData.h
    ErrorHandlingPKG/N_ERH_ErrorMgr.h
    TimeIntegrationPKG/N_TIA_TIAParams.h
    UtilityPKG/ExpressionSrc/astRandEnum.h
    UtilityPKG/ExpressionSrc/expressionParamTypes.h
    UtilityPKG/N_UTL_Dump.h
    UtilityPKG/N_UTL_ExtendedString.h
    UtilityPKG/N_UTL_Factory.h
    UtilityPKG/N_UTL_Interpolators.h
    UtilityPKG/N_UTL_Listener.h
    UtilityPKG/N_UTL_Timer.h
    UtilityPKG/N_UTL_WallTime.h
  DESTINATION include)
endif()

# valid Cmake config package & version names are:
# Xyce-config.cmake
# Xyce-config-version.cmake
# or
# XyceConfig.cmake
# XyceConfigVersion.cmake
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/cmake/XyceConfigVersion.cmake"
  VERSION ${Xyce_VERSION}
  COMPATIBILITY AnyNewerVersion
)

export(EXPORT XyceLibTarget
  FILE "${CMAKE_CURRENT_BINARY_DIR}/cmake/XyceTargets.cmake"
  NAMESPACE Xyce::
)

set( Xyce_Trilinos_DIR ${${TriBITS_prefix}_DIR})
configure_package_config_file(
  cmake/XyceConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/XyceConfig.cmake
  INSTALL_DESTINATION lib/cmake
  PATH_VARS Xyce_Trilinos_DIR
)

set(ConfigPackageLocation lib/cmake/Xyce)
install(EXPORT XyceLibTarget
  FILE
    XyceTargets.cmake
  NAMESPACE
    Xyce::
  DESTINATION
    ${ConfigPackageLocation}
)

install(
  FILES
    ${CMAKE_CURRENT_BINARY_DIR}/cmake/XyceConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/cmake/XyceConfigVersion.cmake
  DESTINATION
    ${ConfigPackageLocation}
  COMPONENT
    Devel
)
