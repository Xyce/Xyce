#
# CMakeLists.txt for UtilityPKG
#
# Author: Aaron Gibson (asgibso@sandia.gov)
#
include_directories( "${CMAKE_CURRENT_BINARY_DIR}" )

#
# Handle FFT options
#
set( FFT_SOURCES )
if( Xyce_USE_FFT )
    if( Xyce_USE_INTEL_FFT )
        set( FFT_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/N_UTL_IntelFFT_Interface.C )
    elseif( Xyce_USE_FFTW )
        set( FFT_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/N_UTL_FFTW_Interface.C )
    endif()
endif()

# Set the source files for this package.
set( Utility_SOURCES
    ${FFT_SOURCES}
    ${CMAKE_CURRENT_SOURCE_DIR}/N_UTL_BreakPoint.C
    ${CMAKE_CURRENT_SOURCE_DIR}/N_UTL_Demangle.C
    ${CMAKE_CURRENT_SOURCE_DIR}/N_UTL_DeviceNameConverters.C
    ${CMAKE_CURRENT_SOURCE_DIR}/N_UTL_Diagnostic.C
    ${CMAKE_CURRENT_SOURCE_DIR}/N_UTL_Expression.C
    ${CMAKE_CURRENT_SOURCE_DIR}/N_UTL_ExpressionData.C
    ${CMAKE_CURRENT_SOURCE_DIR}/N_UTL_ExpressionInternals.C
    ${CMAKE_CURRENT_SOURCE_DIR}/N_UTL_APFT.C
    ${CMAKE_CURRENT_SOURCE_DIR}/N_UTL_FFTInterface.C
    ${CMAKE_CURRENT_SOURCE_DIR}/N_UTL_LogStream.C
    ${CMAKE_CURRENT_SOURCE_DIR}/N_UTL_Misc.C
    ${CMAKE_CURRENT_SOURCE_DIR}/N_UTL_MathSpecialFunctions.C
    ${CMAKE_CURRENT_SOURCE_DIR}/N_UTL_ExtendedString.C
    ${CMAKE_CURRENT_SOURCE_DIR}/N_UTL_NoCase.C
    ${CMAKE_CURRENT_SOURCE_DIR}/N_UTL_Op.C
    ${CMAKE_CURRENT_SOURCE_DIR}/N_UTL_OpBuilder.C
    ${CMAKE_CURRENT_SOURCE_DIR}/N_UTL_OptionBlock.C
    ${CMAKE_CURRENT_SOURCE_DIR}/N_UTL_OutputAPIHelpers.C
    ${CMAKE_CURRENT_SOURCE_DIR}/N_UTL_Param.C
    ${CMAKE_CURRENT_SOURCE_DIR}/N_UTL_Marshal.C
    ${CMAKE_CURRENT_SOURCE_DIR}/N_UTL_JSON.C
    ${CMAKE_CURRENT_SOURCE_DIR}/N_UTL_SendCURL.C
    ${CMAKE_CURRENT_SOURCE_DIR}/N_UTL_RandomNumbers.C
    ${CMAKE_CURRENT_SOURCE_DIR}/N_UTL_ReportHandler.C
    ${CMAKE_CURRENT_SOURCE_DIR}/N_UTL_Registry.C
    ${CMAKE_CURRENT_SOURCE_DIR}/N_UTL_Timer.C
    ${CMAKE_CURRENT_SOURCE_DIR}/N_UTL_CPUTime.C
    ${CMAKE_CURRENT_SOURCE_DIR}/N_UTL_MallocUsed.C
    ${CMAKE_CURRENT_SOURCE_DIR}/N_UTL_FormatMemorySize.C
    ${CMAKE_CURRENT_SOURCE_DIR}/N_UTL_FormatTime.C
    ${CMAKE_CURRENT_SOURCE_DIR}/N_UTL_StatMetricTraits.C
    ${CMAKE_CURRENT_SOURCE_DIR}/N_UTL_Stats.C
    ${CMAKE_CURRENT_SOURCE_DIR}/N_UTL_Platform.C
    ${CMAKE_CURRENT_SOURCE_DIR}/N_UTL_PrintStats.C
    ${CMAKE_CURRENT_SOURCE_DIR}/N_UTL_PrintTable.C
    ${CMAKE_CURRENT_SOURCE_DIR}/N_UTL_WallTime.C
    ${CMAKE_CURRENT_SOURCE_DIR}/N_UTL_Version.C
    ${CMAKE_CURRENT_SOURCE_DIR}/N_UTL_NameLevelKey.C
    ${CMAKE_CURRENT_SOURCE_DIR}/N_UTL_NetlistLocation.C
    ${CMAKE_CURRENT_SOURCE_DIR}/N_UTL_NodeSymbols.C
    ${CMAKE_CURRENT_SOURCE_DIR}/N_UTL_CheckIfValidFile.C
)

add_to_xyce_library( NAME Utility
    FILES ${Utility_SOURCES}
)

# Add the appropriate utility headers for the install.
set( Utility_Headers
    N_UTL_BreakPoint.h
    N_UTL_Diagnostic.h
    N_UTL_Expression.h
    N_UTL_ExpressionSymbolTable.h
    N_UTL_ExpressionData.h
    N_UTL_fwd.h
    N_UTL_FeatureTest.h
    N_UTL_FormatTime.h
    N_UTL_Interface_Enum_Types.h
    N_UTL_JSON.h
    N_UTL_MachDepParams.h
    N_UTL_Marshal.h
    N_UTL_IndexPair.h
    N_UTL_Misc.h
    N_UTL_Math.h
    N_UTL_NameLevelKey.h
    N_UTL_NetlistLocation.h
    N_UTL_NoCase.h
    N_UTL_NodeSymbols.h
    N_UTL_Op.h
    N_UTL_Pack.h
    N_UTL_Param.h
    N_UTL_ReportHandler.h
    N_UTL_Stats.h
    N_UTL_StatMetricTraits.h
    N_UTL_TypeIndex.h
    N_UTL_Xyce.h
)
install( FILES ${Utility_Headers}
    DESTINATION include
)

