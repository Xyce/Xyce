#
# add sub directories
add_subdirectory(ExpressionSrc)

# generate a timestamp.h
string(TIMESTAMP CONFIGURE_TIME %Y%m%d%H%M)
execute_process(
     COMMAND git -C ${CMAKE_SOURCE_DIR} log --pretty=format:%h -1
     OUTPUT_VARIABLE GITSHA  ERROR_QUIET)
if (NOT ${GITSHA} STREQUAL "")
     set(GITSTRING "-g-${GITSHA}")
endif ()
set(XYCEBUILDTIMESTAMP "${CONFIGURE_TIME}${GITSTRING}")
configure_file(timestamp.h.cmake timestamp.h)


# -- process files for this library ------------------------------------------

# process source file options
if ( Xyce_USE_INTEL_FFT )
  set ( FFT_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/N_UTL_IntelFFT_Interface.C )
endif ( Xyce_USE_INTEL_FFT )

if ( Xyce_USE_FFTW )
  set ( FFTW_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/N_UTL_FFTW_Interface.C )
endif ( Xyce_USE_FFTW )


# class source list
set ( Utility_SOURCES
      ${FFT_SOURCES}
      ${FFTW_SOURCES}
      ${CMAKE_CURRENT_SOURCE_DIR}/N_UTL_BreakPoint.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_UTL_Demangle.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_UTL_DeviceNameConverters.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_UTL_Diagnostic.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_UTL_Expression.C
      ${CMAKE_CURRENT_SOURCE_DIR}/N_UTL_ExpressionData.C
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
      ${CMAKE_CURRENT_SOURCE_DIR}/N_UTL_RFparams.C
      CACHE INTERNAL "X_S" )

set ( XYCE_SOURCES ${XYCE_SOURCES} ${Utility_SOURCES} PARENT_SCOPE)
# add header files to project
set ( XYCE_INCLUDE_DIRECTORIES ${XYCE_INCLUDE_DIRECTORIES}
      ${CMAKE_CURRENT_SOURCE_DIR}
      ${CMAKE_CURRENT_BINARY_DIR} 
      PARENT_SCOPE )




# -- build requested libraries -----------------------------------------------
set ( XYCE_PACKAGE_LIBRARIES ${XYCE_PACKAGE_LIBRARIES}
      Utility
      CACHE INTERNAL "X_L" )
