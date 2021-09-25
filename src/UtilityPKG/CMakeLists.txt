#
# add sub directories
add_subdirectory(ExpressionSrc)
add_subdirectory(ExpressionTest)

# generate a timestamp.h
string(TIMESTAMP CONFIGURE_TIME %Y%m%d%H%M)
execute_process(
     COMMAND git -C ${CMAKE_SOURCE_DIR} describe --dirty
     OUTPUT_VARIABLE GITSHA  ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
if (NOT ${GITSHA} STREQUAL "")
     set(GITSTRING "-(${GITSHA})")
endif ()
set(XYCEBUILDTIMESTAMP "${CONFIGURE_TIME}${GITSTRING}")
configure_file(timestamp.h.cmake timestamp.h)


# -- process files for this library ------------------------------------------

# process source file options
if ( Xyce_USE_INTEL_FFT )
  target_sources(XyceLib PRIVATE N_UTL_IntelFFT_Interface.C )
endif ( Xyce_USE_INTEL_FFT )

if ( Xyce_USE_FFTW )
  target_sources(XyceLib PRIVATE N_UTL_FFTW_Interface.C )
endif ( Xyce_USE_FFTW )


# class source list
target_sources(XyceLib PRIVATE
      N_UTL_BreakPoint.C
      N_UTL_Demangle.C
      N_UTL_DeviceNameConverters.C
      N_UTL_Diagnostic.C
      N_UTL_Expression.C
      N_UTL_ExpressionData.C
      N_UTL_APFT.C
      N_UTL_FFTInterface.C
      N_UTL_LogStream.C
      N_UTL_Misc.C
      N_UTL_MathSpecialFunctions.C
      N_UTL_ExtendedString.C
      N_UTL_NoCase.C
      N_UTL_Op.C
      N_UTL_OpBuilder.C
      N_UTL_OptionBlock.C
      N_UTL_OutputAPIHelpers.C
      N_UTL_Param.C
      N_UTL_Marshal.C
      N_UTL_JSON.C
      N_UTL_SendCURL.C
      N_UTL_RandomNumbers.C
      N_UTL_ReportHandler.C
      N_UTL_Timer.C
      N_UTL_CPUTime.C
      N_UTL_MallocUsed.C
      N_UTL_FormatMemorySize.C
      N_UTL_FormatTime.C
      N_UTL_StatMetricTraits.C
      N_UTL_Stats.C
      N_UTL_Platform.C
      N_UTL_PrintStats.C
      N_UTL_PrintTable.C
      N_UTL_WallTime.C
      N_UTL_Version.C
      N_UTL_NameLevelKey.C
      N_UTL_NetlistLocation.C
      N_UTL_NodeSymbols.C
      N_UTL_CheckIfValidFile.C
      N_UTL_RFparams.C)

# add header files to project
target_include_directories(XyceLib PUBLIC
      ${CMAKE_CURRENT_SOURCE_DIR}
      ${CMAKE_CURRENT_BINARY_DIR} )




# -- build requested libraries -----------------------------------------------
set ( XYCE_PACKAGE_LIBRARIES ${XYCE_PACKAGE_LIBRARIES}
      Utility
      CACHE INTERNAL "X_L" )
