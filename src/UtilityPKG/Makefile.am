
AM_CPPFLAGS = @Xyce_INCS@ -DXYCEBUILDTIMESTAMP=\"`$(top_srcdir)/utils/XyceDatestamp.sh $(top_srcdir)`\"

SUBDIRS = ExpressionSrc

noinst_LTLIBRARIES = libUtility.la

if BUILD_INTEL_FFT_SUPPORT
intelfft_SOURCES = \
  N_UTL_IntelFFT_Interface.C \
  N_UTL_IntelFFT_Interface.hpp
else
intelfft_SOURCES =
endif

if BUILD_FFTW_SUPPORT
fftw_SOURCES = \
  N_UTL_FFTW_Interface.C \
  N_UTL_FFTW_Interface.hpp
else
fftw_SOURCES =
endif

if THREADS 
threads_LIB = -lpthread 
threads_CFLAGS = -pthread 
else 
threads_LIB = 
threads_CFLAGS = 
endif 

libUtility_la_SOURCES = \
  $(intelfft_SOURCES) \
  $(fftw_SOURCES) \
  N_UTL_BreakPoint.C \
  N_UTL_Demangle.C \
  N_UTL_DeviceNameConverters.C \
  N_UTL_Diagnostic.C \
  N_UTL_NoCase.C \
  N_UTL_Expression.C \
  N_UTL_ExpressionData.C \
  N_UTL_APFT.C \
  N_UTL_FFTInterface.C \
  N_UTL_LogStream.C \
  N_UTL_ReportHandler.C \
  N_UTL_Marshal.C \
  N_UTL_JSON.C \
  N_UTL_SendCURL.C \
  N_UTL_ByThreadStreambuf.C \
  N_UTL_Misc.C \
  N_UTL_ExtendedString.C \
  N_UTL_NetlistLocation.C \
  N_UTL_NodeSymbols.C \
  N_UTL_OptionBlock.C \
  N_UTL_Op.C \
  N_UTL_OpBuilder.C \
  N_UTL_OutputAPIHelpers.C \
  N_UTL_Param.C \
  N_UTL_RandomNumbers.C \
  N_UTL_RFparams.C \
  N_UTL_Timer.C \
  N_UTL_CPUTime.C \
  N_UTL_MallocUsed.C \
  N_UTL_FormatMemorySize.C \
  N_UTL_FormatTime.C \
  N_UTL_StatMetricTraits.C \
  N_UTL_Stats.C \
  N_UTL_Platform.C \
  N_UTL_PrintTable.C \
  N_UTL_PrintStats.C \
  N_UTL_WallTime.C \
  N_UTL_Version.C \
  N_UTL_Algorithm.C \
  N_UTL_NameLevelKey.C \
  N_UTL_CheckIfValidFile.C \
  N_UTL_MathSpecialFunctions.C \
  N_UTL_AssemblyTypes.h \
  N_UTL_TypeIndex.h \
  N_UTL_NameLevelKey.h \
  N_UTL_DeleteList.h \
  N_UTL_Interface_Enum_Types.h \
  N_UTL_Interpolators.h \
  N_UTL_ExpressionSymbolTable.h \
  N_UTL_Expression.h \
  N_UTL_ExpressionData.h \
  N_UTL_FFTInterface.hpp \
  N_UTL_RandomNumbers.h \
  N_UTL_RFparams.h \
  N_UTL_Demangle.h \
  N_UTL_DeviceNameConverters.h \
  N_UTL_Diagnostic.h \
  N_UTL_IndentStreamBuf.h \
  N_UTL_TeeStreamBuf.h \
  N_UTL_PThread.h \
  N_UTL_ByThreadStreambuf.h \
  N_UTL_LogStream.h \
  N_UTL_ReportHandler.h \
  N_UTL_Marshal.h \
  N_UTL_JSON.h \
  N_UTL_SendCURL.h \
  N_UTL_Op.h \
  N_UTL_OpBuilder.h \
  N_UTL_OutputAPIHelpers.h \
  N_UTL_APFT.h \
  N_UTL_DFTInterfaceDecl.hpp \
  N_UTL_FFTInterfaceDecl.hpp \
  N_UTL_FixedQueue.h \
  N_UTL_Functors.h \
  N_UTL_MachDepParams.h \
  N_UTL_Listener.h \
  N_UTL_Dump.h \
  N_UTL_Misc.h \
  N_UTL_ExtendedString.h \
  N_UTL_NetlistLocation.h \
  N_UTL_NodeSymbols.h \
  N_UTL_SaveIOSState.h \
  N_UTL_BreakPoint.h \
  N_UTL_Factory.h \
  N_UTL_OptionBlock.h \
  N_UTL_Pack.h \
  N_UTL_Param.h \
  N_UTL_Timer.h \
  N_UTL_CPUTime.h \
  N_UTL_FormatMemorySize.h \
  N_UTL_FormatTime.h \
  N_UTL_MallocUsed.h \
  N_UTL_StatMetricTraits.h \
  N_UTL_Stats.h \
  N_UTL_Platform.h \
  N_UTL_PrintStats.h \
  N_UTL_PrintTable.h \
  N_UTL_WallTime.h \
  N_UTL_Graph.h \
  N_UTL_Xyce.h \
  N_UTL_FeatureTest.h \
  N_UTL_IndentStreamBuf.h \
  N_UTL_NoCase.h \
  N_UTL_Version.h \
  N_UTL_Algorithm.h \
  N_UTL_CheckIfValidFile.h \
  N_UTL_Math.h \
  N_UTL_MathSpecialFunctions.h  \
	N_UTL_HspiceBools.h 

BUILT_SOURCES = timestamp-generator

CLEANFILES = timestamp.h

timestamp.h: timestamp-generator

timestamp-generator:
	echo "/* TIMESTAMP APPLIED */" > timestamp.h

