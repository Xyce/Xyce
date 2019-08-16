
AM_CPPFLAGS = @Xyce_INCS@

SUBDIRS = \
  AnalysisPKG \
  CircuitPKG \
  DakotaLinkPKG \
  DeviceModelPKG \
  ErrorHandlingPKG \
  IOInterfacePKG \
  LinearAlgebraServicesPKG \
  LoaderServicesPKG \
  NonlinearSolverPKG \
  ParallelDistPKG \
  TimeIntegrationPKG \
  TopoManagerPKG \
  UtilityPKG \
  MultiTimePDEPKG \
  test

XYCELINK = \
  $(top_builddir)/src/CircuitPKG/libCircuit.la \
  $(top_builddir)/src/DakotaLinkPKG/libDakotaLink.la \
  $(top_builddir)/src/DeviceModelPKG/Core/libDeviceCore.la \
  $(top_builddir)/src/DeviceModelPKG/OpenModels/libOpenModels.la \
  $(top_builddir)/src/DeviceModelPKG/TCADModels/libTCADModels.la \
  $(top_builddir)/src/DeviceModelPKG/IBISModels/libIBISModels.la \
  $(top_builddir)/src/DeviceModelPKG/EXTSC/libEXT.la \
  $(top_builddir)/src/ErrorHandlingPKG/libErrorHandling.la \
  $(top_builddir)/src/IOInterfacePKG/libIOInterface.la \
  $(top_builddir)/src/IOInterfacePKG/Output/libOutput.la \
  $(top_builddir)/src/LinearAlgebraServicesPKG/libLinearAlgebraServices.la \
  $(top_builddir)/src/LoaderServicesPKG/libLoaderServices.la \
  $(top_builddir)/src/NonlinearSolverPKG/libNonlinearSolver.la \
  $(top_builddir)/src/ParallelDistPKG/libParallelDist.la \
  $(top_builddir)/src/AnalysisPKG/libAnalysis.la \
  $(top_builddir)/src/TimeIntegrationPKG/libTimeIntegration.la \
  $(top_builddir)/src/TopoManagerPKG/libTopoManager.la \
  $(top_builddir)/src/UtilityPKG/libUtility.la \
  $(top_builddir)/src/MultiTimePDEPKG/libMPDE.la


# conditionally link in radiation-aware device models.
if RADMODELS
    RADLD = $(top_builddir)/src/DeviceModelPKG/SandiaModels/libSandiaModels.la
else
    RADLD =
endif

# conditionally link in ADMS device models.
if ADMSMODELS
    ADMSLD = $(top_builddir)/src/DeviceModelPKG/ADMS/libADMS.la 
else
    ADMSLD =
endif

# conditionally link in Neuron device models.
if NEURONMODELS
    NEURONLD = $(top_builddir)/src/DeviceModelPKG/NeuronModels/libNeuronModels.la 
else
    NEURONLD =
endif

# conditionally link in radiation-aware device models.
if NONFREEMODELS
    NONFREELD = $(top_builddir)/src/DeviceModelPKG/Xyce_NonFree/libNonFree.la
else
    NONFREELD =
endif

# conditionally link in Berkeley device models.
if MODSPECMODELS
    MODSPECLD =   $(top_builddir)/src/DeviceModelPKG/BERKELEY/libBERKELEY.la \
                  $(top_builddir)/src/DeviceModelPKG/BERKELEY/Xyce_ModSpec_Interface/libModSpec.la 
else
    MODSPECLD =
endif



# build library (static & shared) for linking (w/libtool) into other codes
lib_LTLIBRARIES = libxyce.la
libxyce_la_SOURCES = $(XYCESOURCES)
libxyce_la_LIBADD = $(XYCELINK)
libxyce_la_LDFLAGS = $(AM_LDFLAGS) @XYCELIBS@

if BUILD_XYCE_BINARY

if BUILD_XYCE_SHAREABLE

# standalone Xyce executable (shareable object)
bin_PROGRAMS = Xyce
Xyce_SOURCES = Xyce.C
Xyce_LDADD = libxyce.la $(NONFREELD) $(RADLD) $(MODSPECLD) $(ADMSLD) $(NEURONLD)
Xyce_LDFLAGS = $(AM_LDFLAGS)

else

# standalone Xyce executable (monolithic)
bin_PROGRAMS = Xyce
Xyce_SOURCES = Xyce.C
Xyce_LDADD = libxyce.la $(NONFREELD) $(RADLD) $(MODSPECLD) $(ADMSLD) $(NEURONLD)
Xyce_LDFLAGS = -static $(AM_LDFLAGS)

endif

endif

# install headers needed to compile other codes with libxyce.*
include_HEADERS = \
  $(srcdir)/AnalysisPKG/N_ANP_fwd.h \
  $(srcdir)/AnalysisPKG/N_ANP_NoiseData.h \
  $(srcdir)/CircuitPKG/N_CIR_Xyce.h \
  $(srcdir)/CircuitPKG/N_CIR_Xygra.h \
  $(srcdir)/CircuitPKG/N_CIR_GenCouplingSimulator.h \
  $(srcdir)/DeviceModelPKG/Core/N_DEV_Configuration.h \
  $(srcdir)/DeviceModelPKG/Core/N_DEV_Const.h \
  $(srcdir)/DeviceModelPKG/Core/N_DEV_Device.h \
  $(srcdir)/DeviceModelPKG/Core/N_DEV_DeviceBlock.h \
  $(srcdir)/DeviceModelPKG/Core/N_DEV_DeviceEntity.h \
  $(srcdir)/DeviceModelPKG/Core/N_DEV_DeviceInstance.h \
  $(srcdir)/DeviceModelPKG/Core/N_DEV_DeviceMaster.h \
  $(srcdir)/DeviceModelPKG/Core/N_DEV_DeviceModel.h \
  $(srcdir)/DeviceModelPKG/Core/N_DEV_DeviceOptions.h \
  $(srcdir)/DeviceModelPKG/Core/N_DEV_DeviceSupport.h \
  $(srcdir)/DeviceModelPKG/Core/N_DEV_InstanceName.h \
  $(srcdir)/DeviceModelPKG/Core/N_DEV_ExternalSimulationData.h \
  $(srcdir)/DeviceModelPKG/Core/N_DEV_ExternData.h \
  $(srcdir)/DeviceModelPKG/Core/N_DEV_fwd.h \
  $(srcdir)/DeviceModelPKG/Core/N_DEV_Message.h \
  $(srcdir)/DeviceModelPKG/Core/N_DEV_Param.h \
  $(srcdir)/DeviceModelPKG/Core/N_DEV_Pars.h \
  $(srcdir)/DeviceModelPKG/Core/N_DEV_SolverState.h \
  $(srcdir)/DeviceModelPKG/Core/N_DEV_Units.h \
  $(srcdir)/DeviceModelPKG/Core/N_DEV_VectorComputeInterface.h \
  $(srcdir)/DeviceModelPKG/OpenModels/N_DEV_ADC.h \
  $(srcdir)/DeviceModelPKG/OpenModels/N_DEV_DAC.h \
  $(srcdir)/DeviceModelPKG/OpenModels/N_DEV_Resistor.h \
  $(srcdir)/DeviceModelPKG/OpenModels/N_DEV_Capacitor.h \
  $(srcdir)/DeviceModelPKG/OpenModels/N_DEV_BJT.h \
  $(srcdir)/DeviceModelPKG/OpenModels/N_DEV_Diode.h \
  $(srcdir)/DeviceModelPKG/OpenModels/N_DEV_MOSFET1.h \
  $(srcdir)/ErrorHandlingPKG/N_ERH_fwd.h \
  $(srcdir)/ErrorHandlingPKG/N_ERH_Message.h \
  $(srcdir)/IOInterfacePKG/N_IO_CmdParse.h \
  $(srcdir)/IOInterfacePKG/N_IO_ExtOutInterface.h \
  $(srcdir)/IOInterfacePKG/N_IO_HangingResistor.h \
  $(srcdir)/IOInterfacePKG/N_IO_OutputTypes.h \
  $(srcdir)/IOInterfacePKG/N_IO_fwd.h \
  $(srcdir)/LinearAlgebraServicesPKG/N_LAS_fwd.h \
  $(srcdir)/LinearAlgebraServicesPKG/N_LAS_Vector.h \
  $(srcdir)/LinearAlgebraServicesPKG/N_LAS_Matrix.h \
  $(srcdir)/LinearAlgebraServicesPKG/N_LAS_MultiVector.h \
  $(srcdir)/LoaderServicesPKG/N_LOA_fwd.h \
  $(srcdir)/LoaderServicesPKG/N_LOA_Loader.h \
  $(srcdir)/NonlinearSolverPKG/N_NLS_fwd.h \
  $(srcdir)/NonlinearSolverPKG/N_NLS_TwoLevelEnum.h \
  $(srcdir)/ParallelDistPKG/N_PDS_Manager.h \
  $(srcdir)/ParallelDistPKG/N_PDS_fwd.h \
  $(srcdir)/ParallelDistPKG/N_PDS_MPI.h \
  $(srcdir)/ParallelDistPKG/N_PDS_ParallelMachine.h \
  $(srcdir)/ParallelDistPKG/N_PDS_Serial.h \
  $(srcdir)/TimeIntegrationPKG/N_TIA_TwoLevelError.h \
  $(srcdir)/TimeIntegrationPKG/N_TIA_fwd.h \
  $(srcdir)/TopoManagerPKG/N_TOP_fwd.h \
  $(srcdir)/UtilityPKG/N_UTL_BreakPoint.h \
  $(srcdir)/UtilityPKG/N_UTL_Diagnostic.h \
  $(srcdir)/UtilityPKG/N_UTL_Expression.h \
  $(srcdir)/UtilityPKG/N_UTL_ExpressionSymbolTable.h \
  $(srcdir)/UtilityPKG/N_UTL_ExpressionData.h \
  $(srcdir)/UtilityPKG/N_UTL_fwd.h \
  $(srcdir)/UtilityPKG/N_UTL_FeatureTest.h \
  $(srcdir)/UtilityPKG/N_UTL_FormatTime.h \
  $(srcdir)/UtilityPKG/N_UTL_Interface_Enum_Types.h \
  $(srcdir)/UtilityPKG/N_UTL_JSON.h \
  $(srcdir)/UtilityPKG/N_UTL_MachDepParams.h \
  $(srcdir)/UtilityPKG/N_UTL_Marshal.h \
  $(srcdir)/UtilityPKG/N_UTL_Misc.h \
  $(srcdir)/UtilityPKG/N_UTL_Math.h \
  $(srcdir)/UtilityPKG/N_UTL_NameLevelKey.h \
  $(srcdir)/UtilityPKG/N_UTL_NetlistLocation.h \
  $(srcdir)/UtilityPKG/N_UTL_NoCase.h \
  $(srcdir)/UtilityPKG/N_UTL_NodeSymbols.h \
  $(srcdir)/UtilityPKG/N_UTL_Op.h \
  $(srcdir)/UtilityPKG/N_UTL_Pack.h \
  $(srcdir)/UtilityPKG/N_UTL_Param.h \
  $(srcdir)/UtilityPKG/N_UTL_ReportHandler.h \
  $(srcdir)/UtilityPKG/N_UTL_Stats.h \
  $(srcdir)/UtilityPKG/N_UTL_StatMetricTraits.h \
  $(srcdir)/UtilityPKG/N_UTL_TypeIndex.h \
  $(srcdir)/UtilityPKG/N_UTL_Xyce.h \
  $(top_builddir)/src/Xyce_config.h

xmldir = $(datadir)/xml
xml_DATA= $(top_srcdir)/utils/ADMS/adms.implicit.xml \
	$(top_srcdir)/utils/ADMS/xyceVersion.xml \
	$(top_srcdir)/utils/ADMS/xyceBasicTemplates.xml \
	$(top_srcdir)/utils/ADMS/xyceImplementationFile.xml \
	$(top_srcdir)/utils/ADMS/xyceHeaderFile.xml \
	$(top_srcdir)/utils/ADMS/xyceVersion_nosac.xml \
	$(top_srcdir)/utils/ADMS/xyceBasicTemplates_nosac.xml \
	$(top_srcdir)/utils/ADMS/xyceAnalogFunction_nosac.xml \
	$(top_srcdir)/utils/ADMS/xyceImplementationFile_nosac.xml \
	$(top_srcdir)/utils/ADMS/xyceHeaderFile_nosac.xml \
	$(top_srcdir)/utils/ADMS/xyceMakefile.am.xml \
	$(top_srcdir)/utils/ADMS/xyceBootstrapFile.xml \
	$(top_srcdir)/utils/ADMS/html_params.xml 
