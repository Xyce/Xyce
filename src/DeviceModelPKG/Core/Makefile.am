AM_CPPFLAGS = @Xyce_INCS@

noinst_LTLIBRARIES = libDeviceCore.la

# these are files built by bison and flex
if BUILD_REACTION_PARSER
BUILT_SOURCES= N_DEV_ReactionParser.cxx \
            N_DEV_ReactionParser.hxx \
            N_DEV_ReactionLexer.cxx \
            location.hh \
            position.hh \
	    stack.hh
else
BUILT_SOURCES=
endif

CLEANFILES= $(BUILT_SOURCES)

PARSER_HEADER_DEPS =
if BUILD_REACTION_PARSER
PARSER_HEADER_DEPS += N_DEV_ReactionParser.hxx
endif

#these are extra parser-related files that we want to include conditionally
PARSER_SRCS= N_DEV_ReactionLexer.h

libDeviceCore_la_SOURCES = \
  N_DEV_RegisterDevices.C \
  N_DEV_Algorithm.C \
  N_DEV_Print.C \
  N_DEV_CompositeParam.C \
  N_DEV_Configuration.C \
  N_DEV_DeviceBlock.C \
  N_DEV_DeviceEntity.C \
  N_DEV_DeviceInstance.C \
  N_DEV_DeviceMgr.C \
  N_DEV_DeviceModel.C \
  N_DEV_DeviceOptions.C \
  N_DEV_DeviceSensitivities.C \
  N_DEV_DeviceState.C \
  N_DEV_DeviceSupport.C \
  N_DEV_DeviceMaster.C \
  N_DEV_Dump.C \
  N_DEV_InstanceName.C \
  N_DEV_MatrixLoadData.C \
  N_DEV_Message.C \
  N_DEV_NumericalJacobian.C \
  N_DEV_LaTexDoc.C \
  N_DEV_Op.C \
  N_DEV_OpBuilders.C \
  N_DEV_Param.C \
  N_DEV_Pars.C \
  N_DEV_RateConstantCalculators.C \
  N_DEV_Reaction.C \
  N_DEV_ReactionNetwork.C \
  N_DEV_Region.C \
  N_DEV_RegionData.C \
  N_DEV_RxnSet.C \
  N_DEV_ScalingVars.C \
  N_DEV_SolverState.C \
  N_DEV_Source.C \
  N_DEV_SourceData.C \
  N_DEV_Specie.C \
  N_DEV_SpecieSource.C \
  N_DEV_Units.C \
  N_DEV_fwd.h \
  N_DEV_RegisterDevices.h \
  N_DEV_ExpressionGroupWrapper.h \
  N_DEV_Algorithm.h \
  N_DEV_Print.h \
  N_DEV_CompositeParam.h \
  N_DEV_Configuration.h \
  N_DEV_Const.h \
  N_DEV_Device.h \
  N_DEV_DeviceBlock.h \
  N_DEV_DeviceEntity.h \
  N_DEV_DeviceInstance.h \
  N_DEV_DeviceMgr.h \
  N_DEV_DeviceModel.h \
  N_DEV_DeviceOptions.h \
  N_DEV_DeviceSensitivities.h \
  N_DEV_DeviceState.h \
  N_DEV_DeviceSupport.h \
  N_DEV_DeviceMaster.h \
  N_DEV_Dump.h \
  N_DEV_InstanceName.h \
  N_DEV_ExternData.h \
  N_DEV_ExternalSimulationData.h \
  N_DEV_MatrixLoadData.h \
  N_DEV_Message.h \
  N_DEV_NumericalJacobian.h \
  N_DEV_LaTexDoc.h \
  N_DEV_Op.h \
  N_DEV_OpBuildersg.h \
  N_DEV_Param.h \
  N_DEV_Pars.h \
  N_DEV_RateConstantCalculators.h \
  N_DEV_Reaction.h \
  N_DEV_ReactionNetwork.h \
  N_DEV_Region.h \
  N_DEV_RegionData.h \
  N_DEV_RxnSet.h \
  N_DEV_ScalingVars.h \
  N_DEV_SolverState.h \
  N_DEV_Source.h \
  N_DEV_SourceData.h \
  N_DEV_Specie.h \
  N_DEV_SpecieSource.h \
  N_DEV_TransportHelper.h \
  N_DEV_Units.h \
  N_DEV_VectorComputeInterface.h

# these can't be in "SOURCES" above or automake screams about YACC.
EXTRA_DIST=N_DEV_ReactionLexer.l \
  N_DEV_ReactionParser.yxx

if BUILD_REACTION_PARSER
libDeviceCore_la_SOURCES += $(BUILT_SOURCES) $(PARSER_SRCS)
endif

#We have to do it this way, because automake appears to assume it can use
# yacc and bison interchangably, and we can't.

N_DEV_ReactionParser.cxx: $(srcdir)/N_DEV_ReactionParser.yxx
	$(BISON) -o$@ $(srcdir)/N_DEV_ReactionParser.yxx

N_DEV_ReactionLexer.cxx: $(srcdir)/N_DEV_ReactionLexer.l
	$(LEX) -+ -o$@ $(srcdir)/N_DEV_ReactionLexer.l

N_DEV_ReactionParser.hxx: N_DEV_ReactionParser.cxx
location.hh: N_DEV_ReactionParser.cxx
position.hh: N_DEV_ReactionParser.cxx
stack.hh: N_DEV_ReactionParser.cxx

