
AM_CPPFLAGS = @Xyce_INCS@

# needed for Dakota 4.x not 5.0
if DAKOTA_OBJ_NEEDED 
  DAKOTA_OBJS = 
endif

# conditionally link in radiation-aware device models.
if RADMODELS
    RADLD = $(top_builddir)/src/DeviceModelPKG/SandiaModels/libSandiaModels.la
else
    RADLD = 
endif

# conditionally link in radiation-aware device models.
if NONFREEMODELS
    NONFREELD = $(top_builddir)/src/DeviceModelPKG/Xyce_NonFree/libNonFree.la
else
    NONFREELD = 
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

# standalone XyceLibTest executable
check_PROGRAMS = twoLevelNewtonLinearDC twoLevelNewtonLinearTran
#twoLevelNewton

#twoLevelNewton_SOURCES = twoLevelNewton.C
#twoLevelNewton_LDADD = $(top_builddir)/src/libxyce.la $(RADLD) $(NONFREELD) $(ADMSLD) $(NEURONLD)
#twoLevelNewton_LDFLAGS = -static $(AM_LDFLAGS) $(DAKOTA_OBJS)

twoLevelNewtonLinearDC_SOURCES = twoLevelNewtonLinearDC.C baseNewton.h twoLevelNewtonLinear.h
twoLevelNewtonLinearDC_LDADD = $(top_builddir)/src/libxyce.la $(RADLD) $(NONFREELD) $(ADMSLD) $(NEURONLD)
twoLevelNewtonLinearDC_LDFLAGS = -static $(AM_LDFLAGS) $(DAKOTA_OBJS)

twoLevelNewtonLinearTran_SOURCES = twoLevelNewtonLinearTran.C baseNewton.h twoLevelNewtonLinear.h
twoLevelNewtonLinearTran_LDADD = $(top_builddir)/src/libxyce.la $(RADLD) $(NONFREELD) $(ADMSLD) $(NEURONLD)
twoLevelNewtonLinearTran_LDFLAGS = -static $(AM_LDFLAGS) $(DAKOTA_OBJS)

