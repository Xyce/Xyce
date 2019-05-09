
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

EXTRA_DIST = XygraFile.C

# standalone XyceLibTest executable
check_PROGRAMS = XygraSource Xygra XygraFile
XygraSource_SOURCES = XygraSource.C
XygraSource_LDADD = $(top_builddir)/src/libxyce.la $(RADLD) $(NONFREELD) $(ADMSLD) $(NEURONLD)
XygraSource_LDFLAGS = -static $(AM_LDFLAGS) $(DAKOTA_OBJS)

XygraFile_SOURCES = XygraFile.C
XygraFile_LDADD = $(top_builddir)/src/libxyce.la $(RADLD) $(NONFREELD)
XygraFile_LDFLAGS = -static $(AM_LDFLAGS) $(DAKOTA_OBJS)

Xygra_SOURCES = Xygra.C
Xygra_LDADD = $(top_builddir)/src/libxyce.la $(RADLD) $(NONFREELD)
Xygra_LDFLAGS = -static $(AM_LDFLAGS) $(DAKOTA_OBJS)



