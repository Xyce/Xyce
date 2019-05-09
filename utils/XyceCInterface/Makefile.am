
AM_CPPFLAGS = @Xyce_INCS@

lib_LTLIBRARIES = libxycecinterface.la
libxycecinterface_la_SOURCES = N_CIR_XyceCInterface.C \
    N_CIR_XyceCInterface.h 


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

libxycecinterface_la_LIBADD = $(XYCELINK) $(top_builddir)/src/libxyce.la $(NONFREELD) $(RADLD) $(MODSPECLD) $(ADMSLD) $(NEURONLD)
libxycecinterface_la_LDFLAGS = $(AM_LDFLAGS) @XYCELIBS@

include_HEADERS = \
  $(srcdir)/N_CIR_XyceCInterface.h 


