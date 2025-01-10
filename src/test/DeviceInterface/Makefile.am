
AM_CPPFLAGS = @Xyce_INCS@


# standalone DeviceInterfaceTest executable
check_PROGRAMS = DeviceInterfaceTest ADC_DACRunTest

DeviceInterfaceTESTSOURCES = \
  DeviceInterfaceTest.C 

# standalone DeviceInterface  executable
DeviceInterfaceTest_SOURCES = $(DeviceInterfaceTESTSOURCES)
DeviceInterfaceTest_LDADD = $(top_builddir)/src/libxyce.la
DeviceInterfaceTest_LDFLAGS = -static $(AM_LDFLAGS)

ADC_DACRunTestTESTSOURCES = \
  ADC_DACRunTest.C 

# standalone ADC_DACRunTest executable
ADC_DACRunTest_SOURCES = $(ADC_DACRunTestTESTSOURCES)
ADC_DACRunTest_LDADD = $(top_builddir)/src/libxyce.la
ADC_DACRunTest_LDFLAGS = -static $(AM_LDFLAGS)
 
 
