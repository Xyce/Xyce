
AM_CPPFLAGS = @Xyce_INCS@


# standalone DeviceInterfaceTest executable
check_PROGRAMS = DeviceInterfaceTest

DeviceInterfaceTESTSOURCES = \
  DeviceInterfaceTest.C 

# standalone DeviceInterface  executable
DeviceInterfaceTest_SOURCES = $(DeviceInterfaceTESTSOURCES)
DeviceInterfaceTest_LDADD = $(top_builddir)/src/libxyce.la
DeviceInterfaceTest_LDFLAGS = -static $(AM_LDFLAGS)
 
