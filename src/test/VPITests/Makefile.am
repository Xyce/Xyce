
AM_CPPFLAGS = @Xyce_INCS@ -I$(VPI_BASE)/include/iverilog
AM_CFLAGS = @Xyce_INCS@ -I$(VPI_BASE)/include/iverilog


MKL_EXTRA = -lmkl_rt -lmkl_avx
# since the build targets all require a special build procedure (a compile followed by two link stages)
# use the standard tool to make the intermediate object xxxxInt then do the iverilog link step 
# a cusotom build rule will be needed for each these targets
# standalone executable
#check_PROGRAMS = runXyceInt runXyceInStepsInt runXyceWithDACInt

runXyce : runXyce.c
	$(CC) $(AM_CFLAGS) -c -fpic -o runXyce.o $(srcdir)/runXyce.c
	$(CC) -shared -L$(VPI_BASE)/lib -lvpi -L$(prefix)/lib -lxycecinterface \
		$(XYCELIBS) $(MKL_EXTRA) -Wl,-rpath=$(prefix)/lib -Wl,-rpath=$(top_builddir)/utils/XyceCInterface/.libs  -o runXyce.vpi runXyce.o
	iverilog -o runXyce.vvp $(srcdir)/runXyce.v

runXyceInSteps : runXyceInSteps.c
	$(CC) $(AM_CFLAGS) -c -fpic -o runXyceInSteps.o $(srcdir)/runXyceInSteps.c
	$(CC) -shared -L$(VPI_BASE)/lib -lvpi -L$(prefix)/lib -lxycecinterface \
		$(XYCELIBS) $(MKL_EXTRA) -Wl,-rpath=$(prefix)/lib -Wl,-rpath=$(top_builddir)/utils/XyceCInterface/.libs  -o runXyceInSteps.vpi runXyceInSteps.o
	iverilog -o runXyceInSteps.vvp $(srcdir)/runXyceInSteps.v


runXyceWithDAC : runXyceWithDAC.c
	$(CC) $(AM_CFLAGS) -c -fpic -o runXyceWithDAC.o $(srcdir)/runXyceWithDAC.c
	$(CC) -shared -L$(VPI_BASE)/lib -lvpi -L$(prefix)/lib -lxycecinterface \
		$(XYCELIBS) $(MKL_EXTRA) -Wl,-rpath=$(prefix)/lib -Wl,-rpath=$(top_builddir)/utils/XyceCInterface/.libs  -o runXyceWithDAC.vpi runXyceWithDAC.o
	iverilog -o runXyceWithDAC.vvp $(srcdir)/runXyceWithDAC.v



 
