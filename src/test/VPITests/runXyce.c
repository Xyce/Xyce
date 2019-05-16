#include <vpi_user.h>
#include <stdio.h>
#include <N_CIR_XyceCInterface.h>

static int runXyce_compiletf(char*user_data)
{ 
      return 0;
}

static int runXyce_calltf(char*user_data)
{
      printf("In calltf for runXyce\n");

      // Used as a pointer to a pointer to an N_CIR_Xyce object.
      // This somewhat convoluted syntax is needed to stop p from
      // pointing at the same address as the VPI system task.
      void** p = (void **) malloc( sizeof(void* [1]) );  

      // Turn the desired Xyce command line invocation into an int and
      // char** pointer that can used to initialize an N_CIR_ Xyce object.  
      // This is hard-coded for initial testing purposes.
      //char *argList[5] = {"Xyce","-quiet","-o", "testOutput","test.cir"};
      char *argList[] = {
	  (char*)("Xyce"),
          (char*)("runXyce.cir")
      };
      int argc = sizeof(argList)/sizeof(argList[0]);
      char** argv = argList;

      // Test the methods in utils/XyceCInterface/N_CIR_XyceCInterface.C
      xyce_open(p);
      xyce_initialize(p,argc,argv);
      xyce_runSimulation(p);
      xyce_close(p);

      // pointer clean-up
      free( p );

      vpi_printf("Exiting calltf for runXyce\n"); 
      return 0;
}

void runXyce_register()
{
      s_vpi_systf_data tf_data;

      tf_data.type      = vpiSysTask;
      tf_data.tfname    = "$runXyce";
      tf_data.calltf    = runXyce_calltf;
      tf_data.compiletf = runXyce_compiletf;
      tf_data.sizetf    = 0;
      tf_data.user_data = 0;
      vpi_register_systf(&tf_data);
}

void (*vlog_startup_routines[])() = {
    runXyce_register,
    0 /* final entry must be zero */
};

