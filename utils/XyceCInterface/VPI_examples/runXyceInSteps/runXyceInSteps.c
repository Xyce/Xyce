#include <vpi_user.h>
#include <stdio.h>
#include <N_CIR_XyceCInterface.h>

static int runXyceInSteps_compiletf(char*user_data)
{ 
      return 0;
}

static int runXyceInSteps_calltf(char*user_data)
{
      printf("In calltf for runXyceInSteps\n");

      // Used as a pointer to a pointer to an N_CIR_Xyce object.
      // This somewhat convoluted syntax is needed to stop p from
      // pointing at the same address as the VPI system task.
      void** p = new void* [1];

      // Turn the desired Xyce command line invocation into an int and
      // char** pointer that can used to initialize an N_CIR_ Xyce object.  
      // This is hard-coded for initial testing purposes.
      char *argList[] = {
	  (char*)("Xyce"),
          (char*)("runXyceInSteps.cir")
      };
      int argc = sizeof(argList)/sizeof(argList[0]);
      char** argv = argList;

      // Open and initialize the N_CIR_Xyce object
      xyce_open(p);
      xyce_initialize(p,argc,argv);

      // run Xyce in steps
      int i, status;
      double requested_time, actual_time;
      double* actual_time_ptr = &actual_time;
      for (i=0; i<=10; i++)
      {
        requested_time = 0.1*(i+1);
	printf( "Calling simulateUntil with requested_time of %f\n",requested_time);
        //actual_time = 0.0
        status = xyce_simulateUntil(p, requested_time, actual_time_ptr );
	printf( "Return status from simulateUntil = %d and actual_time = %f\n",status, actual_time);
      }

      // Finalize and close the N_CIR_Xyce object
      xyce_close(p);

      // pointer clean-up
      delete[] p;

      vpi_printf("Exiting calltf for runXyceInSteps\n"); 
      return 0;
}

void runXyceInSteps_register()
{
      s_vpi_systf_data tf_data;

      tf_data.type      = vpiSysTask;
      tf_data.tfname    = "$runXyceInSteps";
      tf_data.calltf    = runXyceInSteps_calltf;
      tf_data.compiletf = runXyceInSteps_compiletf;
      tf_data.sizetf    = 0;
      tf_data.user_data = 0;
      vpi_register_systf(&tf_data);
}

void (*vlog_startup_routines[])() = {
    runXyceInSteps_register,
    0 /* final entry must be zero */
};

