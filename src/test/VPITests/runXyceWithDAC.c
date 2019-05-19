#include <vpi_user.h>
#include <stdio.h>
#include <N_CIR_XyceCInterface.h>

static int runXyceWithDAC_compiletf(char*user_data)
{ 
      return 0;
}

static int runXyceWithDAC_calltf(char*user_data)
{
      //printf() displays the string inside quotation
      printf("In calltf for runXyceWithDAC\n");

      // Used as a pointer to a pointer to an N_CIR_Xyce object.
      // This somewhat convoluted syntax is needed to stop p from
      // pointing at the same address as the VPI system task.
      void** p = (void **) malloc( sizeof(void* [1]) );  

      // Turn the desired Xyce command line invocation into an int and
      // char** pointer that can used to initialize an N_CIR_ Xyce object.  
      // This is hard-coded for initial testing purposes.
      char *argList[] = {
	  (char*)("Xyce"),
          (char*)("runXyceWithDAC.cir")
      };
      int argc = sizeof(argList)/sizeof(argList[0]);
      char** argv = argList;

      // Open and initialize the N_CIR_Xyce object
      xyce_open(p);
      xyce_initialize(p,argc,argv);

      // get Device, ADC and DAC names;
      int numDevNames, numDACnames;

      // Initialize array of char array
      const int nameLength = 1000;
      int dataCount = 1000;

      // Initialize arrays of char arrays
      char ** devNames;
      char ** DACnames;
      //devNames = new char*[dataCount];
      devNames = (char **) malloc( dataCount * sizeof(char*));
      int i;
      for (i = 0; i < dataCount; i++)
      {
        //devNames[i] = new char[nameLength];
        devNames[i] = (char *) malloc( nameLength*sizeof(char) );
      }

      //DACnames = new char*[dataCount];
      DACnames = (char **) malloc( dataCount*sizeof(char*) );
      for (i = 0; i < dataCount; i++)
      {
        // DACnames[i] = new char[nameLength];
        DACnames[i] = (char *) malloc( nameLength*sizeof(char) );
      }

      xyce_getDeviceNames(p, (char *)"YADC", &numDevNames, devNames);
      printf("Found %d Memristor devices\n",numDevNames);
      for (i = 0; i < numDevNames; i++)
      {
	printf("Device Name 1: %s\n",devNames[i]);
      }

      xyce_getDACDeviceNames(p, &numDACnames, DACnames);
      printf("Found %d DACs\n",numDACnames);
      for (i = 0; i < numDACnames; i++)
      {
	printf("DAC Name 1: %s\n",DACnames[i]);
      }

      // A bug in the DAC device (put there for Habinero support) only takes
      // the last time in the time voltage pairs list if the current sim time is 0.
      // So simulate a bit first.
      int status;
      double requested_time = 1.0e-10;
      double actual_time, value;
   
      status = xyce_simulateUntil(p, requested_time, &actual_time );
      printf( "Return status from simulateUntil = %d and actual_time = %f\n",status, actual_time);

      double voltageArray [9] = { 0.0, 0.0, 3.0, 3.0, 0.0, 0.0, 3.0, 3.0, 0.0 }; 
      double timeArrayBase [9] = { 0.0, 0.1e-4, 0.2e-4, 0.4e-4, 0.5e-4, 0.7e-4, 0.8e-4, 1.0e-4, 1.1e-4 };
      double timeArray [9];
      for (i=0; i< 9; i++)
      { 
        timeArray[i] = timeArrayBase[i];
      }
      
      double total_sim_time = 20.0e-4;
      for (i=0; i<=10; i++)
      { 
        xyce_updateTimeVoltagePairs( p, DACnames[0], 9, timeArray, voltageArray );
        requested_time = 0.0 + (i+1) * 0.1 * total_sim_time;
        printf( "Calling simulateUntil for requested time %f\n", requested_time );
        status = xyce_simulateUntil(p, requested_time, &actual_time );
        printf( "Return status from simulateUntil = %d and actual_time = %f\n",status, actual_time);
  
	// get some result from the ciruit
        xyce_obtainResponse(p,(char *)"YMEMRISTORRES",&value);
        printf( "Result = %f\n", value);

        // update timeArray to repeat pulse 
        int j;
        for (j=0; j<9; j++)
	{
          timeArray[j] = timeArrayBase[j] + requested_time;
        }
      }

      // Finalize and close the N_CIR_Xyce object
      xyce_close(p);

      // pointer clean-up
      free( p );
      for (i = 0; i < dataCount; i++)
      {
        free( devNames[i] );
      }
      free( devNames );

      for (i = 0; i < dataCount; i++)
      {
        free( DACnames[i] );
      }
      free( DACnames );

      vpi_printf("Exiting calltf for runXyceWithDAC\n"); 
      return 0;
}

void runXyceWithDAC_register()
{
      s_vpi_systf_data tf_data;

      tf_data.type      = vpiSysTask;
      tf_data.tfname    = "$runXyceWithDAC";
      tf_data.calltf    = runXyceWithDAC_calltf;
      tf_data.compiletf = runXyceWithDAC_compiletf;
      tf_data.sizetf    = 0;
      tf_data.user_data = 0;
      vpi_register_systf(&tf_data);
}

void (*vlog_startup_routines[])() = {
    runXyceWithDAC_register,
    0 /* final entry must be zero */
};

