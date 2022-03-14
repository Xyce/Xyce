//-------------------------------------------------------------------------
//   Copyright 2002-2020 National Technology & Engineering Solutions of
//   Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
//   NTESS, the U.S. Government retains certain rights in this software.
//
//   This file is part of the Xyce(TM) Parallel Electrical Simulator.
//
//   Xyce(TM) is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   Xyce(TM) is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with Xyce(TM).
//   If not, see <http://www.gnu.org/licenses/>.
//-------------------------------------------------------------------------


#include "xyce2nLv_sfunction.h"

#include "Xyce_config.h"
#include "N_CIR_Xyce.h"
#include "N_CIR_SecondLevelSimulator.h"
#include "N_DEV_ExternalSimulationData.h"
#include "N_TIA_TwoLevelError.h"

#define S_FUNCTION_LEVEL 2
#define S_FUNCTION_NAME  xyce2nLv_sfunction

// Need to include simstruc.h for the definition of the SimStruct and
// its associated macro definitions.
#include "simstruc.h"

// needed to process data between Matlab and xyce
// This is a Matlab data type used by Matlab
#include "matrix.h"

#include <sstream>
#include <string>
#include <vector>
#include <map>


#define IS_PARAM_DOUBLE(pVal) (mxIsNumeric(pVal) && !mxIsLogical(pVal) &&\
!mxIsEmpty(pVal) && !mxIsSparse(pVal) && !mxIsComplex(pVal) && mxIsDouble(pVal))

// Reserve place for C++ objects and data in Matlab's Work array.
// This enum avoids hardcoding indicies 
// 0 -> pointer to the Xyce object
// 1 -> pointer to current time step.
// 3 -> pointer to the current 
enum WorkArray:int  {XycePtr, CurrentTimeStepPtr, NumberVConnectsPtr, TwoLevelDataPtr, WorkSize};

//
// class to hold the data required for control of inner Xyce
// problem from static functions.  This class is only needed by
// the static functions in this file, so it will be local to this file.
//
class TwoLevelData
{
  public:
    TwoLevelData()
    {
      InnerProblemSize = 0;
      voltageMapKeys = new std::vector<std::string>();
      voltageInputMap = new std::map<std::string,double>();
      outputVector = new std::vector<double>();
      innerJacobianStamp = new std::vector< std::vector<double> >();
      dcOpDone = false;
      initJctFlag = false;
      tlError = new Xyce::TimeIntg::TwoLevelError();
      extSimData = new Xyce::Device::ExternalSimulationData();
    };
    
    ~TwoLevelData()
    {
      delete voltageInputMap;
      delete voltageMapKeys;
      delete outputVector;
      delete innerJacobianStamp;
      delete tlError;
      delete extSimData;
    };
    
    void setUpInputsAndOutputs( int numIn, int numOut)
    {
      
      for( int i=0; i< numIn; i++)
      {
        std::stringstream keystream;
        keystream.clear();
        keystream << "vconnect";
        keystream.width(4);
        keystream.fill('0');
        keystream <<  i;
        std::string stringRep( keystream.str());
        
        mexPrintf(" key is \"%s\"\n", stringRep.c_str());
        voltageMapKeys->push_back(stringRep);
        voltageInputMap->emplace(stringRep, 0.0);
      }
      // set up output ports 
      outputVector->resize(numOut,0.0);
   
      // set up jacobian 
      innerJacobianStamp->resize(numOut);
      for( int i=0; i<numIn; i++)
      {
        innerJacobianStamp->at(i).resize(numIn,0.0);
      }
    };
    
    int InnerProblemSize;
    std::vector<std::string> *voltageMapKeys;
    std::map<std::string,double> *voltageInputMap;
    std::vector<double> *outputVector;
    std::vector< std::vector<double> > *innerJacobianStamp;

    Xyce::TimeIntg::TwoLevelError *tlError;  // data structure related to global error control
    Xyce::Device::ExternalSimulationData *extSimData;  // data on time step control
    bool dcOpDone;
    bool initJctFlag;
    
    

};


//
// Function: mdlInitializeSizes ===============================================
//
// Abstract:
//    The sizes information is used by Simulink to determine the S-function
//    block's characteristics (number of inputs, outputs, states, etc.).
//
static void mdlInitializeSizes(SimStruct *S)
{
  // No expected parameters
  ssSetNumSFcnParams(S, 1);

  // Parameter mismatch will be reported by Simulink
  if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) 
  {
    return;
  }
  
  // Specify I/O
  const int numIOPorts = 2;
  if (!ssSetNumInputPorts(S, numIOPorts)) 
  {
    return;
  }
  for( int i=0; i<numIOPorts; i++)
  {
    ssSetInputPortWidth(S, i, DYNAMICALLY_SIZED);
    ssSetInputPortDirectFeedThrough(S, i, 1);
  }
  if (!ssSetNumOutputPorts(S, numIOPorts))
  { 
    return;
  }
  for( int i=0; i<numIOPorts; i++)
  {
    ssSetOutputPortWidth(S, i, DYNAMICALLY_SIZED);
  }
  
  ssSetNumSampleTimes(S, 1);

  // allocate space in Matlab's work array
  // this is where instance specific data must be stored as 
  // all of the functions that interact with Simulink are static.
  ssSetNumPWork(S, WorkSize);
  
  // Simulink's terminology for Operating point is any valid solution 
  // at a given time.  It's how Simulink handled roll back to earlier times
  // by just giving a model a new operating point from which to start.
  ssSetOperatingPointCompliance(S, USE_CUSTOM_OPERATING_POINT);

  /* Set this S-function as runtime thread-safe for multicore execution */
  ssSetRuntimeThreadSafetyCompliance(S, RUNTIME_THREAD_SAFETY_COMPLIANCE_TRUE);

  ssSetOptions(S,
         SS_OPTION_WORKS_WITH_CODE_REUSE |
         SS_OPTION_EXCEPTION_FREE_CODE |
         SS_OPTION_DISALLOW_CONSTANT_SAMPLE_TIME);
}


//
// Function: mdlInitializeSampleTimes =========================================
// Abstract:
//   This function is used to specify the sample time(s) for your
//   S-function. You must register the same number of sample times as
//   specified in ssSetNumSampleTimes.
//
static void mdlInitializeSampleTimes(SimStruct *S)
{
  ssSetSampleTime(S, 0, INHERITED_SAMPLE_TIME);
  ssSetOffsetTime(S, 0, 0.0);
  ssSetModelReferenceSampleTimeDefaultInheritance(S); 
}

//
// Function: mdlStart =======================================================
// Abstract:
//   This function is called once at start of model execution. If you
//   have states that should be initialized once, this is the place
//   to do it.
//
#define MDL_START
static void mdlStart(SimStruct *S)
{
  // Store an N_CIR_Xyce object in the workspace
  Xyce::Circuit::SecondLevelSimulator *xyce = new Xyce::Circuit::SecondLevelSimulator();
  ssGetPWork(S)[XycePtr] = xyce;

  // get pointer to first param
  const mxArray * firstParamPtr = ssGetSFcnParam(S, 0);
  int numElements = mxGetNumberOfElements(firstParamPtr);
  int numFields = mxGetNumberOfFields(firstParamPtr);
  // just for debugging to ensure read of circuit name from parameter list worked.
  //mexPrintf("NumElements = %i\n", numElements);
  //mexPrintf("NumChars = %i\n", numFields);
  char * tempBuf = new char[numElements+1];

  mxGetString(firstParamPtr, tempBuf, numElements+1);

  std::string theFilename( tempBuf );
  mexPrintf( "Input Filename is \"%s\"\n", theFilename.c_str() );

  // initialize Xyce with this circuit 
  char ** tempCBuf = new char *[3];
  tempCBuf[0] = const_cast<char *>(std::string("Xyce").c_str());
  tempCBuf[1] = const_cast<char *>(theFilename.c_str());
  tempCBuf[2] = 0;
  
  int retVal = xyce->initialize(2, tempCBuf);
  if ( retVal != Xyce::Circuit::Simulator::RunStatus::SUCCESS )
  {
    mexPrintf("Xyce failed on initialization. Return Code: %i\n", retVal );
  }
  xyce->startupSolvers();

  // clean up temporary allocated storage.
  delete [] tempCBuf;
  
  // allocate space for the number of vconnect ports
  int numVConnects = 0;
  int * nVCptr = new int;
  ssGetPWork(S)[NumberVConnectsPtr] = nVCptr;
  
  // query Xyce for number of "vconnect" devices
  const std::string baseName("v");
  std::vector<std::string> potentialConnectDevices;
  xyce->getDeviceNames(baseName, potentialConnectDevices);
  if( potentialConnectDevices.empty())
  {
    mexPrintf("No v-source devices in the circuit.\n");
  }
  else
  {
    for( int i=0; i<potentialConnectDevices.size(); i++)
    {
      if(potentialConnectDevices[i].find("VCONNECT") != std::string::npos )
      {
        numVConnects++;
      }
    }
  }
  
  // save the number of vconnect ports for later use
  *nVCptr = numVConnects;

  // set up storage of Xyce runtime data
  TwoLevelData * xyceSimData = new TwoLevelData();
  int_T numberOfInputs = ssGetNumInputPorts(S);
  int_T numberOfOutputs = ssGetNumOutputPorts(S);
  if( (numberOfInputs > numVConnects) )
  {
    mexPrintf("Warning: There are more Simulink inputs than VCONNECT ports in the Xyce circuit.\n Inputs after %d will be ignored.\n", numVConnects);
  }
  xyceSimData->setUpInputsAndOutputs( static_cast<int>(numVConnects), static_cast<int>(numVConnects));
  
  // get data on transient simulation from Xyce
  //Xyce::Device::ExternalSimulationData  xyceSimData;

  // this has Xyce populate the external simulator data structure.
  xyce->endTimeStep(*(xyceSimData->extSimData));
  
  // this is the final time Simulink expects.
  time_T simulinkMaxTime = ssGetTFinal(S);
  mexPrintf("Simulink Max time = %g, Xyce Max time = %g", simulinkMaxTime, xyceSimData->extSimData->finalTime);
  
  // store the sim data
  ssGetPWork(S)[TwoLevelDataPtr] = xyceSimData;

  // Store Xyce simulaiton time so we can tell if we need 
  // to advance time
  double * XyceSimTime = new double;
  *XyceSimTime = 0.0;
  ssGetPWork(S)[CurrentTimeStepPtr] = XyceSimTime;
  // Store new C++ object in the pointers vector
  //DoubleAdder *da  = new DoubleAdder();
  //ssGetPWork(S)[1] = da;
}


//
// Function: mdlOutputs =======================================================
// Abstract:
//   In this function, you compute the outputs of your S-function
//   block.
//
static void mdlOutputs(SimStruct *S, int_T tid)
{  
  mexPrintf("Entering mdlOutputs \n");
  // Get data addresses of I/O
  InputRealPtrsType  u = ssGetInputPortRealSignalPtrs(S,0);
  int_T numberOfInputs = ssGetNumInputPorts(S);
  mexPrintf("In mdlOutputs with input ports at %d\n", static_cast<int>(numberOfInputs));
  real_T *y = ssGetOutputPortRealSignal(S, 0);
  int_T numberOfOutputs = ssGetNumOutputPorts(S);

  // get instance data from the work space
  Xyce::Circuit::SecondLevelSimulator *xyce = static_cast<Xyce::Circuit::SecondLevelSimulator *>(ssGetPWork(S)[XycePtr]);
  double * XyceSimTime = static_cast<double *>(ssGetPWork(S)[CurrentTimeStepPtr]);
  TwoLevelData * xyceSimData = static_cast<TwoLevelData *>(ssGetPWork(S)[TwoLevelDataPtr]);
  int * nVCptr = static_cast<int *>(ssGetPWork(S)[NumberVConnectsPtr]);
  
  // update inputs from Simulink to Xyce
  for( int i=0; i<numberOfInputs; i++)
  {
    mexPrintf("Setting xyce node \"%s\" to %g\n", (xyceSimData->voltageMapKeys->at(i)).c_str(), (*u[i]));
    if( i < *nVCptr)
    {
      xyceSimData->voltageInputMap->at( xyceSimData->voltageMapKeys->at(i))  = *u[i];
    }
  }
  // There is a concept of system time which is the global time that
  // Simulink is working at and task time which can be different from Simulink's time
  // if a given model is implemented to take a task time step.  I'm not sure 
  time_T systemTime = ssGetT(S);
  time_T taskTime = ssGetTaskTime(S, tid);
  time_T finalTime = ssGetTFinal(S);
  
  
  if( !(xyce->simulationComplete())  )
  { 
    // mexPrintf("simulationComplete() was false\n");
    // first check if we need to do a DC operating point in Xyce 
    if( !(xyceSimData->dcOpDone))
    {
      mexPrintf("Doing DC OP calculation.\n");
      // set up data structure that controls Xyce for a DC operating 
      // point calculation that precedes a transient simulation.
      xyceSimData->initJctFlag = true;
      xyceSimData->extSimData->is_transient = true;
      xyceSimData->extSimData->current_time = 0.0;
      xyceSimData->extSimData->final_time = finalTime;
      xyceSimData->extSimData->current_time_step_size = xyceSimData->extSimData->nextTimeStep;
      xyceSimData->extSimData->previous_time_step_size = xyceSimData->extSimData->currTimeStep;
      xyceSimData->extSimData->time_step_number = 0;
      xyceSimData->extSimData->forceOrder=true;
      xyceSimData->extSimData->imposedTimeIntegrationOrder=1;
      xyceSimData->extSimData->forceBeginningIntegration=true;
      xyceSimData->extSimData->imposedBeginningIntegration=true;
      
      bool xyceReturnFlag = false;
      xyceReturnFlag = xyce->startTimeStep(*(xyceSimData->extSimData));
      
      // set up inputs to the inner problem
      
      xyceReturnFlag = xyce->simulateStep(xyceSimData->initJctFlag, *(xyceSimData->voltageInputMap), *(xyceSimData->outputVector), *(xyceSimData->innerJacobianStamp), *(xyceSimData->tlError));
      if (!xyceReturnFlag)
      {
        // may need to try homotopy at this point if Xyce doesn't do that
        // automatically.  Need to verify this.
        mexPrintf("Xyce problem failed to solve DC operating point.  Exiting.\n");
        exit(0);
      }
      // check convergence 
      xyce->stepSuccess(Xyce::Analysis::TWO_LEVEL_MODE_TRANSIENT_DCOP);
      xyce->endTimeStep(*(xyceSimData->extSimData));
      xyceSimData->dcOpDone = true;
    }
    else
    {
      mexPrintf("Taking Transient Step with simulink time %g\n",systemTime );
      // Next have Xyce take as many small time steps as needed to 
      // get to the time point requested by Simulink 
      
      // is this still true AFTER we step out of the DC OP?
      xyceSimData->extSimData->beginIntegrationFlag=true;
      int iTimeStep=0;
      bool timeSteppingFinished = false;
      // set maximum time step for Xyce to be the difference between what 
      // Simulink wants to take and where Xyce currently is 
      //double maxTimeStepForXyce = (0.1*(systemTime - (*XyceSimTime)) > 0.1) ?  0.1*(systemTime - (*XyceSimTime)): 0.1;
      double maxTimeStepForXyce = xyceSimData->extSimData->current_time_step_size;
      // time loop
      while(!timeSteppingFinished)
      {
        // before resetting the sim data, check against maxTimeStepForXyce.
        if ( xyceSimData->extSimData->nextTimeStep > systemTime)
        {
          xyceSimData->extSimData->nextTime -= xyceSimData->extSimData->nextTimeStep;
          xyceSimData->extSimData->nextTime += maxTimeStepForXyce;
          xyceSimData->extSimData->nextTimeStep = maxTimeStepForXyce;
        }

        //resetSimData(extSimData);
        xyceSimData->extSimData->is_transient = true;
        xyceSimData->extSimData->current_time = xyceSimData->extSimData->nextTime;
        xyceSimData->extSimData->final_time = systemTime;
        xyceSimData->extSimData->current_time_step_size = maxTimeStepForXyce;
        xyceSimData->extSimData->previous_time_step_size = maxTimeStepForXyce;
        xyceSimData->extSimData->time_step_number = xyceSimData->extSimData->timeStepNumber;
        xyceSimData->extSimData->forceOrder=true;
        xyceSimData->extSimData->imposedTimeIntegrationOrder= xyceSimData->extSimData->currentOrder;
        xyceSimData->extSimData->forceBeginningIntegration=true;
        xyceSimData->extSimData->imposedBeginningIntegration = xyceSimData->extSimData->beginIntegrationFlag;
        
        //mexPrintf("Begin Time step (outer solve) # %d, current, next, step, final = %g, %g, %g, %g\n", iTimeStep+1 , 
        //  xyceSimData->extSimData->current_time, xyceSimData->extSimData->nextTime, xyceSimData->extSimData->nextTimeStep, xyceSimData->extSimData->final_time);
        
        xyceSimData->initJctFlag = false;
        
        // update inputs from Simulink 
        // update the sinewave source:
        // double time = (extSimData.nextTime);
        // vconnect000 = V0 + VA * std::sin(2.0*mpi*(FREQ*time));
        bool xyceReturnFlag = false;
        xyceReturnFlag = xyce->startTimeStep(*(xyceSimData->extSimData));
        
        // set up inputs to the inner problem
        // xyceReturnFlag = newtonSolve();
        xyceReturnFlag = xyce->simulateStep(xyceSimData->initJctFlag, *(xyceSimData->voltageInputMap), *(xyceSimData->outputVector), *(xyceSimData->innerJacobianStamp), *(xyceSimData->tlError));
        if (!xyceReturnFlag)
        {
          // may need to try homotopy at this point if Xyce doesn't do that
          // automatically.  Need to verify this.
          mexPrintf("Xyce problem failed to solve in transient step at time %g.\n", (xyceSimData->extSimData->nextTime));
          exit(0);
        }

        if (xyceReturnFlag)
        {
          // process "successful" time step
          xyce->stepSuccess(Xyce::Analysis::TWO_LEVEL_MODE_TRANSIENT); 
        }
        else
        {
          Xyce::lout() << "Time step failed to converge. Exiting" << std::endl;
          exit(0);
        }

        xyce->endTimeStep(*(xyceSimData->extSimData));

        // check if we have reached the time step Simulink requires
        
        if ((xyceSimData->extSimData->current_time) >= (xyceSimData->extSimData->final_time)) 
        { 
          timeSteppingFinished = true; 
          mexPrintf( "==>Finishing time step: current - final = %g \n", ( (xyceSimData->extSimData->current_time) - (xyceSimData->extSimData->final_time)) );
        }

        //Xyce::lout() << extSimData;

        xyceSimData->extSimData->beginIntegrationFlag=false;
        ++iTimeStep;
      }
      // set XyceSimTime to where Xyce reported leaving off.
      *XyceSimTime = xyceSimData->extSimData->final_time;
      mexPrintf("Finished a time step with Xyce's final time at %g and finalTime at %g\n", 
        *XyceSimTime, xyceSimData->extSimData->finalTime);
    }
    

  }
  else
  {
    mexPrintf( "Simulink time, %g, is greater than Xyce's simulation time, %g.\n", systemTime, *XyceSimTime);
    mexPrintf( "Xyce will no longer update!  Please update simulation time in the Xyce input file.\n");
  }
  
  
  // Call AddTo method and return peak value
  for( int i=0; i<numberOfOutputs; i++)
  {
    mexPrintf( "Returning %d value %g\n", i, (xyceSimData->outputVector->at(i)));
    if( i < *nVCptr) 
    {
      y[i] = xyceSimData->outputVector->at(i);
    }
    else
    {
      y[i] = 0.0;
    }
  }
  
}

#ifdef MATLAB_MEX_FILE
/* Define to indicate that this S-Function has the mdlG[S]etOperatingPoint methods */
#define MDL_OPERATING_POINT

//
// Function: mdlGetOperatingPoint ==================================================
// Abstract:
//    Save the operating point of this block and return it to Simulink 
//
static mxArray* mdlGetOperatingPoint(SimStruct* S)
{
  //DoubleAdder *da = static_cast<DoubleAdder*>(ssGetPWork(S)[1]);
  //return mxCreateDoubleScalar(da->GetPeak());
  return mxCreateDoubleScalar(0.0);
}


//
// Function: mdlSetOperatingPoint =================================================
// Abstract:
//   Restore the operating point of this block based on the provided data (ma)
//    The data was saved by mdlGetOperatingPoint
//
static void mdlSetOperatingPoint(SimStruct* S, const mxArray* ma)
{
  // Retrieve C++ object from the pointers vector
  //DoubleAdder *da = static_cast<DoubleAdder*>(ssGetPWork(S)[1]);
  //da->SetPeak(mxGetPr(ma)[0]);
}
#endif // MATLAB_MEX_FILE

//
// Function: mdlTerminate =====================================================
// Abstract:
//   In this function, you should perform any actions that are necessary
//   at the termination of a simulation.  For example, if memory was
//   allocated in mdlStart, this is the place to free it.
//
static void mdlTerminate(SimStruct *S)
{
  Xyce::Circuit::SecondLevelSimulator *xyce = static_cast<Xyce::Circuit::SecondLevelSimulator *>(ssGetPWork(S)[XycePtr]);
  xyce->finishSolvers();
  int retVal = xyce->finalize();
  if( retVal != Xyce::Circuit::Simulator::RunStatus::SUCCESS )
  {
	mexPrintf("ERROR: Xyce failed to close without error, code=%d\n", retVal);
  }
  
  delete xyce;
  // Retrieve and destroy C++ object
  //DoubleAdder *da = static_cast<DoubleAdder *>(ssGetPWork(S)[1]);
  //delete da;
  double * XyceSimTime = static_cast<double *>(ssGetPWork(S)[CurrentTimeStepPtr]);
  delete XyceSimTime;
  
  // deallocate number of connections 
  int * nVCptr  = static_cast<int *>(ssGetPWork(S)[NumberVConnectsPtr]);
  delete nVCptr;
}


// Required S-function trailer
#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif
