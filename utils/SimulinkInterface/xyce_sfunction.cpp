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


#include "xyce_sfunction.h"

#include "Xyce_config.h"
#include "N_CIR_MixedSignalSimulator.h"


#define S_FUNCTION_LEVEL 2
#define S_FUNCTION_NAME  xyce_sfunction

// Need to include simstruc.h for the definition of the SimStruct and
// its associated macro definitions.
#include "simstruc.h"

// needed to process data between Matlab and xyce
// This is a Matlab data type used by Matlab
#include "matrix.h"

#include <utility>
#include <fstream>
#include <sstream>

#define IS_PARAM_DOUBLE(pVal) (mxIsNumeric(pVal) && !mxIsLogical(pVal) &&\
!mxIsEmpty(pVal) && !mxIsSparse(pVal) && !mxIsComplex(pVal) && mxIsDouble(pVal))

// Reserve place for C++ objects and data in Matlab's Work array.
// This enum avoids hardcoding indicies 
// 0 -> pointer to the Xyce object
// 1 -> pointer to current time step. 
enum WorkArray:int  {InputNamesVecPtr, OutputNamesVecPtr, XycePtr, CurrentTimeStepPtr, DACVecPtr, ADCMapPtr, WorkSize};

// 
// Utility functions.  
// Matlab and Simulink rely on C character arrays for strings.  These are utility functions
// to help convert char * [] to strings and to parse up a string into sub-strings
//
std::string ConvertMxArrayToString( const mxArray * aSimulinkCharArray)
{
  int numElements = mxGetNumberOfElements(aSimulinkCharArray);
  char * tempBuf = new char[numElements+1];
  mxGetString(aSimulinkCharArray, tempBuf, numElements+1);
  std::string fullStringValue( tempBuf );
  return fullStringValue;
}

std::vector<std::string> ParseCommaSeparatedString( const std::string inputString)
{
  std::vector<std::string> returnVector;
  std::stringstream inputAsStream( inputString );
  while( inputAsStream.good())
  {
    std::string aSubString;
    getline( inputAsStream, aSubString, ',' );
    returnVector.push_back( aSubString);
  }
  return returnVector;
}



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
  ssSetNumSFcnParams(S, 3);

  // Parameter mismatch will be reported by Simulink
  if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) 
  {
    return;
  }
  
  // Specify I/O
  if (!ssSetNumInputPorts(S, 1))
  {
    return;
  }
  ssSetInputPortWidth(S, 0, DYNAMICALLY_SIZED);
  
  // input ports are used on calculations of outputs
  ssSetInputPortDirectFeedThrough(S, 0, 1);
  
  if (!ssSetNumOutputPorts(S,1)) 
  {
    return;
  }
  ssSetOutputPortWidth(S, 0, DYNAMICALLY_SIZED);

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
  Xyce::Circuit::MixedSignalSimulator *xyce = new Xyce::Circuit::MixedSignalSimulator();
  ssGetPWork(S)[XycePtr] = xyce;

  // get pointer to first param
  const mxArray * firstParamPtr = ssGetSFcnParam(S, 0);
  std::string theFilename = ConvertMxArrayToString(firstParamPtr); 
  mexPrintf( "Input Filename is \"%s\"\n", theFilename.c_str() );
  
  const mxArray * inputNamesParamPtr = ssGetSFcnParam(S, 1);
  std::string inputNames = ConvertMxArrayToString(inputNamesParamPtr); 
  std::vector<std::string> inputNamesVec = ParseCommaSeparatedString(inputNames);
  std::vector<std::string> * inputNamesVecPtr = new std::vector<std::string>(inputNamesVec);
  ssGetPWork(S)[InputNamesVecPtr] = inputNamesVecPtr;
  mexPrintf( "Input names num elements %s\n", inputNames.c_str());
  for( auto i=0; i<inputNamesVecPtr->size(); i++)
  {
    mexPrintf("%i, %s\n", i, inputNamesVecPtr->at(i).c_str());
  }
  
  const mxArray * outputNamesParamPtr = ssGetSFcnParam(S, 2);
  std::string outputNames = ConvertMxArrayToString(outputNamesParamPtr); 
  std::vector<std::string> outputNamesVec = ParseCommaSeparatedString(outputNames);
  std::vector<std::string> * outputNamesVecPtr = new std::vector<std::string>(outputNamesVec);
  ssGetPWork(S)[OutputNamesVecPtr] = outputNamesVecPtr;
  mexPrintf( "Output names num elements %s\n", outputNames.c_str());
  for( auto i=0; i<outputNamesVecPtr->size(); i++)
  {
    mexPrintf("%i, %s\n", i, outputNamesVecPtr->at(i).c_str());
  }
  
  // check if file exists because Xyce will exit if it does not.
  std::ifstream theInputFileStream( theFilename.c_str() );
  if( !theInputFileStream.good())
  {
    // report error
    static char errorString[256];
    sprintf( errorString, "Input file, \"%s\", does not exist.  Ending simulation.",  theFilename.c_str());
    ssSetLocalErrorStatus(S, errorString);
    mexPrintf( errorString );
    ssSetStopRequested( S, 1);
    return;
  }

  // initialize Xyce with this circuit 
  char ** tempCBuf = new char *[2];
  tempCBuf[0] = const_cast<char *>(std::string("Xyce").c_str());
  tempCBuf[1] = const_cast<char *>(theFilename.c_str());
  int retVal = xyce->initialize(2, tempCBuf);
  if ( retVal != Xyce::Circuit::Simulator::RunStatus::SUCCESS )
  {
    static char errorString[256];
    sprintf( errorString, "Xyce failed on initialization. Return Code: %i\n", retVal );
    ssSetLocalErrorStatus(S, errorString);
    mexPrintf( errorString );
    ssSetStopRequested( S, 1);
  }

  // clean up temporary allocated storage.
  delete [] tempCBuf;

  // Store Xyce simulaiton time so we can tell if we need 
  // to advance time
  double * XyceSimTime = new double;
  *XyceSimTime = 0.0;
  ssGetPWork(S)[CurrentTimeStepPtr] = XyceSimTime;
  
  std::vector< std::string > * dacNamesPtr = new std::vector< std::string >();
  ssGetPWork(S)[DACVecPtr] = dacNamesPtr;
  
  // populate the map from Xyce 
  xyce->getDACDeviceNames(*dacNamesPtr);
  if( dacNamesPtr->empty())
  {
    mexPrintf("No DAC devcies found in the circuit.\n");
  }

  // check that input names match a name from Xyce.
  for( auto i=0; i<inputNamesVec.size(); i++)
  {
    // locate the Xyce device name that matches inputNamesVec[i]
    bool found=false;
    for( auto j=0; (j<dacNamesPtr->size() && !found); j++)
    {
      if( inputNamesVec[i].compare(dacNamesPtr->at(j) ) )
      {
        mexPrintf( "Found matching input name %s\n", inputNamesVec[i].c_str());
        found=true;
      }
    }
    if( !found )
    {
      mexPrintf("Input name %s was not found in the Xyce netlist\n", inputNamesVec[i].c_str());
    }
  }
  
  std::map<std::string,std::map<std::string,double> > * adcStrMapPtr = new  std::map<std::string,std::map<std::string,double> >();
  ssGetPWork(S)[ADCMapPtr] = adcStrMapPtr;
  // populate the ADC map from xyce
  xyce->getADCMap( *adcStrMapPtr );
  if( adcStrMapPtr->empty() )
  {
    mexPrintf("No ADC devices found in the circuit.\n");
  }
  // set up a map from the Simulink outputNamesVec index to the names as ordered by Xyce 
  
  for( auto i=0; i<outputNamesVec.size(); i++)
  {
    // locate the Xyce device name that matches inputNamesVec[i]
    bool found=false;
    for( auto mapItr=adcStrMapPtr->begin(); ((mapItr!=adcStrMapPtr->end()) && !found); mapItr++)
    {
      if( outputNamesVec[i].compare(mapItr->first ) )
      {
        mexPrintf( "Found matching output name %s\n", outputNamesVec[i].c_str());
        found=true;
      }
    }
    if( !found )
    {
      mexPrintf("Output name %s was not found in the Xyce netlist\n", outputNamesVec[i].c_str());
    }
  }
  
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
  // Retrieve C++ object from the pointers vector
  //DoubleAdder *da = static_cast<DoubleAdder *>(ssGetPWork(S)[1]);
  
  int_T nInputPorts = ssGetNumInputPorts(S);
  int_T nPortWidth = 0;
  if( nInputPorts == 1)
  {
    for( int i=0; i< nInputPorts; i++)
    {
      nPortWidth = ssGetInputPortWidth(S,i);
      mexPrintf("InputPort has a width of %i\n", nPortWidth);
    }
  }
  
  // Get data addresses of I/O
  InputRealPtrsType  u = ssGetInputPortRealSignalPtrs(S,0);
  //double * u = static_cast<double *>(const_cast<void *>(ssGetInputPortSignal(S,0)));
  
  //real_T *y = ssGetOutputPortRealSignal(S, 0);
  real_T * y = (real_T *)ssGetOutputPortSignal(S,0);
  
  int_T numberOfOutputs = ssGetNumOutputPorts(S); 
  int_T outputPortWidth = 0;
  if( numberOfOutputs == 1) 
  {
    for( int i=0; i< numberOfOutputs; i++)
    {
      outputPortWidth = ssGetOutputPortWidth(S,i);
      mexPrintf("OutputPort has a width of %i\n", outputPortWidth);
    }
  }

  double * XyceSimTime = static_cast<double *>(ssGetPWork(S)[CurrentTimeStepPtr]);  
  Xyce::Circuit::MixedSignalSimulator *xyce = static_cast<Xyce::Circuit::MixedSignalSimulator *>(ssGetPWork(S)[XycePtr]);
  // There is a concept of system time which is the global time that
  // Simulink is working at and task time which can be different from Simulink's time
  // if a given model is implemented to take a task time step.  I'm not sure 
  time_T systemTime = ssGetT(S);
  time_T taskTime = ssGetTaskTime(S, tid);
  
  // update inputs from Simulink to Xyce
  std::vector< std::string > * dacNamesPtr = static_cast<std::vector< std::string > *>(ssGetPWork(S)[DACVecPtr]);
  
  if( !dacNamesPtr->empty())
  {
    std::map< std::string, std::vector< std::pair<double,double> >* > timeVoltageUpdateMap;
    int numDacs = dacNamesPtr->size();
    
    int_T simulinkInput = 0;
    
    for( auto dacNamesItr = dacNamesPtr->begin(); dacNamesItr != dacNamesPtr->end(); dacNamesItr++ )
    {
      std::vector< std::pair<double,double> > *tvpVec = new std::vector< std::pair<double,double> >();  
      if ( simulinkInput < nPortWidth)
      {
        double val = *(u[simulinkInput++]);
        tvpVec->push_back( std::pair<double,double>(systemTime, val ));
        mexPrintf("DAC %s time %g  value %g\n", (*dacNamesItr).c_str(), systemTime, val );
      }
      else
      {
        tvpVec->push_back( std::pair<double,double>(systemTime, 0.0 ));
        mexPrintf("DAC %s time %g  value %g\n", (*dacNamesItr).c_str(), systemTime, 0.0 );
      }
      timeVoltageUpdateMap.emplace( *dacNamesItr, tvpVec);
      
    }
    
    bool retVal = xyce->updateTimeVoltagePairs(timeVoltageUpdateMap);  
    if( !retVal )
    {
      static char errorString[256];
      sprintf( errorString, "Xyce returned an error on call to updateTimeVoltagePairs at time %g\n", systemTime );
      ssSetLocalErrorStatus(S, errorString);
      mexPrintf( errorString );
      ssSetStopRequested( S, 1);
    }
    
    //
    // need to delete old tvpVec objects held in timeVoltageUpdateMap;
    //
    auto mapItr = timeVoltageUpdateMap.begin();
    auto endItr = timeVoltageUpdateMap.end();
    while( mapItr != endItr)
    {
      // delete the std::vector< std::pair<double,double> >* object allocated for this entry
      delete (*mapItr).second;
      // move on to the next entry
      mapItr++;
    }

    
  }
  
  
  // needs two corrections
  // 1. use take provisional step in the MixedSignalSimulator interface
  // 2. use output from that step to set values on Simulink output
  
  // storage for updates from Xyce
  std::map< std::string, std::vector< std::pair<double,double> > > timeVoltageUpdateMap;
  
  if( !(xyce->simulationComplete()) )
  {  
    if( systemTime >= *XyceSimTime )
    {
      // Simulink system time is ahead of the Xyce time.
      // have Xyce try to integrate forward in time
 
      bool timeSteppingFinished = false;
      double maxTimeStepForXyce = (0.1*(systemTime - (*XyceSimTime)) > 1e-6) ?  0.1*(systemTime - (*XyceSimTime)): 1e-6;
      double timeStepTaken = 0;
      while(!timeSteppingFinished)
      {
        bool stepStatus = xyce->provisionalStep( maxTimeStepForXyce, timeStepTaken, timeVoltageUpdateMap);
        if( stepStatus )
        {
          // three cases that need to be handled.  
          // 1. Simulink wanted the solution at time=0.  This is just a DC OP for Xyce.
          // 2. Short step was taken by Xyce but Xyce needs to continue so that it gets to simulink's time 
          // 3. There was a change in output on Xyce's DAC's and results should be sent to 
          //    Simulink before we get to the end of the simulink step.  Not sure if Simulink 
          //    can handle this case.
          if( (systemTime == 0) && (*XyceSimTime == 0))
          {   
            // for the DC OP we don't accept the provisional step or advance time
            // but time stepping is finished. 
            timeSteppingFinished = true;
          }
          else
          {
            xyce->acceptProvisionalStep();
            *XyceSimTime = *XyceSimTime + timeStepTaken;
            if( *XyceSimTime >= systemTime)
            {
              timeSteppingFinished = true;
              // horrible hack to see if what happens when Simulink forces Xyce to take the steps it needs.
              if( *XyceSimTime > systemTime)
              {
                *XyceSimTime = systemTime;
              }
            }
          }
        }
      }
    }
    else
    {
      mexPrintf( "Xyce jumped past simulink time. Simulink: %g,  Xyce: %g. diff %g\n", systemTime, *XyceSimTime, (*XyceSimTime-systemTime));
      mexPrintf( "Reusing last time results from Xyce.\n");
      //bool stepStatus = xyce->getTimeVoltagePairs(timeVoltageUpdateMap);
      //if( !stepStatus)
      //{
      //  mexPrintf( "Within simulation, repeat call to getTimeVoltagePairs failed\n");
      //}
    }
  }
  else
  {
    mexPrintf( "Simulink time, %g, is greater than Xyce's simulation time, %g.\n", systemTime, *XyceSimTime);
    mexPrintf( "Xyce will no longer update!  Please update simulation time in the Xyce input file.\n");
    bool stepStatus = xyce->getTimeVoltagePairs(timeVoltageUpdateMap);
    if( !stepStatus)
    {
      mexPrintf( "Repeat call to getTimeVoltagePairs failed\n");
    }
  }

  if( timeVoltageUpdateMap.empty())  
  {
    // no results from Xyce (either there are not any ADC's or Xyce didn't update )
    // just return zeros 
    for( int i=0; i<outputPortWidth; i++)
    {
      y[i] = 0.0;
    }
  }
  else
  {
    int simulinkOutputIndex = 0;
    
    // this loop is for debugging only.
    for(auto tVUpdateMapItr = timeVoltageUpdateMap.begin(); tVUpdateMapItr != timeVoltageUpdateMap.end(); tVUpdateMapItr++ )
    {
      std::vector< std::pair<double,double> >  timeVoltagePairs = tVUpdateMapItr->second;
      std::string devName = tVUpdateMapItr->first;
      // may need better logic here, but for now just grap the last pair 
      if( ! timeVoltagePairs.empty())
      {
        auto lastElementItr = timeVoltagePairs.rbegin();
        double lastTime = lastElementItr->first;
        double lastValue = lastElementItr->second;
        mexPrintf("DEBUG For devcie %s values (%g, %g) in time vec pair\n", devName.c_str(), lastTime, lastValue);
      }
    }
    
    for(auto tVUpdateMapItr = timeVoltageUpdateMap.begin(); tVUpdateMapItr != timeVoltageUpdateMap.end(); tVUpdateMapItr++ )
    {
      if( simulinkOutputIndex < outputPortWidth)
      {
        std::vector< std::pair<double,double> >  timeVoltagePairs = tVUpdateMapItr->second;
        // may need better logic here, but for now just grap the last pair 
        if( timeVoltagePairs.empty())
        {
          y[simulinkOutputIndex] = 0.0;
        }
        else
        {
          auto lastElementItr = timeVoltagePairs.rbegin();
          double lastTime = lastElementItr->first;
          double lastValue = lastElementItr->second;
          mexPrintf("Found values (%g, %g) in time vec pair\n", lastTime, lastValue);
          y[simulinkOutputIndex] = lastValue;
        }
        simulinkOutputIndex++;
      }
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
  Xyce::Circuit::MixedSignalSimulator *xyce = static_cast<Xyce::Circuit::MixedSignalSimulator *>(ssGetPWork(S)[XycePtr]);
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
  
  std::vector< std::string > * dacNamesPtr = static_cast<std::vector< std::string > *>(ssGetPWork(S)[DACVecPtr]);
  delete dacNamesPtr;
   
  std::map<std::string,std::map<std::string,double> > * adcStrMapPtr = static_cast<std::map<std::string,std::map<std::string,double> > *>(ssGetPWork(S)[ADCMapPtr]);
  delete adcStrMapPtr;
  
  
}


// Required S-function trailer
#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif
