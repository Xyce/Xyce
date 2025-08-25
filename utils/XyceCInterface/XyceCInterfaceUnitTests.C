//-------------------------------------------------------------------------
//   Copyright 2002-2025 National Technology & Engineering Solutions of
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

//-------------------------------------------------------------------------
//
// Purpose        : Unit tests for top-level Xyce::Simulator class
//
// Special Notes  :
//
// Creator        : Richard Schiek, SNL, Parallel Computational Sciences
//
// Creation Date  : 2024/2022
//
//
//
//
//-------------------------------------------------------------------------

#include "Xyce_config.h"
#include <N_UTL_Math.h>
#include <N_CIR_XyceCInterface.h>

#include <gtest/gtest.h>
#include <string.h>
#include <math.h>
#if (__cplusplus>=201703L)
#include <filesystem>
#endif

//
// Xyce::Circuit::Simulator functions that need to be tested
//
// x constructor
// x initialize()
// x finalize()
// x getTime()
// x getDACDeviceNames()
// x getADCMap()
// x getTimeVoltagePairs()
// x updateTimeVoltagePairs()
// x simulationComplete()
// x getFinalTime()
// x setCircuitParameter()
// x getCircuitValue()
// x simulateUntil() in stead of provisionalStep() & acceptProvisionalStep()
// 


TEST ( XyceCInterface, Open)
{
  void * xycePtr = NULL;
  xyce_open( & xycePtr);
  EXPECT_TRUE( *((long *)xycePtr) != 0 );
  
}

TEST ( XyceCInterface, OpenAndClose)
{
  void * xycePtr = NULL;
  xyce_open( & xycePtr);
  EXPECT_TRUE( ((long *)xycePtr) != 0 );
  xyce_close( & xycePtr );
  EXPECT_TRUE( ((long *)xycePtr) == 0 );
  
}

TEST ( XyceCInterface, InitializeNetlist1)
{
  void * xycePtr = NULL;
  xyce_open( & xycePtr);
  EXPECT_TRUE( ((long *)xycePtr) != 0 );
  const char argv0[] = "Xyce";
  const char argv1[] = "TestNetlist1.cir"; 
  char * argvarray[2];
  argvarray[0] = (char *) &argv0[0];
  argvarray[1] = (char *) &argv1[0];
  int return_status = xyce_initialize( & xycePtr, 2, argvarray);
  EXPECT_EQ( return_status, 1);
  xyce_close( & xycePtr );
  EXPECT_TRUE( ((long *)xycePtr) == 0 );
}


TEST ( XyceCInterface, RunNetlist1)
{
  void * xycePtr = NULL;
  xyce_open( & xycePtr);
  EXPECT_TRUE( ((long *)xycePtr) != 0 );
  const char argv0[] = "Xyce";
  const char argv1[] = "TestNetlist1.cir"; 
  char * argvarray[2];
  argvarray[0] = (char *) &argv0[0];
  argvarray[1] = (char *) &argv1[0];
  int return_status = xyce_initialize( & xycePtr, 2, argvarray);
  EXPECT_EQ( return_status, 1);
  
  return_status = xyce_runSimulation( & xycePtr );
  EXPECT_EQ( return_status, 1);
  
  xyce_close( & xycePtr );
  EXPECT_TRUE( ((long *)xycePtr) == 0 );
}


TEST ( XyceCInterface, GetTimeNetlist1)
{
  void * xycePtr = NULL;
  xyce_open( & xycePtr);
  EXPECT_TRUE( ((long *)xycePtr) != 0 );
  const char argv0[] = "Xyce";
  const char argv1[] = "TestNetlist1.cir"; 
  char * argvarray[2];
  argvarray[0] = (char *) &argv0[0];
  argvarray[1] = (char *) &argv1[0];
  int return_status = xyce_initialize( & xycePtr, 2, argvarray);
  EXPECT_EQ( return_status, 1);
  
  return_status = xyce_runSimulation( & xycePtr );
  EXPECT_EQ( return_status, 1);
  
  double simTime = xyce_getTime( & xycePtr );
  EXPECT_EQ( simTime, 1.0 );
  
  xyce_close( & xycePtr );
  EXPECT_TRUE( ((long *)xycePtr) == 0 );
}


TEST ( XyceCInterface, GetFinalTimeNetlist1)
{
  void * xycePtr = NULL;
  xyce_open( & xycePtr);
  EXPECT_TRUE( ((long *)xycePtr) != 0 );
  const char argv0[] = "Xyce";
  const char argv1[] = "TestNetlist1.cir"; 
  char * argvarray[2];
  argvarray[0] = (char *) &argv0[0];
  argvarray[1] = (char *) &argv1[0];
  int return_status = xyce_initialize( & xycePtr, 2, argvarray);
  EXPECT_EQ( return_status, 1);
  
  double finalTime = xyce_getFinalTime( & xycePtr );
  EXPECT_EQ( finalTime, 1.0);
  
  xyce_close( & xycePtr );
  EXPECT_TRUE( ((long *)xycePtr) == 0 );
}



TEST ( XyceCInterface, MultiStepNetlist1)
{
  void * xycePtr = NULL;
  xyce_open( & xycePtr);
  EXPECT_TRUE( ((long *)xycePtr) != 0 );
  const char argv0[] = "Xyce";
  const char argv1[] = "TestNetlist1.cir"; 
  char * argvarray[2];
  argvarray[0] = (char *) &argv0[0];
  argvarray[1] = (char *) &argv1[0];
  int return_status = xyce_initialize( & xycePtr, 2, argvarray);
  EXPECT_EQ( return_status, 1);
  
  double finalTime = xyce_getFinalTime( & xycePtr );
  EXPECT_EQ( finalTime, 1.0);
  // run this simulation in several sub-stesps
  const int numSteps=100;
  for(auto i = 0; i<(numSteps+1); i++)
  {
    // simToTime must be greater than zero.  passing zero in will cause Xyce to abort.
    double simToTime = (i)*finalTime/numSteps;
    double actualTime=0.0;
    bool stepResult = xyce_simulateUntil( & xycePtr, simToTime, & actualTime);
    EXPECT_TRUE(stepResult);
    // for this simple circuit we expect the simToTime and actualTime to be equal 
    EXPECT_EQ( simToTime, actualTime );
    double reportedTime = xyce_getTime( & xycePtr );
    EXPECT_EQ( reportedTime, actualTime );
    bool isSimComplete = xyce_simulationComplete( & xycePtr );
    // should be false on all but the last step
    if( i==numSteps )
    {
      EXPECT_TRUE( isSimComplete );
    }
    else
    {
      EXPECT_FALSE( isSimComplete );
    }
  }
  xyce_close( & xycePtr );
  EXPECT_TRUE( ((long *)xycePtr) == 0 );
}

TEST ( XyceCInterface, MultiStepNetlist3)
{
  void * xycePtr = NULL;
  xyce_open( & xycePtr);
  EXPECT_TRUE( ((long *)xycePtr) != 0 );
  const char argv0[] = "Xyce";
  const char argv1[] = "TestNetlist3.cir"; 
  char * argvarray[2];
  argvarray[0] = (char *) &argv0[0];
  argvarray[1] = (char *) &argv1[0];
  int return_status = xyce_initialize( & xycePtr, 2, argvarray);
  EXPECT_EQ( return_status, 1);
  
  double finalTime = xyce_getFinalTime( & xycePtr );
  EXPECT_EQ( finalTime, 2.0e-3);
  // run this simulation in several sub-stesps
  const int numSteps=100;
  for(auto i = 0; i<(numSteps+1); i++)
  {
    // simToTime must be greater than zero.  passing zero in will cause Xyce to abort.
    double simToTime = (i)*finalTime/numSteps;
    double actualTime=0.0;
    bool stepResult = xyce_simulateUntil( & xycePtr, simToTime, & actualTime);
    EXPECT_TRUE(stepResult);
    // for this simple circuit we expect the simToTime and actualTime to be equal 
    EXPECT_EQ( simToTime, actualTime );
    double reportedTime = xyce_getTime( & xycePtr );
    EXPECT_EQ( reportedTime, actualTime );
    bool isSimComplete = xyce_simulationComplete( & xycePtr );
    // should be false on all but the last step
    if( i==numSteps )
    {
      EXPECT_TRUE( isSimComplete );
    }
    else
    {
      EXPECT_FALSE( isSimComplete );
    }
  }
  xyce_close( & xycePtr );
  EXPECT_TRUE( ((long *)xycePtr) == 0 );
}

TEST ( XyceCInterface, DACDeviceNamesTestNetlist1 )
{
  void * xycePtr = NULL;
  xyce_open( & xycePtr);
  EXPECT_TRUE( ((long *)xycePtr) != 0 );
  const char argv0[] = "Xyce";
  const char argv1[] = "TestNetlist1.cir"; 
  char * argvarray[2];
  argvarray[0] = (char *) &argv0[0];
  argvarray[1] = (char *) &argv1[0];
  int return_status = xyce_initialize( & xycePtr, 2, argvarray);
  EXPECT_EQ( return_status, 1);

  char * dacNames[10];
  int numDacs;
  bool circuitHasDACs = xyce_getDACDeviceNames( &xycePtr, & numDacs, dacNames);
  EXPECT_FALSE( circuitHasDACs);
  EXPECT_EQ( numDacs, 0);

  xyce_close( & xycePtr );
  EXPECT_TRUE( ((long *)xycePtr) == 0 );
}


TEST ( XyceCInterface, DACDeviceNamesTestNetlist3 )
{
  void * xycePtr = NULL;
  xyce_open( & xycePtr);
  EXPECT_TRUE( ((long *)xycePtr) != 0 );
  const char argv0[] = "Xyce";
  const char argv1[] = "TestNetlist3.cir"; 
  char * argvarray[2];
  argvarray[0] = (char *) &argv0[0];
  argvarray[1] = (char *) &argv1[0];
  int return_status = xyce_initialize( & xycePtr, 2, argvarray);
  EXPECT_EQ( return_status, 1);
  
  const int maxNumNames = 10;
  const int maxNameLength = 256;
  char * dacNames[maxNumNames];
  for( int i=0; i<maxNumNames; i++)
  {
    dacNames[i] = (char *) malloc(maxNameLength*sizeof(char));
  }
  int numDacs;
  bool circuitHasDACs = xyce_getDACDeviceNames( &xycePtr, & numDacs, dacNames);
  EXPECT_TRUE( circuitHasDACs);
  EXPECT_EQ( numDacs, 1);

  xyce_close( & xycePtr );
  EXPECT_TRUE( ((long *)xycePtr) == 0 );
  for( int i=0; i<maxNumNames; i++)
  {
    free(dacNames[i]);
  }
}


TEST ( XyceCInterface, CheckDeviceParamTestNetlist2 )
{
  void * xycePtr = NULL;
  xyce_open( & xycePtr);
  EXPECT_TRUE( ((long *)xycePtr) != 0 );
  const char argv0[] = "Xyce";
  const char argv1[] = "TestNetlist2.cir"; 
  char * argvarray[2];
  argvarray[0] = (char *) &argv0[0];
  argvarray[1] = (char *) &argv1[0];
  int return_status = xyce_initialize( & xycePtr, 2, argvarray);
  EXPECT_EQ( return_status, 1);
  
  // this should fail as TOMP is not a simulation parameter in this netlist 
  char param1[] = "TOMP";
  bool setParamResult = xyce_checkCircuitParameterExists( & xycePtr, param1);
  EXPECT_FALSE( setParamResult );
  
  // this should pass 
  char param2[] = "TEMP";
  setParamResult = xyce_checkCircuitParameterExists( & xycePtr, param2);
  EXPECT_TRUE( setParamResult );
  
  // this should passs
  char param3[] = "R1:R";
  setParamResult = xyce_checkCircuitParameterExists( & xycePtr, param3);
  EXPECT_TRUE( setParamResult );
  
  xyce_close( & xycePtr );
  EXPECT_TRUE( ((long *)xycePtr) == 0 );
}


TEST ( XyceCInterface, GetCircuitValuesTestNetlist2 )
{
  void * xycePtr = NULL;
  xyce_open( & xycePtr);
  EXPECT_TRUE( ((long *)xycePtr) != 0 );
  const char argv0[] = "Xyce";
  const char argv1[] = "TestNetlist2.cir"; 
  char * argvarray[2];
  argvarray[0] = (char *) &argv0[0];
  argvarray[1] = (char *) &argv1[0];
  int return_status = xyce_initialize( & xycePtr, 2, argvarray);
  EXPECT_EQ( return_status, 1);

  double circuitParamValue = 0.0;
  char param1[] = "R1:R";
  circuitParamValue = xyce_getCircuitValue( & xycePtr, param1 );
  EXPECT_EQ( circuitParamValue, 1.0 );
  char param2[] = "TEMP";
  circuitParamValue = xyce_getCircuitValue( & xycePtr, param2);
  EXPECT_EQ( circuitParamValue, 27.0 );
  char param3[] = "V(1)";
  circuitParamValue = xyce_getCircuitValue( & xycePtr, param3);
  EXPECT_NEAR( circuitParamValue, -1.0e-7,  1.0e-7 );
  
  
  double finalTime = xyce_getFinalTime( & xycePtr );
  EXPECT_EQ( finalTime, 1.0);
  
  // run this simulation in several sub-stesps
  const int numSteps=100;
  for(auto i = 0; i<numSteps; i++)
  {
    // simToTime must be greater than zero.  passing zero in will cause Xyce to abort.
    double simToTime = (i+1)*finalTime/numSteps;
    double actualTime=0.0;
    bool stepResult = xyce_simulateUntil( & xycePtr, simToTime, & actualTime);
    EXPECT_TRUE(stepResult);
    // for this simple circuit we expect the simToTime and actualTime to be equal 
    EXPECT_EQ( simToTime, actualTime );
    double reportedTime = xyce_getTime(& xycePtr);
    EXPECT_EQ( reportedTime, actualTime );
    if( i==(numSteps/2))
    {
      // in the middle of the run change the resistance on resistor R1
      bool setParamResult = xyce_setCircuitParameter( & xycePtr, param1, 1000.0);
      EXPECT_TRUE( setParamResult );
      setParamResult = xyce_setCircuitParameter( & xycePtr, param2, -50.0);
      EXPECT_TRUE( setParamResult );
    }
    if( i==((numSteps/2)+1))
    {
      // in the middle of the run change the resistance on resistor R1
      circuitParamValue = 0.0;
      circuitParamValue = xyce_getCircuitValue( & xycePtr, param1);
      EXPECT_EQ( circuitParamValue, 1000.0 );
      circuitParamValue = xyce_getCircuitValue( & xycePtr, param2);
      EXPECT_EQ( circuitParamValue, -50.0 );
    }
    // set this to value out of the range of V(1) during this simulation 
    // V(1) goes from -1.0 to +1.0
    circuitParamValue = -1.001;
    circuitParamValue = xyce_getCircuitValue( & xycePtr, param3);
    EXPECT_TRUE( (circuitParamValue >= -1.0) && (circuitParamValue <= 1.0));
    
    bool isSimComplete = xyce_simulationComplete(& xycePtr );
    // should be false on all but the last step
    if( i==(numSteps-1))
    {
      EXPECT_TRUE( isSimComplete );
    }
    else
    {
      EXPECT_FALSE( isSimComplete );
    }
  }
  
  // try to get values that were in measure functions.
  char param4[] = "MAXV1";
  circuitParamValue = xyce_getCircuitValue( & xycePtr, param4);
  EXPECT_NEAR( circuitParamValue, 1.0, 1.0e-7 );
  char param5[] = "MINV1";
  circuitParamValue = xyce_getCircuitValue( & xycePtr, param5);
  EXPECT_NEAR( circuitParamValue, -1.0,  1.0e-7 );
  
  xyce_close( & xycePtr );
  EXPECT_TRUE( ((long *)xycePtr) == 0 );
}



TEST ( XyceCInterface, GetADCMapTestNetlist1 )
{
  void * xycePtr = NULL;
  xyce_open( & xycePtr);
  EXPECT_TRUE( ((long *)xycePtr) != 0 );
  const char argv0[] = "Xyce";
  const char argv1[] = "TestNetlist1.cir"; 
  char * argvarray[2];
  argvarray[0] = (char *) &argv0[0];
  argvarray[1] = (char *) &argv1[0];
  int return_status = xyce_initialize( & xycePtr, 2, argvarray);
  EXPECT_EQ( return_status, 1);

  const int maxNumADC = 10;
  const int maxNameLength = 256;
  char * adcNames[maxNumADC];
  for( int i=0; i<maxNumADC; i++)
  {
    adcNames[i] = (char *) malloc(maxNameLength*sizeof(char));
  }
  int numAdcs;
  int adcWidths[maxNumADC];
  double adcRes[maxNumADC];
  double adcUpperVLim[maxNumADC];
  double adcLowerVLim[maxNumADC];
  double adcSettlingTime[maxNumADC];
  bool circuitHasADCs = xyce_getADCMap( &xycePtr, & numAdcs, adcNames, 
    adcWidths, adcRes, adcUpperVLim, adcLowerVLim, adcSettlingTime );
	           
  EXPECT_FALSE( circuitHasADCs);
  EXPECT_EQ( numAdcs, 0);

  xyce_close( & xycePtr );
  EXPECT_TRUE( ((long *)xycePtr) == 0 );
  for( int i=0; i<maxNumADC; i++)
  {
    free( adcNames[i]);
  }
}

TEST ( XyceCInterface, GetADCMapTestNetlist3 )
{
  void * xycePtr = NULL;
  xyce_open( & xycePtr);
  EXPECT_TRUE( ((long *)xycePtr) != 0 );
  const char argv0[] = "Xyce";
  const char argv1[] = "TestNetlist3.cir"; 
  char * argvarray[2];
  argvarray[0] = (char *) &argv0[0];
  argvarray[1] = (char *) &argv1[0];
  int return_status = xyce_initialize( & xycePtr, 2, argvarray);
  EXPECT_EQ( return_status, 1);

  const int maxNumADC = 10;
  const int maxNameLength = 256;
  char * adcNames[maxNumADC];
  for( int i=0; i<maxNumADC; i++)
  {
    adcNames[i] = (char *) malloc(maxNameLength*sizeof(char));
  }
  int numAdcs;
  int adcWidths[maxNumADC];
  double adcRes[maxNumADC];
  double adcUpperVLim[maxNumADC];
  double adcLowerVLim[maxNumADC];
  double adcSettlingTime[maxNumADC];
  bool circuitHasADCs = xyce_getADCMap( &xycePtr, & numAdcs, adcNames, 
    adcWidths, adcRes, adcUpperVLim, adcLowerVLim, adcSettlingTime );
	           
  EXPECT_TRUE( circuitHasADCs);
  EXPECT_EQ( numAdcs, 2);

  xyce_close( & xycePtr );
  EXPECT_TRUE( ((long *)xycePtr) == 0 );
  for( int i=0; i<maxNumADC; i++)
  {
    free( adcNames[i]);
  }
}


TEST ( XyceCInterface, GetTimeVoltagePairsNetlist3 )
{
  void * xycePtr = NULL;
  xyce_open( & xycePtr);
  EXPECT_TRUE( ((long *)xycePtr) != 0 );
  const char argv0[] = "Xyce";
  const char argv1[] = "TestNetlist3.cir"; 
  char * argvarray[2];
  argvarray[0] = (char *) &argv0[0];
  argvarray[1] = (char *) &argv1[0];
  int return_status = xyce_initialize( & xycePtr, 2, argvarray);
  EXPECT_EQ( return_status, 1);

  const int maxNumADC = 10;
  const int maxNameLength = 256;
  char * adcNames[maxNumADC];
  for( int i=0; i<maxNumADC; i++)
  {
    adcNames[i] = (char *) malloc(maxNameLength*sizeof(char));
  }
  int numAdcs;
  int adcWidths[maxNumADC];
  double adcRes[maxNumADC];
  double adcUpperVLim[maxNumADC];
  double adcLowerVLim[maxNumADC];
  double adcSettlingTime[maxNumADC];
  bool circuitHasADCs = xyce_getADCMap( &xycePtr, & numAdcs, adcNames, 
    adcWidths, adcRes, adcUpperVLim, adcLowerVLim, adcSettlingTime );
	           
  EXPECT_TRUE( circuitHasADCs);
  EXPECT_EQ( numAdcs, 2);
  
  // simulate forward in time about half way through the simulation
  double finalTime = xyce_getFinalTime( & xycePtr );
  double simToTime = 0.5*finalTime;
  double actualTime=0.0;
  bool stepResult = xyce_simulateUntil( & xycePtr, simToTime, & actualTime);
  EXPECT_TRUE(stepResult);
  // for this simple circuit we expect the simToTime and actualTime to be equal 
  EXPECT_EQ( simToTime, actualTime );
  double reportedTime = xyce_getTime( & xycePtr );
  EXPECT_EQ( reportedTime, actualTime );
  bool isSimComplete = xyce_simulationComplete( & xycePtr );
  // should be false on all but the last step
  EXPECT_FALSE( isSimComplete );
  
  int numADCret = 0;
  int numPoints = 0;
  const int maxDataPoints=10000;
  double * timeArray[maxNumADC];
  double * vArray[maxNumADC];
  for( int i=0; i<maxNumADC; i++)
  {
    timeArray[i] = (double *) malloc(maxDataPoints*sizeof(double));
    vArray[i] = (double *) malloc(maxDataPoints*sizeof(double));
  }
  int retResult = xyce_getTimeVoltagePairsADC( & xycePtr, & numADCret, adcNames, & numPoints, timeArray, vArray );
  
  EXPECT_EQ( retResult, 1);
  EXPECT_EQ(numADCret, 2);
  EXPECT_EQ(numPoints, 2);

  xyce_close( & xycePtr );
  EXPECT_TRUE( ((long *)xycePtr) == 0 );
  for( int i=0; i<maxNumADC; i++)
  {
    free( adcNames[i]);
    free( timeArray[i]);
    free( vArray[i]);
  }
}

TEST ( XyceCInterface, UpdateTimeVoltagePairsNetlist3 )
{
  void * xycePtr = NULL;
  xyce_open( & xycePtr);
  EXPECT_TRUE( ((long *)xycePtr) != 0 );
  const char argv0[] = "Xyce";
  const char argv1[] = "TestNetlist3.cir"; 
  char * argvarray[2];
  argvarray[0] = (char *) &argv0[0];
  argvarray[1] = (char *) &argv1[0];
  int return_status = xyce_initialize( & xycePtr, 2, argvarray);
  EXPECT_EQ( return_status, 1);
  
  const int maxNumDev = 10;
  const int maxNameLength = 256;
  char * dacNames[maxNumDev];
  for( int i=0; i<maxNumDev; i++)
  {
    dacNames[i] = (char *) malloc(maxNameLength*sizeof(char));
  }
  int numDacs;
  bool circuitHasDACs = xyce_getDACDeviceNames( &xycePtr, & numDacs, dacNames);
  EXPECT_TRUE( circuitHasDACs);
  EXPECT_EQ( numDacs, 1);
  
  char * adcNames[maxNumDev];
  for( int i=0; i<maxNumDev; i++)
  {
    adcNames[i] = (char *) malloc(maxNameLength*sizeof(char));
  }
  int numAdcs;
  int adcWidths[maxNumDev];
  double adcRes[maxNumDev];
  double adcUpperVLim[maxNumDev];
  double adcLowerVLim[maxNumDev];
  double adcSettlingTime[maxNumDev];
  bool circuitHasADCs = xyce_getADCMap( &xycePtr, & numAdcs, adcNames, 
    adcWidths, adcRes, adcUpperVLim, adcLowerVLim, adcSettlingTime );
	           
  EXPECT_TRUE( circuitHasADCs);
  EXPECT_EQ( numAdcs, 2);
  
  // simulate forward in time about half way through the simulation
  double finalTime = xyce_getFinalTime( & xycePtr );
  double simToTime = 0.5*finalTime;
  double actualTime=0.0;
  bool stepResult = xyce_simulateUntil( & xycePtr, simToTime, & actualTime);
  EXPECT_TRUE(stepResult);
  // for this simple circuit we expect the simToTime and actualTime to be equal 
  EXPECT_EQ( simToTime, actualTime );
  double reportedTime = xyce_getTime( & xycePtr );
  EXPECT_EQ( reportedTime, actualTime );
  bool isSimComplete = xyce_simulationComplete( & xycePtr );
  // should be false on all but the last step
  EXPECT_FALSE( isSimComplete );
  
  int numADCret = 0;
  int numPoints = 0;
  const int maxDataPoints=100;
  double * timeArray[maxNumDev];
  double * vArray[maxNumDev];
  for( int i=0; i<maxNumDev; i++)
  {
    timeArray[i] = (double *) malloc(maxDataPoints*sizeof(double));
    vArray[i] = (double *) malloc(maxDataPoints*sizeof(double));
  }
  int retResult = xyce_getTimeVoltagePairsADC( & xycePtr, & numADCret, adcNames, & numPoints, timeArray, vArray );
  
  EXPECT_EQ( retResult, 1);
  EXPECT_EQ(numADCret, 2);
  EXPECT_EQ(numPoints, 2);
  
  const int numTVPairs = 2;
  double timeData[numTVPairs];
  double dacVData[numTVPairs];
  timeData[0] = actualTime + 1e-7;
  timeData[1] = actualTime + 1e-3;
  dacVData[0] = 0.5;
  dacVData[1] = 0.5;
  // set up data for the DAC
  retResult = xyce_updateTimeVoltagePairs(& xycePtr, dacNames[0], numTVPairs, timeData, dacVData);
  EXPECT_EQ( retResult, 1);
  
  return_status = xyce_runSimulation( & xycePtr );
  EXPECT_EQ( return_status, 1);
  
  char param3[] = "V(N1)";
  double circuitParamValue = xyce_getCircuitValue( & xycePtr, param3);
  EXPECT_NEAR( circuitParamValue, 0.5,  1.0e-7 );

  xyce_close( & xycePtr );
  EXPECT_TRUE( ((long *)xycePtr) == 0 );
  for( int i=0; i<maxNumDev; i++)
  {
    free( adcNames[i]);
    free( timeArray[i]);
    free( vArray[i]);
    free(dacNames[i]);
  }
}


TEST ( XyceCInterface, MultiPortGetTimeVoltagePairs )
{
  void * xycePtr = NULL;
  xyce_open( & xycePtr);
  EXPECT_TRUE( ((long *)xycePtr) != 0 );
  const char argv0[] = "Xyce";
  const char argv1[] = "MultiPort.cir"; 
  char * argvarray[2];
  argvarray[0] = (char *) &argv0[0];
  argvarray[1] = (char *) &argv1[0];
  int return_status = xyce_initialize( & xycePtr, 2, argvarray);
  EXPECT_EQ( return_status, 1);

  const int maxNumADC = 10;
  const int maxNameLength = 256;
  char * adcNames[maxNumADC];
  for( int i=0; i<maxNumADC; i++)
  {
    adcNames[i] = (char *) malloc(maxNameLength*sizeof(char));
  }
  int numAdcs;
  int adcWidths[maxNumADC];
  double adcRes[maxNumADC];
  double adcUpperVLim[maxNumADC];
  double adcLowerVLim[maxNumADC];
  double adcSettlingTime[maxNumADC];
  bool circuitHasADCs = xyce_getADCMap( &xycePtr, & numAdcs, adcNames, 
    adcWidths, adcRes, adcUpperVLim, adcLowerVLim, adcSettlingTime );
	           
  EXPECT_TRUE( circuitHasADCs);
  EXPECT_EQ( numAdcs, 2);
  
  // simulate forward in time about half way through the simulation
  double finalTime = xyce_getFinalTime( & xycePtr );
  // run this simulation in several sub-stesps
  const int numSteps=100;
  double * timeArray[maxNumADC];
  double * vArray[maxNumADC];
  const int maxDataPoints=10000;  
  for( int i=0; i<maxNumADC; i++)
  {
    timeArray[i] = (double *) malloc(maxDataPoints*sizeof(double));
    vArray[i] = (double *) malloc(maxDataPoints*sizeof(double));
  }
  
  for(auto i = 0; i<numSteps; i++)
  {
    // simToTime must be greater than zero.  passing zero in will cause Xyce to abort.
    double simToTime = (i+1)*finalTime/numSteps;
    double actualTime=0.0;
    bool stepResult = xyce_simulateUntil( & xycePtr, simToTime, & actualTime);
    EXPECT_TRUE(stepResult);
    // for this simple circuit we expect the simToTime and actualTime to be equal 
    EXPECT_EQ( simToTime, actualTime );
    double reportedTime = xyce_getTime( & xycePtr );
    EXPECT_EQ( reportedTime, actualTime );
    bool isSimComplete = xyce_simulationComplete( & xycePtr );
    // should be false on all but the last step
    if( i < (numSteps-1))
    {
      EXPECT_FALSE( isSimComplete );
    }
    else
    {
      EXPECT_TRUE( isSimComplete );
    }
    
    int numADCret = 0;
    int numPoints = 0;
    int retResult = xyce_getTimeVoltagePairsADC( & xycePtr, & numADCret, adcNames, & numPoints, timeArray, vArray );
    
    EXPECT_EQ( retResult, 1);
    EXPECT_EQ(numADCret, 2);
    for( auto adcID=0; adcID<numADCret; adcID++)
    {
      if( strncmp("YADC!ADCHIGH", adcNames[adcID], 12) == 0 )
      {
        for( auto ts=0; ts<numPoints; ts++)
        {
          if( fabs(timeArray[adcID][ts] - reportedTime) < 1e-8 )
          {
            EXPECT_NEAR(vArray[adcID][ts], 5.0, 1e-7 );
            break;
          }
        }
      }
      if( strncmp("YADC!ADCLOW", adcNames[adcID], 12) == 0 )
      {
        for( auto ts=0; ts<numPoints; ts++)
        {
          if( fabs(timeArray[adcID][ts] - reportedTime) < 1e-8 )
          {
            EXPECT_NEAR(vArray[adcID][ts], -3.0, 1e-7 );
            break;
          }
        }
      }
    }
  }
  
  xyce_close( & xycePtr );
  EXPECT_TRUE( ((long *)xycePtr) == 0 );
  for( int i=0; i<maxNumADC; i++)
  {
    free( adcNames[i]);
    free( timeArray[i]);
    free( vArray[i]);
  }
}

TEST ( XyceCInterface, MultiPortGetTimeVoltagePairsTimeVariant )
{
  void * xycePtr = NULL;
  xyce_open( & xycePtr);
  EXPECT_TRUE( ((long *)xycePtr) != 0 );
  const char argv0[] = "Xyce";
  const char argv1[] = "MultiPort2.cir"; 
  char * argvarray[2];
  argvarray[0] = (char *) &argv0[0];
  argvarray[1] = (char *) &argv1[0];
  int return_status = xyce_initialize( & xycePtr, 2, argvarray);
  EXPECT_EQ( return_status, 1);

  const int maxNumADC = 10;
  const int maxNameLength = 256;
  char * adcNames[maxNumADC];
  for( int i=0; i<maxNumADC; i++)
  {
    adcNames[i] = (char *) malloc(maxNameLength*sizeof(char));
  }
  int numAdcs;
  int adcWidths[maxNumADC];
  double adcRes[maxNumADC];
  double adcUpperVLim[maxNumADC];
  double adcLowerVLim[maxNumADC];
  double adcSettlingTime[maxNumADC];
  bool circuitHasADCs = xyce_getADCMap( &xycePtr, & numAdcs, adcNames, 
    adcWidths, adcRes, adcUpperVLim, adcLowerVLim, adcSettlingTime );
	           
  EXPECT_TRUE( circuitHasADCs);
  EXPECT_EQ( numAdcs, 2);
  
  // simulate forward in time about half way through the simulation
  double finalTime = xyce_getFinalTime( & xycePtr );
  // run this simulation in several sub-stesps
  const int numSteps=100;
  double * timeArray[maxNumADC];
  double * vArray[maxNumADC];
  const int maxDataPoints=10000;  
  for( int i=0; i<maxNumADC; i++)
  {
    timeArray[i] = (double *) malloc(maxDataPoints*sizeof(double));
    vArray[i] = (double *) malloc(maxDataPoints*sizeof(double));
  }
  
  for(auto i = 0; i<numSteps; i++)
  {
    // simToTime must be greater than zero.  passing zero in will cause Xyce to abort.
    double simToTime = (i+1)*finalTime/numSteps;
    double actualTime=0.0;
    bool stepResult = xyce_simulateUntil( & xycePtr, simToTime, & actualTime);
    EXPECT_TRUE(stepResult);
    // for this simple circuit we expect the simToTime and actualTime to be equal 
    EXPECT_EQ( simToTime, actualTime );
    double reportedTime = xyce_getTime( & xycePtr );
    EXPECT_EQ( reportedTime, actualTime );
    bool isSimComplete = xyce_simulationComplete( & xycePtr );
    // should be false on all but the last step
    if( i < (numSteps-1))
    {
      EXPECT_FALSE( isSimComplete );
    }
    else
    {
      EXPECT_TRUE( isSimComplete );
    }
    
    int numADCret = 0;
    int numPoints = 0;
    int retResult = xyce_getTimeVoltagePairsADC( & xycePtr, & numADCret, adcNames, & numPoints, timeArray, vArray );
    
    EXPECT_EQ( retResult, 1);
    EXPECT_EQ(numADCret, 2);
    for( auto adcID=0; adcID<numADCret; adcID++)
    {
      if( strncmp("YADC!ADCHIGH", adcNames[adcID], 12) == 0 )
      {
        for( auto ts=0; ts<numPoints; ts++)
        {
          if( fabs(timeArray[adcID][ts] - reportedTime) < 1e-8 )
          {
            double ans = 5.0*sin(2*M_PI*5e3*simToTime);
            //printf("%s: %g: %g %g diff %g\n", adcNames[adcID], simToTime, vArray[adcID][ts], ans, (vArray[adcID][ts]-ans));
            EXPECT_NEAR(vArray[adcID][ts], ans, 1e-7 );
            break;
          }
        }
      }
      if( strncmp("YADC!ADCLOW", adcNames[adcID], 12) == 0 )
      {
        for( auto ts=0; ts<numPoints; ts++)
        {
          if( fabs(timeArray[adcID][ts] - reportedTime) < 1e-8 )
          {
            double ans = 3.0*cos(2*M_PI*2e3*simToTime);
            //printf("%s: %g: %g %g diff %g\n", adcNames[adcID], simToTime, vArray[adcID][ts], ans, (vArray[adcID][ts]-ans));
            EXPECT_NEAR(vArray[adcID][ts], ans, 1e-7 );
            break;
          }
        }
      }
    }
  }
  
  xyce_close( & xycePtr );
  EXPECT_TRUE( ((long *)xycePtr) == 0 );
  for( int i=0; i<maxNumADC; i++)
  {
    free( adcNames[i]);
    free( timeArray[i]);
    free( vArray[i]);
  }
}


TEST ( XyceCInterface, SetSimulationDirectory )
{
  // This is a test to show that Xyce API function setWorkingDirectory() does change the working directory
#if (__cplusplus>=201703L)
  void * xycePtr = NULL;
  xyce_open( & xycePtr);
  std::filesystem::path startingDirectory = std::filesystem::current_path();
  std::filesystem::create_directory("TestDir");
  const char * dirName = "TestDir";
  xyce_set_working_directory( &xycePtr, dirName);
  std::filesystem::path currentDirectory = std::filesystem::current_path();
  EXPECT_NE( startingDirectory.string(), currentDirectory.string());
  xyce_close( & xycePtr );
#endif
}

TEST ( XyceCInterface, InitFromWorkingDirectory )
{
  // this test ensures that using Xyce API function setWorkingDirectory() allows 
  // Xyce it initialize itself with a netlist in that subdirectory.
#if (__cplusplus>=201703L)
  void * xycePtr = NULL;
  xyce_open( & xycePtr);
  
  std::filesystem::path startingDirectory = std::filesystem::current_path();
  std::filesystem::path workingDir("TestDir2");
  std::filesystem::create_directory(workingDir);
  std::string netlist("TestNetlist1.cir");
  std::filesystem::path originalFile(startingDirectory);
  originalFile /= netlist;
  std::filesystem::path newFile(startingDirectory);
  newFile /= workingDir;
  newFile /=netlist;
  std::filesystem::copy_file( originalFile, newFile, std::filesystem::copy_options::overwrite_existing);
  const char * dirName = "TestDir2";
  xyce_set_working_directory( &xycePtr, dirName);
  
  const int numArgs = 4;
  const char argv0[] = "Xyce";
  const char argv1[] = "TestNetlist1.cir"; 
  const char argv2[] = "-l";
  const char argv3[] = "TestNetlist2.cir.log";
  char * argvarray[numArgs];
  argvarray[0] = (char *) &argv0[0];
  argvarray[1] = (char *) &argv1[0];
  argvarray[2] = (char *) &argv2[0];
  argvarray[3] = (char *) &argv3[0];
  int return_status = xyce_initialize( & xycePtr, 4, argvarray);
  EXPECT_EQ( return_status, 1);
  xyce_close( & xycePtr );
#endif
}

//-------------------------------------------------------------------------------
int main (int argc, char **argv)
{
  testing::InitGoogleTest(&argc, argv);

  return RUN_ALL_TESTS();
}

