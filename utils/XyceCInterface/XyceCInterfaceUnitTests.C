//-------------------------------------------------------------------------
//   Copyright 2002-2022 National Technology & Engineering Solutions of
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

#include <gtest/gtest.h>
#include "Xyce_config.h"
#include <N_CIR_XyceCInterface.h>

//
// Xyce::Circuit::Simulator functions that need to be tested
//
// x constructor
// x initialize()
// x finalize()
// x getTime()
// x getDACDeviceNames()
// x getADCMap()
// getTimeVoltagePairs()
// updateTimeVoltagePairs()
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

TEST ( XyceSimulator, DACDeviceNamesTestNetlist1 )
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


TEST ( XyceSimulator, DACDeviceNamesTestNetlist3 )
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
  
  char * dacNames[10];
  int numDacs;
  bool circuitHasDACs = xyce_getDACDeviceNames( &xycePtr, & numDacs, dacNames);
  EXPECT_TRUE( circuitHasDACs);
  EXPECT_EQ( numDacs, 1);

  xyce_close( & xycePtr );
  EXPECT_TRUE( ((long *)xycePtr) == 0 );
}


TEST ( XyceSimulator, CheckDeviceParamTestNetlist2 )
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


TEST ( XyceSimulator, GetCircuitValuesTestNetlist2 )
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


//-------------------------------------------------------------------------------
int main (int argc, char **argv)
{
  testing::InitGoogleTest(&argc, argv);

  return RUN_ALL_TESTS();
}

