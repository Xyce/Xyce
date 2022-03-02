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
#include <N_CIR_Xyce.h>

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
// setCircuitParameter()
// getCircuitValue()
// x simulateUntil() in stead of provisionalStep() & acceptProvisionalStep()
// 

TEST ( XyceSimulator, create)
{
  Xyce::Circuit::Simulator * xycePtr = NULL;
  xycePtr = new Xyce::Circuit::Simulator();
  EXPECT_TRUE( xycePtr != NULL );
  delete xycePtr;
}

TEST ( XyceSimulator, InitializeTestNetlist1 )
{
  Xyce::Circuit::Simulator * xycePtr = NULL;
  xycePtr = new Xyce::Circuit::Simulator();
  EXPECT_TRUE( xycePtr != NULL );
  int numArgs = 4;
  char * cmdLineArgs[numArgs];
  std::string xyceBin = "XyceTests";
  std::string netlist = "TestNetlist1.cir";
  std::string logOption = "-l";
  std::string logFile = netlist + ".log";
  cmdLineArgs[0] = const_cast<char *>(xyceBin.c_str());
  cmdLineArgs[1] = const_cast<char *>(logOption.c_str());
  cmdLineArgs[2] = const_cast<char *>(logFile.c_str());
  cmdLineArgs[3] = const_cast<char *>(netlist.c_str());
  
  Xyce::Circuit::Simulator::RunStatus status = xycePtr->initialize( numArgs, cmdLineArgs);
  EXPECT_EQ( status, Xyce::Circuit::Simulator::RunStatus::SUCCESS );
  delete xycePtr;
}

TEST ( XyceSimulator, InitializeTestNetlist2 )
{
  Xyce::Circuit::Simulator * xycePtr = NULL;
  xycePtr = new Xyce::Circuit::Simulator();
  EXPECT_TRUE( xycePtr != NULL );
  int numArgs = 2;
  char * cmdLineArgs[numArgs];
  std::string xyceBin = "XyceTests";
  std::string netlist = "TestNetlist2.cir";
  cmdLineArgs[0] = const_cast<char *>(xyceBin.c_str());
  cmdLineArgs[1] = const_cast<char *>(netlist.c_str());
  Xyce::Circuit::Simulator::RunStatus status = xycePtr->initialize( numArgs, cmdLineArgs);
  EXPECT_EQ( status, Xyce::Circuit::Simulator::RunStatus::SUCCESS );
  delete xycePtr;
}

TEST ( XyceSimulator, InitializeTestNetlist3 )
{
  Xyce::Circuit::Simulator * xycePtr = NULL;
  xycePtr = new Xyce::Circuit::Simulator();
  EXPECT_TRUE( xycePtr != NULL );
  int numArgs = 2;
  char * cmdLineArgs[numArgs];
  std::string xyceBin = "XyceTests";
  std::string netlist = "TestNetlist3.cir";
  cmdLineArgs[0] = const_cast<char *>(xyceBin.c_str());
  cmdLineArgs[1] = const_cast<char *>(netlist.c_str());
  Xyce::Circuit::Simulator::RunStatus status = xycePtr->initialize( numArgs, cmdLineArgs);
  EXPECT_EQ( status, Xyce::Circuit::Simulator::RunStatus::SUCCESS );
  delete xycePtr;
}

TEST ( XyceSimulator, RunTestNetlist1 )
{
  Xyce::Circuit::Simulator * xycePtr = NULL;
  xycePtr = new Xyce::Circuit::Simulator();
  EXPECT_TRUE( xycePtr != NULL );
  int numArgs = 2;
  char * cmdLineArgs[numArgs];
  std::string xyceBin = "XyceTests";
  std::string netlist = "TestNetlist1.cir";
  cmdLineArgs[0] = const_cast<char *>(xyceBin.c_str());
  cmdLineArgs[1] = const_cast<char *>(netlist.c_str());
  Xyce::Circuit::Simulator::RunStatus status = xycePtr->initialize( numArgs, cmdLineArgs);
  EXPECT_EQ( status, Xyce::Circuit::Simulator::RunStatus::SUCCESS );
  status = xycePtr->runSimulation();
  EXPECT_EQ( status, Xyce::Circuit::Simulator::RunStatus::SUCCESS );
  delete xycePtr;
}


TEST ( XyceSimulator, FinalizeTestNetlist1 )
{
  Xyce::Circuit::Simulator * xycePtr = NULL;
  xycePtr = new Xyce::Circuit::Simulator();
  EXPECT_TRUE( xycePtr != NULL );
  int numArgs = 2;
  char * cmdLineArgs[numArgs];
  std::string xyceBin = "XyceTests";
  std::string netlist = "TestNetlist1.cir";
  cmdLineArgs[0] = const_cast<char *>(xyceBin.c_str());
  cmdLineArgs[1] = const_cast<char *>(netlist.c_str());
  Xyce::Circuit::Simulator::RunStatus status = xycePtr->initialize( numArgs, cmdLineArgs);
  EXPECT_EQ( status, Xyce::Circuit::Simulator::RunStatus::SUCCESS );
  status = xycePtr->runSimulation();
  EXPECT_EQ( status, Xyce::Circuit::Simulator::RunStatus::SUCCESS );
  status = xycePtr->finalize();
  EXPECT_EQ( status, Xyce::Circuit::Simulator::RunStatus::SUCCESS );
  delete xycePtr;
}

TEST ( XyceSimulator, GetTimeTestNetlist1 )
{
  Xyce::Circuit::Simulator * xycePtr = NULL;
  xycePtr = new Xyce::Circuit::Simulator();
  EXPECT_TRUE( xycePtr != NULL );
  int numArgs = 2;
  char * cmdLineArgs[numArgs];
  std::string xyceBin = "XyceTests";
  std::string netlist = "TestNetlist1.cir";
  cmdLineArgs[0] = const_cast<char *>(xyceBin.c_str());
  cmdLineArgs[1] = const_cast<char *>(netlist.c_str());
  Xyce::Circuit::Simulator::RunStatus status = xycePtr->initialize( numArgs, cmdLineArgs);
  EXPECT_EQ( status, Xyce::Circuit::Simulator::RunStatus::SUCCESS );
  status = xycePtr->runSimulation();
  EXPECT_EQ( status, Xyce::Circuit::Simulator::RunStatus::SUCCESS );
  double simTime = xycePtr->getTime();
  EXPECT_EQ( simTime, 1.0 );
  status = xycePtr->finalize();
  EXPECT_EQ( status, Xyce::Circuit::Simulator::RunStatus::SUCCESS );
  delete xycePtr;
}


TEST ( XyceSimulator, MultiTimeStepTestNetlist1 )
{
  Xyce::Circuit::Simulator * xycePtr = NULL;
  xycePtr = new Xyce::Circuit::Simulator();
  EXPECT_TRUE( xycePtr != NULL );
  int numArgs = 2;
  char * cmdLineArgs[numArgs];
  std::string xyceBin = "XyceTests";
  std::string netlist = "TestNetlist1.cir";
  cmdLineArgs[0] = const_cast<char *>(xyceBin.c_str());
  cmdLineArgs[1] = const_cast<char *>(netlist.c_str());
  Xyce::Circuit::Simulator::RunStatus status = xycePtr->initialize( numArgs, cmdLineArgs);
  EXPECT_EQ( status, Xyce::Circuit::Simulator::RunStatus::SUCCESS );
  double finalSimTime = xycePtr->getFinalTime();
  EXPECT_EQ( finalSimTime, 1.0 );
  // run this simulation in several sub-stesps
  const int numSteps=100;
  for(auto i = 0; i<numSteps; i++)
  {
    // simToTime must be greater than zero.  passing zero in will cause Xyce to abort.
    double simToTime = (i+1)*finalSimTime/numSteps;
    double actualTime=0.0;
    bool stepResult = xycePtr->simulateUntil(simToTime, actualTime);
    EXPECT_TRUE(stepResult);
    // for this simple circuit we expect the simToTime and actualTime to be equal 
    EXPECT_EQ( simToTime, actualTime );
    double reportedTime = xycePtr->getTime();
    EXPECT_EQ( reportedTime, actualTime );
    bool isSimComplete = xycePtr->simulationComplete();
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
  status = xycePtr->finalize();
  EXPECT_EQ( status, Xyce::Circuit::Simulator::RunStatus::SUCCESS );
  delete xycePtr;
}

TEST ( XyceSimulator, MultiTimeStepTestNetlist2 )
{
  Xyce::Circuit::Simulator * xycePtr = NULL;
  xycePtr = new Xyce::Circuit::Simulator();
  EXPECT_TRUE( xycePtr != NULL );
  int numArgs = 2;
  char * cmdLineArgs[numArgs];
  std::string xyceBin = "XyceTests";
  std::string netlist = "TestNetlist2.cir";
  cmdLineArgs[0] = const_cast<char *>(xyceBin.c_str());
  cmdLineArgs[1] = const_cast<char *>(netlist.c_str());
  Xyce::Circuit::Simulator::RunStatus status = xycePtr->initialize( numArgs, cmdLineArgs);
  EXPECT_EQ( status, Xyce::Circuit::Simulator::RunStatus::SUCCESS );
  double finalSimTime = xycePtr->getFinalTime();
  EXPECT_EQ( finalSimTime, 1.0 );
  // run this simulation in several sub-stesps
  const int numSteps=100;
  for(auto i = 0; i<numSteps; i++)
  {
    // simToTime must be greater than zero.  passing zero in will cause Xyce to abort.
    double simToTime = (i+1)*finalSimTime/numSteps;
    double actualTime=0.0;
    bool stepResult = xycePtr->simulateUntil(simToTime, actualTime);
    EXPECT_TRUE(stepResult);
    // for this simple circuit we expect the simToTime and actualTime to be equal 
    EXPECT_EQ( simToTime, actualTime );
    double reportedTime = xycePtr->getTime();
    EXPECT_EQ( reportedTime, actualTime );
    bool isSimComplete = xycePtr->simulationComplete();
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
  status = xycePtr->finalize();
  EXPECT_EQ( status, Xyce::Circuit::Simulator::RunStatus::SUCCESS );
  delete xycePtr;
}

TEST ( XyceSimulator, MultiTimeStepTestNetlist3 )
{
  Xyce::Circuit::Simulator * xycePtr = NULL;
  xycePtr = new Xyce::Circuit::Simulator();
  EXPECT_TRUE( xycePtr != NULL );
  int numArgs = 2;
  char * cmdLineArgs[numArgs];
  std::string xyceBin = "XyceTests";
  std::string netlist = "TestNetlist3.cir";
  cmdLineArgs[0] = const_cast<char *>(xyceBin.c_str());
  cmdLineArgs[1] = const_cast<char *>(netlist.c_str());
  Xyce::Circuit::Simulator::RunStatus status = xycePtr->initialize( numArgs, cmdLineArgs);
  EXPECT_EQ( status, Xyce::Circuit::Simulator::RunStatus::SUCCESS );
  double finalSimTime = xycePtr->getFinalTime();
  EXPECT_EQ( finalSimTime, 20.0e-4 );
  // run this simulation in several sub-stesps
  const int numSteps=100;
  for(auto i = 0; i<numSteps; i++)
  {
    // simToTime must be greater than zero.  passing zero in will cause Xyce to abort.
    double simToTime = (i+1)*finalSimTime/numSteps;
    double actualTime=0.0;
    bool stepResult = xycePtr->simulateUntil(simToTime, actualTime);
    EXPECT_TRUE(stepResult);
    // for this simple circuit we expect the simToTime and actualTime to be equal 
    EXPECT_EQ( simToTime, actualTime );
    double reportedTime = xycePtr->getTime();
    EXPECT_EQ( reportedTime, actualTime );
    bool isSimComplete = xycePtr->simulationComplete();
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
  status = xycePtr->finalize();
  EXPECT_EQ( status, Xyce::Circuit::Simulator::RunStatus::SUCCESS );
  delete xycePtr;
}


TEST ( XyceSimulator, DACDeviceNamesTestNetlist1 )
{
  Xyce::Circuit::Simulator * xycePtr = NULL;
  xycePtr = new Xyce::Circuit::Simulator();
  EXPECT_TRUE( xycePtr != NULL );
  int numArgs = 2;
  char * cmdLineArgs[numArgs];
  std::string xyceBin = "XyceTests";
  std::string netlist = "TestNetlist1.cir";
  cmdLineArgs[0] = const_cast<char *>(xyceBin.c_str());
  cmdLineArgs[1] = const_cast<char *>(netlist.c_str());
  Xyce::Circuit::Simulator::RunStatus status = xycePtr->initialize( numArgs, cmdLineArgs);
  EXPECT_EQ( status, Xyce::Circuit::Simulator::RunStatus::SUCCESS );
  double finalSimTime = xycePtr->getFinalTime();
  EXPECT_EQ( finalSimTime, 1.0 );
  std::vector<std::string> dacNames;
  bool circuitHasDACs = xycePtr->getDACDeviceNames(dacNames);
  EXPECT_FALSE( circuitHasDACs);
  EXPECT_EQ( dacNames.size(), 0);
  /*
  // run this simulation in several sub-stesps
  const int numSteps=100;
  for(auto i = 0; i<numSteps; i++)
  {
    // simToTime must be greater than zero.  passing zero in will cause Xyce to abort.
    double simToTime = (i+1)*finalSimTime/numSteps;
    double actualTime=0.0;
    bool stepResult = xycePtr->simulateUntil(simToTime, actualTime);
    EXPECT_TRUE(stepResult);
    // for this simple circuit we expect the simToTime and actualTime to be equal 
    EXPECT_EQ( simToTime, actualTime );
    double reportedTime = xycePtr->getTime();
    EXPECT_EQ( reportedTime, actualTime );
    bool isSimComplete = xycePtr->simulationComplete();
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
  */
  status = xycePtr->finalize();
  EXPECT_EQ( status, Xyce::Circuit::Simulator::RunStatus::SUCCESS );
  delete xycePtr;
}


TEST ( XyceSimulator, DACDeviceNamesTestNetlist3 )
{
  Xyce::Circuit::Simulator * xycePtr = NULL;
  xycePtr = new Xyce::Circuit::Simulator();
  EXPECT_TRUE( xycePtr != NULL );
  int numArgs = 2;
  char * cmdLineArgs[numArgs];
  std::string xyceBin = "XyceTests";
  std::string netlist = "TestNetlist3.cir";
  cmdLineArgs[0] = const_cast<char *>(xyceBin.c_str());
  cmdLineArgs[1] = const_cast<char *>(netlist.c_str());
  Xyce::Circuit::Simulator::RunStatus status = xycePtr->initialize( numArgs, cmdLineArgs);
  EXPECT_EQ( status, Xyce::Circuit::Simulator::RunStatus::SUCCESS );
  double finalSimTime = xycePtr->getFinalTime();
  EXPECT_EQ( finalSimTime, 20.0e-4 );
  std::vector<std::string> dacNames;
  bool circuitHasDACs = xycePtr->getDACDeviceNames(dacNames);
  EXPECT_TRUE( circuitHasDACs);
  EXPECT_EQ( dacNames.size(), 1);
  /*
  // run this simulation in several sub-stesps
  const int numSteps=100;
  for(auto i = 0; i<numSteps; i++)
  {
    // simToTime must be greater than zero.  passing zero in will cause Xyce to abort.
    double simToTime = (i+1)*finalSimTime/numSteps;
    double actualTime=0.0;
    bool stepResult = xycePtr->simulateUntil(simToTime, actualTime);
    EXPECT_TRUE(stepResult);
    // for this simple circuit we expect the simToTime and actualTime to be equal 
    EXPECT_EQ( simToTime, actualTime );
    double reportedTime = xycePtr->getTime();
    EXPECT_EQ( reportedTime, actualTime );
    bool isSimComplete = xycePtr->simulationComplete();
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
  */
  status = xycePtr->finalize();
  EXPECT_EQ( status, Xyce::Circuit::Simulator::RunStatus::SUCCESS );
  delete xycePtr;
}


TEST ( XyceSimulator, ADCMapTestNetlist1 )
{
  Xyce::Circuit::Simulator * xycePtr = NULL;
  xycePtr = new Xyce::Circuit::Simulator();
  EXPECT_TRUE( xycePtr != NULL );
  int numArgs = 2;
  char * cmdLineArgs[numArgs];
  std::string xyceBin = "XyceTests";
  std::string netlist = "TestNetlist1.cir";
  cmdLineArgs[0] = const_cast<char *>(xyceBin.c_str());
  cmdLineArgs[1] = const_cast<char *>(netlist.c_str());
  Xyce::Circuit::Simulator::RunStatus status = xycePtr->initialize( numArgs, cmdLineArgs);
  EXPECT_EQ( status, Xyce::Circuit::Simulator::RunStatus::SUCCESS );
  double finalSimTime = xycePtr->getFinalTime();
  EXPECT_EQ( finalSimTime, 1.0 );
  std::vector<std::string> dacNames;
  bool circuitHasDACs = xycePtr->getDACDeviceNames(dacNames);
  EXPECT_FALSE( circuitHasDACs);
  EXPECT_EQ( dacNames.size(), 0);
  std::map<std::string,std::map<std::string,double> > theADCMap;
  bool circuitHasADCs = xycePtr->getADCMap( theADCMap );
  EXPECT_FALSE( circuitHasADCs);
  EXPECT_EQ(theADCMap.size(), 0);
  /*
  // run this simulation in several sub-stesps
  const int numSteps=100;
  for(auto i = 0; i<numSteps; i++)
  {
    // simToTime must be greater than zero.  passing zero in will cause Xyce to abort.
    double simToTime = (i+1)*finalSimTime/numSteps;
    double actualTime=0.0;
    bool stepResult = xycePtr->simulateUntil(simToTime, actualTime);
    EXPECT_TRUE(stepResult);
    // for this simple circuit we expect the simToTime and actualTime to be equal 
    EXPECT_EQ( simToTime, actualTime );
    double reportedTime = xycePtr->getTime();
    EXPECT_EQ( reportedTime, actualTime );
    bool isSimComplete = xycePtr->simulationComplete();
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
  */
  status = xycePtr->finalize();
  EXPECT_EQ( status, Xyce::Circuit::Simulator::RunStatus::SUCCESS );
  delete xycePtr;
}

TEST ( XyceSimulator, ADCMapTestNetlist3 )
{
  Xyce::Circuit::Simulator * xycePtr = NULL;
  xycePtr = new Xyce::Circuit::Simulator();
  EXPECT_TRUE( xycePtr != NULL );
  int numArgs = 2;
  char * cmdLineArgs[numArgs];
  std::string xyceBin = "XyceTests";
  std::string netlist = "TestNetlist3.cir";
  cmdLineArgs[0] = const_cast<char *>(xyceBin.c_str());
  cmdLineArgs[1] = const_cast<char *>(netlist.c_str());
  Xyce::Circuit::Simulator::RunStatus status = xycePtr->initialize( numArgs, cmdLineArgs);
  EXPECT_EQ( status, Xyce::Circuit::Simulator::RunStatus::SUCCESS );
  double finalSimTime = xycePtr->getFinalTime();
  EXPECT_EQ( finalSimTime, 20.0e-4 );
  std::vector<std::string> dacNames;
  bool circuitHasDACs = xycePtr->getDACDeviceNames(dacNames);
  EXPECT_TRUE( circuitHasDACs);
  EXPECT_EQ( dacNames.size(), 1);
  std::map<std::string,std::map<std::string,double> > theADCMap;
  bool circuitHasADCs = xycePtr->getADCMap( theADCMap );
  EXPECT_TRUE( circuitHasADCs);
  EXPECT_EQ(theADCMap.size(), 2);
  /*
  // run this simulation in several sub-stesps
  const int numSteps=100;
  for(auto i = 0; i<numSteps; i++)
  {
    // simToTime must be greater than zero.  passing zero in will cause Xyce to abort.
    double simToTime = (i+1)*finalSimTime/numSteps;
    double actualTime=0.0;
    bool stepResult = xycePtr->simulateUntil(simToTime, actualTime);
    EXPECT_TRUE(stepResult);
    // for this simple circuit we expect the simToTime and actualTime to be equal 
    EXPECT_EQ( simToTime, actualTime );
    double reportedTime = xycePtr->getTime();
    EXPECT_EQ( reportedTime, actualTime );
    bool isSimComplete = xycePtr->simulationComplete();
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
  */
  status = xycePtr->finalize();
  EXPECT_EQ( status, Xyce::Circuit::Simulator::RunStatus::SUCCESS );
  delete xycePtr;
}

TEST ( XyceSimulator, SetDeviceParamTestNetlist2 )
{
  Xyce::Circuit::Simulator * xycePtr = NULL;
  xycePtr = new Xyce::Circuit::Simulator();
  EXPECT_TRUE( xycePtr != NULL );
  int numArgs = 2;
  char * cmdLineArgs[numArgs];
  std::string xyceBin = "XyceTests";
  std::string netlist = "TestNetlist2.cir";
  cmdLineArgs[0] = const_cast<char *>(xyceBin.c_str());
  cmdLineArgs[1] = const_cast<char *>(netlist.c_str());
  Xyce::Circuit::Simulator::RunStatus status = xycePtr->initialize( numArgs, cmdLineArgs);
  EXPECT_EQ( status, Xyce::Circuit::Simulator::RunStatus::SUCCESS );
  double finalSimTime = xycePtr->getFinalTime();
  EXPECT_EQ( finalSimTime, 1.0 );
  
  // this should fail as TOMP is not a simulation parameter in this netlist 
  bool setParamResult = xycePtr->setCircuitParameter( "TOMP", 10.0);
  EXPECT_FALSE( setParamResult );
  
  // this should pass 
  setParamResult = xycePtr->setCircuitParameter( "TEMP", 10.0);
  EXPECT_TRUE( setParamResult );
  
  // this should passs
  setParamResult = xycePtr->setCircuitParameter( "R1:R", 10.0);
  EXPECT_TRUE( setParamResult );
  
  // run this simulation in several sub-stesps
  const int numSteps=100;
  for(auto i = 0; i<numSteps; i++)
  {
    // simToTime must be greater than zero.  passing zero in will cause Xyce to abort.
    double simToTime = (i+1)*finalSimTime/numSteps;
    double actualTime=0.0;
    bool stepResult = xycePtr->simulateUntil(simToTime, actualTime);
    EXPECT_TRUE(stepResult);
    // for this simple circuit we expect the simToTime and actualTime to be equal 
    EXPECT_EQ( simToTime, actualTime );
    double reportedTime = xycePtr->getTime();
    EXPECT_EQ( reportedTime, actualTime );
    if( i==(numSteps/2))
    {
      // in the middle of the run change the resistance on resistor R1
      setParamResult = xycePtr->setCircuitParameter( "R1:R", 1000.0);
      EXPECT_TRUE( setParamResult );
      setParamResult = xycePtr->setCircuitParameter( "TEMP", -50.0);
      EXPECT_TRUE( setParamResult );
    }
    bool isSimComplete = xycePtr->simulationComplete();
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
  status = xycePtr->finalize();
  EXPECT_EQ( status, Xyce::Circuit::Simulator::RunStatus::SUCCESS );
  delete xycePtr;
}


TEST ( XyceSimulator, GetCircuitValuesTestNetlist2 )
{
  Xyce::Circuit::Simulator * xycePtr = NULL;
  xycePtr = new Xyce::Circuit::Simulator();
  EXPECT_TRUE( xycePtr != NULL );
  int numArgs = 2;
  char * cmdLineArgs[numArgs];
  std::string xyceBin = "XyceTests";
  std::string netlist = "TestNetlist2.cir";
  cmdLineArgs[0] = const_cast<char *>(xyceBin.c_str());
  cmdLineArgs[1] = const_cast<char *>(netlist.c_str());
  Xyce::Circuit::Simulator::RunStatus status = xycePtr->initialize( numArgs, cmdLineArgs);
  EXPECT_EQ( status, Xyce::Circuit::Simulator::RunStatus::SUCCESS );
  double finalSimTime = xycePtr->getFinalTime();
  EXPECT_EQ( finalSimTime, 1.0 );

  double circuitParamValue = 0.0;
  bool getParamResult = xycePtr->getCircuitValue( "R1:R",circuitParamValue);
  EXPECT_EQ( circuitParamValue, 50.0 );
  getParamResult = xycePtr->getCircuitValue( "TEMP", circuitParamValue);
  EXPECT_EQ( circuitParamValue, 27.0 );
  
  // run this simulation in several sub-stesps
  const int numSteps=100;
  for(auto i = 0; i<numSteps; i++)
  {
    // simToTime must be greater than zero.  passing zero in will cause Xyce to abort.
    double simToTime = (i+1)*finalSimTime/numSteps;
    double actualTime=0.0;
    bool stepResult = xycePtr->simulateUntil(simToTime, actualTime);
    EXPECT_TRUE(stepResult);
    // for this simple circuit we expect the simToTime and actualTime to be equal 
    EXPECT_EQ( simToTime, actualTime );
    double reportedTime = xycePtr->getTime();
    EXPECT_EQ( reportedTime, actualTime );
    if( i==(numSteps/2))
    {
      // in the middle of the run change the resistance on resistor R1
      bool setParamResult = xycePtr->setCircuitParameter( "R1:R", 1000.0);
      EXPECT_TRUE( setParamResult );
      setParamResult = xycePtr->setCircuitParameter( "TEMP", -50.0);
      EXPECT_TRUE( setParamResult );
    }
    if( i==((numSteps/2)+1))
    {
      // in the middle of the run change the resistance on resistor R1
      circuitParamValue = 0.0;
      getParamResult = xycePtr->getCircuitValue( "R1:R",circuitParamValue);
      EXPECT_EQ( circuitParamValue, 1000.0 );
      getParamResult = xycePtr->getCircuitValue( "TEMP", circuitParamValue);
      EXPECT_EQ( circuitParamValue, -50.0 );
    }
    
    bool isSimComplete = xycePtr->simulationComplete();
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
  getParamResult = xycePtr->getCircuitValue( "MAXV2",circuitParamValue);
  EXPECT_NEAR( circuitParamValue, 6.666667e-01, 1.0e-7 );
  getParamResult = xycePtr->getCircuitValue( "MINV2",circuitParamValue);
  EXPECT_NEAR( circuitParamValue, -6.666667e-01,  1.0e-7 );
  
  status = xycePtr->finalize();
  EXPECT_EQ( status, Xyce::Circuit::Simulator::RunStatus::SUCCESS );
  delete xycePtr;
}

//-------------------------------------------------------------------------------
int main (int argc, char **argv)
{
  testing::InitGoogleTest(&argc, argv);

  return RUN_ALL_TESTS();
}

