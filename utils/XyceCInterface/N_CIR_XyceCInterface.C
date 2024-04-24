//-------------------------------------------------------------------------
//   Copyright 2002-2024 National Technology & Engineering Solutions of
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

//-----------------------------------------------------------------------------
//
// Purpose        : C callable functions to link to Xyce's C++ interface
//
// Special Notes  : These functions expose Xyce's C++ interface to 
//                  other programs (i.e. Python) that can dynamically 
//                  link to a C library (via ctypes in Python)
//
//
//-----------------------------------------------------------------------------


// ---------- Standard Includes ----------
#include <iostream>
#include <ctime>
#include <cstring>

// ----------   Xyce Includes   ----------
#include <Xyce_config.h>
#include <N_CIR_GenCouplingSimulator.h>
#include <N_DEV_Device.h>
//#include <N_DEV_Algorithm.h>
#include <N_DEV_ADC.h>

#include <N_CIR_XyceCInterface.h>

//-----------------------------------------------------------------------------
// Function      : xyce_open
// Purpose       : Create a pointer to an Xyce::Circuit::GenCouplingSimulator object (whose class name 
//                 is Xyce::Circuit::Simulator). This function must be called 
//                 before any of the other functions in this file are used.
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Electrical Models and Simulation
// Creation Date : 3/01/2018
//-----------------------------------------------------------------------------
void xyce_open( void ** ptr)
{
  Xyce::Circuit::GenCouplingSimulator * xycePtr = new Xyce::Circuit::GenCouplingSimulator();
  *ptr = xycePtr;
}

//-----------------------------------------------------------------------------
// Function      : xyce_close
// Purpose       : Call the Xyce::Circuit::Simulator::finalize function via a 
//                 pointer to an Xyce::Circuit::GenCouplingSimulator object. It then cleans up that 
//                 pointer (that typically was created with the xyce_open 
//                 function).
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Electrical Models and Simulation
// Creation Date : 3/01/2018
//-----------------------------------------------------------------------------
void xyce_close( void** ptr)
{
  Xyce::Circuit::GenCouplingSimulator * xycePtr = static_cast<Xyce::Circuit::GenCouplingSimulator *>( *ptr );
  xycePtr->finalize();
  delete xycePtr;
  *ptr = 0;
}

//-----------------------------------------------------------------------------
// Function      : xyce_initialize
// Purpose       : Call the Xyce::Circuit::Simulator::initialize function via a 
//                 pointer to an Xyce::Circuit::GenCouplingSimulator object. A typical use case is to 
//                 create that pointer with the xyce_open() function.  The 
//                 return values are Xyce::Circuit::Simulator::RunStatus enum 
//                 values.
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Electrical Models and Simulation
// Creation Date : 3/01/2018
//-----------------------------------------------------------------------------
int xyce_initialize_early( void** ptr, int argc, char ** argv )
{
  Xyce::Circuit::GenCouplingSimulator * xycePtr = static_cast<Xyce::Circuit::GenCouplingSimulator *>( *ptr );
  return( xycePtr->initializeEarly( argc, argv ) );
}

int xyce_initialize_late( void** ptr )
{
  Xyce::Circuit::GenCouplingSimulator * xycePtr = static_cast<Xyce::Circuit::GenCouplingSimulator *>( *ptr );
  return( xycePtr->initializeLate() );
}

int xyce_initialize( void** ptr, int argc, char ** argv )
{
  int run_status = xyce_initialize_early (ptr, argc, argv);
  if ( run_status != Xyce::Circuit::Simulator::SUCCESS ) {
    return run_status;
  }
  return xyce_initialize_late (ptr);
}
    

//-----------------------------------------------------------------------------
// Function      : xyce_runSimulation
// Purpose       : Call the Xyce::Circuit::Simulator::runSimulation function via 
//                 a pointer to an Xyce::Circuit::GenCouplingSimulator object. A typical use case is to 
//                 create that pointer with the xyce_open() function, and to then 
//                 initialize the Xyce::Circuit::GenCouplingSimulator object with the xyce_initialize() 
//                 function before using this function.
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Electrical Models and Simulation
// Creation Date : 3/01/2018
//-----------------------------------------------------------------------------
int xyce_runSimulation(void ** ptr)
{
  Xyce::Circuit::GenCouplingSimulator * xycePtr = static_cast<Xyce::Circuit::GenCouplingSimulator *>( *ptr );
  return(xycePtr->runSimulation());
}

//-----------------------------------------------------------------------------
// Function      : xyce_simulateUntil
// Purpose       : Call the Xyce::Circuit::Simulator::simulateUntil function
//                 via a pointer to an Xyce::Circuit::GenCouplingSimulator object. A typical use case is 
//                 to create that pointer with the xyce_open() function, and to 
//                 then initialize the Xyce::Circuit::GenCouplingSimulator object with the xyce_initialize() 
//                 function before using this function.  This function will
//                 return 0 if the simulation failed, because of an issue such as
//                 a DCOP failure.  See the "Application Note: Mixed Signal 
//                 Simulation with Xyce" for a discussion of the parameters
//                 requestedUntilTime and completedUntilTime
// Special Notes : For C compatibility, this function uses double*, rather than
//                 double &, for completedUntilTime which essentially returns
//                 that value to the calling program.
// Scope         : public
// Creator       : Rich Schiek, SNL, Electrical Models and Simulation
// Creation Date : 3/01/2018
//-----------------------------------------------------------------------------
int xyce_simulateUntil(  void **ptr, double requestedUntilTime, double* completedUntilTime )
{
  Xyce::Circuit::GenCouplingSimulator * xycePtr = static_cast<Xyce::Circuit::GenCouplingSimulator *>( *ptr );
  double cTime;
  bool retVal = xycePtr->simulateUntil( requestedUntilTime,  cTime);
  *completedUntilTime =cTime;
  if (retVal)
    return 1;
  return 0;
}


//-----------------------------------------------------------------------------
// Function      : xyce_simulationComplete
// Purpose       : Call the Xyce::Circuit::Simulator::simulationComplete()
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Electrical Models and Simulation
// Creation Date : 9/28/2022
//-----------------------------------------------------------------------------
bool xyce_simulationComplete( void ** ptr)
{
  Xyce::Circuit::GenCouplingSimulator * xycePtr = static_cast<Xyce::Circuit::GenCouplingSimulator *>( *ptr );
  bool simCompleteFlag = xycePtr->simulationComplete();
  return simCompleteFlag;
}

//-----------------------------------------------------------------------------
// Function      : xyce_checkCircuitParameterExists
// Purpose       : Call the Xyce::Circuit::Simulator::checkCircuitParameterExists()
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Electrical Models and Simulation
// Creation Date : 9/28/2022
//-----------------------------------------------------------------------------
bool xyce_checkCircuitParameterExists(void **ptr, char * paramName )
{
  Xyce::Circuit::GenCouplingSimulator * xycePtr = static_cast<Xyce::Circuit::GenCouplingSimulator *>( *ptr );
  std::string sParamName( paramName );
  bool result = xycePtr->checkCircuitParameterExists( sParamName);
  return result;
}


//-----------------------------------------------------------------------------
// Function      : xyce_getTime
// Purpose       : Call the Xyce::Circuit::Simulator::getTime()
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Electrical Models and Simulation
// Creation Date : 9/28/2022
//-----------------------------------------------------------------------------
double xyce_getTime(void ** ptr)
{
  Xyce::Circuit::GenCouplingSimulator * xycePtr = static_cast<Xyce::Circuit::GenCouplingSimulator *>( *ptr );
  double simulatorTime = xycePtr->getTime();
  return simulatorTime;
}

//-----------------------------------------------------------------------------
// Function      : xyce_getFinalTime
// Purpose       : Call the Xyce::Circuit::Simulator::getFinalTime()
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Electrical Models and Simulation
// Creation Date : 9/28/2022
//-----------------------------------------------------------------------------
double xyce_getFinalTime(void ** ptr)
{
  Xyce::Circuit::GenCouplingSimulator * xycePtr = static_cast<Xyce::Circuit::GenCouplingSimulator *>( *ptr );
  double finalTime = xycePtr->getFinalTime();
  return finalTime;
}


//-----------------------------------------------------------------------------
// Function      : xyce_getNumDevices
// Purpose       : Call the Xyce::Circuit::Simulator::getDeviceNames function
//                 via a pointer to an Xyce::Circuit::GenCouplingSimulator object. A typical use case is
//                 to create that pointer with the xyce_open() function, and to
//                 then initialize the Xyce::Circuit::GenCouplingSimulator object with the xyce_initialize()
//                 function before using this function.  This function is then
//                 typically used to get the number for devices (of the requested
//                 type) before a subsequent call to xyce_getNumDevices(),
//                 xyce_getADCMap() or xyce_getDACDeviceNames() or is made.
//                 This allows the calling C (or Python) program to allocate
//                 the correct-sized arrays for their returned parameter(s).
//                 See SON Bug 1066 for more details.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical Models and Simulation
// Creation Date : 4/29/2019
//-----------------------------------------------------------------------------
int xyce_getNumDevices(void **ptr, char * modelGroupName, int* numDevNames, int* maxDevNameLength)
{
  Xyce::Circuit::GenCouplingSimulator * xycePtr = static_cast<Xyce::Circuit::GenCouplingSimulator *>( *ptr );
  const std::string modelGroupNameString( modelGroupName );
  std::vector<std::string> deviceNamesVec;
  bool retVal = xycePtr->getDeviceNames(modelGroupNameString , deviceNamesVec);

  *numDevNames = deviceNamesVec.size();
  *maxDevNameLength = 0;

  for( int i=0;i < deviceNamesVec.size(); i++ )
  {
    if (strlen(deviceNamesVec.at(i).c_str()) > *maxDevNameLength )
    {
      *maxDevNameLength = strlen(deviceNamesVec.at(i).c_str());
    }
  }
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : xyce_getTotalNumDevices
// Purpose       : Call the Xyce::Circuit::Simulator::getDeviceNames function
//                 via a pointer to an Xyce::Circuit::GenCouplingSimulator object. A typical use case is
//                 to create that pointer with the xyce_open() function, and to
//                 then initialize the Xyce::Circuit::GenCouplingSimulator object with the xyce_initialize()
//                 function before using this function.  This function is then
//                 typically used to get the total number for devices in the
//                 netlist before a subsequent call to xyce_getAllDeviceNames().
//                 This allows the calling C (or Python) program to allocate
//                 the correct-sized arrays for their returned parameter(s).
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical Models and Simulation
// Creation Date : 12/11/2019
//-----------------------------------------------------------------------------
int xyce_getTotalNumDevices(void **ptr, int* numDevNames, int* maxDevNameLength)
{
  Xyce::Circuit::GenCouplingSimulator * xycePtr = static_cast<Xyce::Circuit::GenCouplingSimulator *>( *ptr );
  std::vector<std::string> deviceNamesVec;
  int retVal = xycePtr->getAllDeviceNames(deviceNamesVec);

  *numDevNames = deviceNamesVec.size();
  *maxDevNameLength = 0;

  for( int i=0;i < deviceNamesVec.size(); i++ )
  {
    if (strlen(deviceNamesVec.at(i).c_str()) > *maxDevNameLength )
    {
      *maxDevNameLength = strlen(deviceNamesVec.at(i).c_str());
    }
  }
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : xyce_getDeviceNames
// Purpose       : Call the Xyce::Circuit::Simulator::getDeviceNames function
//                 via a pointer to an Xyce::Circuit::GenCouplingSimulator object. A typical use case is 
//                 to create that pointer with the xyce_open() function, and to 
//                 then initialize the Xyce::Circuit::GenCouplingSimulator object with the xyce_initialize() 
//                 function before using this function.  The modelGroupName can
//                 be any valid model group name (e.g., "M" for MOSFETs, "Q" for 
//                 BJTs, etc.).  As other examples, for Y devices, "YDAC" and not "DAC"
//                 should be used.  For U devices, "BUF" and not "UBUF" should be used.   
//                 The deviceNames parameter will then contain the fully-qualifed 
//                 names (if any) of the devices from that model group in the 
//                 netlist.  Return 1 if any device of the requested model group
//                 exists in the netlist.  Return 0 otherwise.
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Electrical Models and Simulation
// Creation Date : 3/01/2018
//-----------------------------------------------------------------------------
int xyce_getDeviceNames(void ** ptr, char * modelGroupName, int* numDevNames, char ** deviceNames)
{
  Xyce::Circuit::GenCouplingSimulator * xycePtr = static_cast<Xyce::Circuit::GenCouplingSimulator *>( *ptr );
  const std::string modelGroupNameString( modelGroupName );
  std::vector<std::string> deviceNamesVec;
  bool retVal = xycePtr->getDeviceNames(modelGroupNameString , deviceNamesVec);
  
  *numDevNames = deviceNamesVec.size();
  for( int i=0;i < deviceNamesVec.size(); i++ )
  {
    std::strcpy( (deviceNames)[i], deviceNamesVec.at(i).c_str() );
  }
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : xyce_getAllDeviceNames
// Purpose       : Call the Xyce::Circuit::Simulator::getDeviceNames function
//                 via a pointer to an Xyce::Circuit::GenCouplingSimulator object. A typical use case is
//                 to create that pointer with the xyce_open() function, and to
//                 then initialize the Xyce::Circuit::GenCouplingSimulator object with the xyce_initialize()
//                 function before using this function.  The deviceNames parameter
//                 will then contain the fully-qualifed names of all the devices
//                 in the netlist.  Return 0 if there are no devices in the netlist.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 12/11/2019
//-----------------------------------------------------------------------------
int xyce_getAllDeviceNames(void ** ptr, int* numDevNames, char ** deviceNames)
{
  Xyce::Circuit::GenCouplingSimulator * xycePtr = static_cast<Xyce::Circuit::GenCouplingSimulator *>( *ptr );
  std::vector<std::string> deviceNamesVec;
  bool retVal = xycePtr->getAllDeviceNames(deviceNamesVec);

  *numDevNames = deviceNamesVec.size();
  for( int i=0;i < deviceNamesVec.size(); i++ )
  {
    std::strcpy( (deviceNames)[i], deviceNamesVec.at(i).c_str() );
  }
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : xyce_getDACDeviceNames
// Purpose       : Call the Xyce::Circuit::Simulator::getDACDeviceNames 
//                 function via a pointer to an Xyce::Circuit::GenCouplingSimulator object. So, it is
//                 basically a specialized version of xyce_getDeviceNames. 
//                 A typical use case is to create that pointer with the 
//                 xyce_open() function, and to then initialize the Xyce::Circuit::GenCouplingSimulator 
//                 object with the xyce_initialize() function before using this 
//                 function.  Return 1 if any YDAC devices exist in the netlist. 
//                 Return 0 otherwise.
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Electrical Models and Simulation
// Creation Date : 3/01/2018
//-----------------------------------------------------------------------------
int xyce_getDACDeviceNames(void ** ptr, int* numDevNames, char ** deviceNames)
{
  Xyce::Circuit::GenCouplingSimulator * xycePtr = static_cast<Xyce::Circuit::GenCouplingSimulator *>( *ptr );
  std::vector<std::string> deviceNamesVec;
  bool retVal= xycePtr->getDACDeviceNames( deviceNamesVec );

  *numDevNames = deviceNamesVec.size();
  for( int i=0;i < deviceNamesVec.size(); i++ )
  {
    std::strcpy( (deviceNames)[i], deviceNamesVec.at(i).c_str() );
  }
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : xyce_checkDeviceParamName
// Purpose       : Call the Xyce::Circuit::Simulator::checkDeviceParamName function
//                 via a pointer to an Xyce::Circuit::GenCouplingSimulator object.  A typical use case is
//                 to create that pointer with the xyce_open() function, and to
//                 then initialize the Xyce::Circuit::GenCouplingSimulator object with the xyce_initialize()
//                 function before using this function.  The full_param_name
//                 should be identical to that used on a .PRINT line; so X1:R1:R
//                 for the resistance (R) of device R1 that is in subcircuit X1.
//                 This function returns 1 if full_param__name is a valid
//                 parameter for a device that exists in the netlist. It returns
//                 0, otherwise.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 12/18/2019
//-----------------------------------------------------------------------------
int xyce_checkDeviceParamName(void **ptr, char* full_param_name)
{
  Xyce::Circuit::GenCouplingSimulator * xycePtr = static_cast<Xyce::Circuit::GenCouplingSimulator *>( *ptr );
  std::string full_param_nameStr( full_param_name);
  int result = xycePtr->checkDeviceParamName( full_param_nameStr );
  return result;
}

//-----------------------------------------------------------------------------
// Function      : xyce_getDeviceParamVal
// Purpose       : Call the Xyce::Circuit::Simulator::getDeviceParamVal function
//                 via a pointer to an Xyce::Circuit::GenCouplingSimulator object.  A typical use case is
//                 to create that pointer with the xyce_open() function, and to
//                 then initialize the Xyce::Circuit::GenCouplingSimulator object with the xyce_initialize()
//                 function before using this function.  The full_param_name
//                 should be identical to that used on a .PRINT line; so X1:R1:R
//                 for the resistance (R) of device R1 that is in subcircuit X1.
//                 This function returns 1 if full_param__name is a valid
//                 parameter for a device that exists in the netlist. It returns
//                 0, otherwise. For a return value of 0, the value parameter is
//                 also set to 0.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 12/18/2019
//-----------------------------------------------------------------------------
int xyce_getDeviceParamVal(void **ptr, char* full_param_name, double* value)
{
  Xyce::Circuit::GenCouplingSimulator * xycePtr = static_cast<Xyce::Circuit::GenCouplingSimulator *>( *ptr );
  std::string full_param_nameStr( full_param_name);
  double pVal;
  int result = xycePtr->getDeviceParamVal( full_param_nameStr, pVal );
  *value = pVal;
  if( result == 0)
  {
    // device or parameter wasn't found.  Need to set the return parameter 
    // value to zero as that is the expected behavior of this function.
    *value = 0.0;
  }
  return result;
}

//-----------------------------------------------------------------------------
// Function      : xyce_getNumAdjNodesForDevice
// Purpose       : Call the Xyce::Circuit::Simulator::getNumAdjNodesForDevice function
//                 via a pointer to an Xyce::Circuit::GenCouplingSimulator object.  A typical use case is
//                 to create that pointer with the xyce_open() function, and to
//                 then initialize the Xyce::Circuit::GenCouplingSimulator object with the xyce_initialize()
//                 function before using this function.  The device name should be
//                 fully-qualified; so X1:R1 for device R1 in subcircuit X1.
//                 This function returns 1 if the device that exists in the netlist.
//                 It returns 0, otherwise. For a return value of 0, the numAdjNodes
//                 parameter is also set to 0.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 1/7/2020
//-----------------------------------------------------------------------------
int xyce_getNumAdjNodesForDevice(void **ptr, char * deviceName, int* numAdjNodes)
{
  Xyce::Circuit::GenCouplingSimulator * xycePtr = static_cast<Xyce::Circuit::GenCouplingSimulator *>( *ptr );
  std::string deviceName_str(deviceName);
  int mNumAdjNodes;

  int result = xycePtr->getNumAdjNodesForDevice(deviceName_str, mNumAdjNodes);
  *numAdjNodes = mNumAdjNodes;

  return result;
}

//-----------------------------------------------------------------------------
// Function      : xyce_getAdjGIDsForDevice
// Purpose       : Call the Xyce::Circuit::Simulator::getAdjGIDsForDevice function
//                 via a pointer to an Xyce::Circuit::GenCouplingSimulator object.  A typical use case is
//                 to create that pointer with the xyce_open() function, and to
//                 then initialize the Xyce::Circuit::GenCouplingSimulator object with the xyce_initialize()
//                 function before using this function.  The device name should be
//                 fully-qualified; so X1:R1 for device R1 in subcircuit X1.
//                 This function returns 1 if the device that exists in the netlist.
//                 It returns 0, otherwise. For a return value of 0, the adjGID
//                 parameter is empty.
// Special Notes : The GID of the ground node is -1.  The other GIDs are
//                 non-negative integers.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 1/7/2020
//-----------------------------------------------------------------------------
int xyce_getAdjGIDsForDevice(void **ptr, char * deviceName, int* numAdjNodes, int* adjGIDs)
{
  Xyce::Circuit::GenCouplingSimulator * xycePtr = static_cast<Xyce::Circuit::GenCouplingSimulator *>( *ptr );
  std::string deviceName_str(deviceName);
  std::vector<int> gids;

  int result = xycePtr->getAdjGIDsForDevice(deviceName_str, gids);
  *numAdjNodes = gids.size();

  for( int i=0; i<gids.size(); i++ )
  {
    adjGIDs[i] = gids[i];
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : xyce_getADCmap
// Purpose       : Call the Xyce::Circuit::Simulator::getADCMap function via
//                 a pointer to an Xyce::Circuit::GenCouplingSimulator object. A typical use case is to
//                 create that pointer with the xyce_open() function, and to 
//                 then initialize the Xyce::Circuit::GenCouplingSimulator object with the xyce_initialize() 
//                 function before using this function. This function will 
//                 return 1 if any YADC devices exist in the netlist. It will 
//                 return 0 otherwise. 
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical Models and Simulation
// Creation Date : 12/19/2018
//-----------------------------------------------------------------------------
int xyce_getADCMap(void ** ptr, int* numADCnames, char ** ADCnames, int *  widths, double * resistances,
                   double * upperVLimits, double * lowerVLimits, double * settlingTimes)
{
  Xyce::Circuit::GenCouplingSimulator * xycePtr = static_cast<Xyce::Circuit::GenCouplingSimulator *>( *ptr );
  std::map<std::string, std::map<std::string, double> > ADCMap;

  int status = xycePtr->getADCMap(ADCMap);
  *numADCnames = ADCMap.size();
  
  std::map<std::string, std::map<std::string, double> >::const_iterator it;
  int i=0;
  for(it=ADCMap.begin(); it!=ADCMap.end(); it++, i++ )
  {
    std::string name = (*it).first;
    std::strcpy( (ADCnames)[i], name.c_str() );
    widths[i] = ADCMap[name]["width"];
    resistances[i] = ADCMap[name]["r"];
    upperVLimits[i] = ADCMap[name]["upperVoltageLimit"];
    lowerVLimits[i] = ADCMap[name]["lowerVoltageLimit"];
    settlingTimes[i] = ADCMap[name]["settlingTime"];
  }

  return status;
}

//-----------------------------------------------------------------------------
// Function      : xyce_updateTimeVoltagePairs
// Purpose       : Call the Xyce::Circuit::Simulator::updateTimeVoltagePairs 
//                 function via a pointer to an Xyce::Circuit::GenCouplingSimulator object. A typical 
//                 use case is to create that pointer with the xyce_open() 
//                 function, and to then initialize the Xyce::Circuit::GenCouplingSimulator object with 
//                 the xyce_initialize() function before using this function. 
//                 See the "Application Note: Mixed Signal Simulation with 
//                 Xyce" for a discussion of how to use this function. This 
//                 function method will return 1 if the time-voltage pairs for 
//                 the specified DACname were successfully updated. It will 
//                 return 0, otherwise.
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Electrical Models and Simulation
// Creation Date : 3/01/2018
//-----------------------------------------------------------------------------
int xyce_updateTimeVoltagePairs(void ** ptr, char * DACname, int numPoints, double * timeArray, double * voltageArray)
{
  Xyce::Circuit::GenCouplingSimulator * xycePtr = static_cast<Xyce::Circuit::GenCouplingSimulator *>( *ptr );

  std::string name(DACname);
  std::map< std::string, std::vector<std::pair<double,double> > *> timePointsMap;
  std::vector<std::pair<double,double> > storageVec;
  if( numPoints > 0) 
  {
    for( int i=0 ; i<numPoints; i++ ) 
    {
      storageVec.push_back( std::pair<double,double>( timeArray[i], voltageArray[i] ) );
    }
  }
  timePointsMap[ name ] = &storageVec;
  /*
  std::cout << "name= \"" << name << "\"" << std::endl;
  for( int i=0; i<numPoints; i++ )
    std::cout << "storageVec[ " << i << "] = " << storageVec[i].first << ", " << storageVec[i].second << std::endl;
  */
  return( xycePtr->updateTimeVoltagePairs(timePointsMap) );
}

//-----------------------------------------------------------------------------
// Function      : xyce_getTimeVoltagePairsADC 
// Purpose       : Call the Xyce::Circuit::Simulator::getTimeVoltagePairs
//                 function via a pointer to an Xyce::Circuit::GenCouplingSimulator object. A typical 
//                 use case is to create that pointer with the xyce_open() 
//                 function, and to then initialize the Xyce::Circuit::GenCouplingSimulator object with 
//                 the xyce_initialize() function before using this function. 
//                 See the "Application Note: Mixed Signal Simulation with 
//                 Xyce" for a discussion of how to use this function. This 
//                 function returns 1 if any ADC devices were found in the
//                 simulation.  It returns 0, otherwise.
// Special Notes : This function is the least "mature" of the implemented 
//                 functions.  The Application Note has more details on that.
// Scope         : public
// Creator       : Rich Schiek, SNL, Electrical Models and Simulation
// Creation Date : 3/01/2018
//-----------------------------------------------------------------------------
int xyce_getTimeVoltagePairsADC( void** ptr, int * numADCnames, char ** ADCnames, int * numPoints, double ** timeArray, double ** voltageArray )
{
  Xyce::Circuit::GenCouplingSimulator * xycePtr = static_cast<Xyce::Circuit::GenCouplingSimulator *>( *ptr );
 
  std::map< std::string, std::vector< std::pair<double,double> > > timeVoltageUpdateMap;
  int status = xycePtr->getTimeVoltagePairs( timeVoltageUpdateMap );
  
  std::map< std::string, std::vector< std::pair<double,double> > >::iterator cMap = timeVoltageUpdateMap.begin();
  std::map< std::string, std::vector< std::pair<double,double> > >::iterator eMap = timeVoltageUpdateMap.end();
  int ADCnum = 0;
  int PTnum = 0;
  while( cMap != eMap )
  {
    strcpy( ADCnames[ADCnum], (cMap->first).c_str() );
    std::vector< std::pair<double,double> > dataVec = cMap->second;
    int PTnumThisADC = 0;
    for( int j = 0; j<dataVec.size(); j++ )
    {
      timeArray[ADCnum][j] = dataVec[j].first;
      voltageArray[ADCnum][j] = dataVec[j].second;
      PTnumThisADC++;
      //Xyce::dout() << "timeArray and Voltage Array for ADC " << ADCnum << " for " << j
      //     << " are "  << dataVec[j].first << " and " << dataVec[j].second << std::endl;
    }
    if (PTnumThisADC > PTnum) PTnum = PTnumThisADC;
    ADCnum++;
    ++cMap;
  }
  *numPoints=PTnum;
  *numADCnames=ADCnum;
  
  return status;
}

//-----------------------------------------------------------------------------
// Function      : xyce_getTimeVoltagePairsADCsz 
// Purpose       : Call the Xyce::Circuit::Simulator::getTimeVoltagePairsSz()
//                 function via a pointer to an Xyce::Circuit::GenCouplingSimulator object.
// Scope         : public
// Creator       : 
// Creation Date : 11/11/2021
//-----------------------------------------------------------------------------
int xyce_getTimeVoltagePairsADCsz( void** ptr, int *maxPoints )
{
  Xyce::Circuit::GenCouplingSimulator * xycePtr = static_cast<Xyce::Circuit::GenCouplingSimulator *>( *ptr );

  int maxPts;
  
  int status = xycePtr->getTimeVoltagePairsSz( maxPts );

  *maxPoints = maxPts;

  return status;
}


//-----------------------------------------------------------------------------
// Function      : xyce_getTimeVoltagePairsADCLimitData
// Purpose       : Call the Xyce::Circuit::Simulator::getTimeVoltagePairs
//                 function via a pointer to an Xyce::Circuit::GenCouplingSimulator object but limit the
//                 returned data assuming the maximum dimensions of the arras are:
//                 ADCnamesArray[maxNumADCnames][maxNameLength]
//                 timeArray[maxNumADCnames][maxNumPoints]
//                 voltageArray[maxNumADCnames][maxNumPoints]
//                 numPointsArray[maxNumADCnames]
//
//                 numPointsArray[i] is the number of time and voltage points 
//                 for the i'th ADC in the timeArray[i][] and voltageArray[i][]
//
//                 Using numPointsArray allows the calling program to know the 
//                 number of valid points in the time and voltage arrays. 
//                 A typical use case is to create that pointer with the 
//                 xyce_open() function, and to then initialize the Xyce::Circuit::GenCouplingSimulator 
//                 object with the xyce_initialize() function before using this 
//                 function.  Also, precallocation of the arrays ADCnamesArray,
//                 timeArray, voltageArray and numPointsArray is REQUIRED.
//
//                 See the "Application Note: Mixed Signal Simulation with 
//                 Xyce" for a discussion of how to use this function. This 
//                 function returns 1 if any ADC devices were found in the
//                 simulation.  It returns 0, otherwise.
//
// Scope         : public
// Creator       : Rich Schiek, SNL, Electrical Models and Simulation
// Creation Date : 12/02/2021
//-----------------------------------------------------------------------------
int xyce_getTimeVoltagePairsADCLimitData( void** ptr, 
  const int maxNumADCnames, const int maxNameLength, const int maxNumPoints,
  int * numADCnames, char ** ADCnamesArray, int * numPointsArray, double ** timeArray, double ** voltageArray )
{
  Xyce::Circuit::GenCouplingSimulator * xycePtr = static_cast<Xyce::Circuit::GenCouplingSimulator *>( *ptr );
  
  std::map< std::string, std::vector< std::pair<double,double> > > timeVoltageUpdateMap;
  int status = xycePtr->getTimeVoltagePairs( timeVoltageUpdateMap );
  
  std::map< std::string, std::vector< std::pair<double,double> > >::iterator cMap = timeVoltageUpdateMap.begin();
  std::map< std::string, std::vector< std::pair<double,double> > >::iterator eMap = timeVoltageUpdateMap.end();
  int ADCnum = 0;
  while( (cMap != eMap) && (ADCnum <maxNumADCnames) )
  {
    strncpy( ADCnamesArray[ADCnum], (cMap->first).c_str(), maxNameLength );
    std::vector< std::pair<double,double> > dataVec = cMap->second;
    // return the size of the array or maxNumPoints whichever is smaller.
    if(dataVec.size() < maxNumPoints)
    {
      numPointsArray[ADCnum]= dataVec.size();
      int j=0;
      while((j<dataVec.size()))
      {
        timeArray[ADCnum][j] = dataVec[j].first;
        voltageArray[ADCnum][j] = dataVec[j].second;
        j++;
      }
    }
    else
    {
      // number of points is greater than allocated space.  
      // just return the LAST points as they include the end time 
      // step and steps closest to that time. 
      numPointsArray[ADCnum]= maxNumPoints;
      int j=dataVec.size()-1;
      int k= maxNumPoints-1;
      while(k>=0)
      {
        timeArray[ADCnum][k] = dataVec[j].first;
        voltageArray[ADCnum][k] = dataVec[j].second;
        j--;
        k--;
      }
    }    
    ADCnum++;
    ++cMap;
  }
  
  *numADCnames=ADCnum;
  return status;
}

//-----------------------------------------------------------------------------
// Function      : xyce_getTimeStatePairsADC 
// Purpose       : Call the Xyce::Circuit::Simulator::getTimeStatePairs
//                 function via a pointer to an Xyce::Circuit::GenCouplingSimulator object. A typical 
//                 use case is to create that pointer with the xyce_open() 
//                 function, and to then initialize the Xyce::Circuit::GenCouplingSimulator object with 
//                 the xyce_initialize() function before using this function. 
//                 See the "Application Note: Mixed Signal Simulation with 
//                 Xyce" for a discussion of how to use this function. This 
//                 function returns 1 if any ADC devices were found in the
//                 simulation.  It returns 0, otherwise.
// Special Notes : This function is the least "mature" of the implemented 
//                 functions.  The Application Note has more details on that.
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical Models and Simulation
// Creation Date : 11/13/2018
//-----------------------------------------------------------------------------
int xyce_getTimeStatePairsADC( void** ptr, int * numADCnames, char ** ADCnames, int * numPoints, double ** timeArray, int ** stateArray )
{
  Xyce::Circuit::GenCouplingSimulator * xycePtr = static_cast<Xyce::Circuit::GenCouplingSimulator *>( *ptr );
 
  std::map< std::string, std::vector< std::pair<double,int> > > timeStateUpdateMap;
  int status = xycePtr->getTimeStatePairs( timeStateUpdateMap );
  
  std::map< std::string, std::vector< std::pair<double,int> > >::iterator cMap = timeStateUpdateMap.begin();
  std::map< std::string, std::vector< std::pair<double,int> > >::iterator eMap = timeStateUpdateMap.end();
  int ADCnum = 0;
  int PTnum = 0;
  while( cMap != eMap )
  {
    strcpy( ADCnames[ADCnum], (cMap->first).c_str() );
    std::vector< std::pair<double,int> > dataVec = cMap->second;
    int PTnumThisADC = 0;
    for( int j = 0; j<dataVec.size(); j++ )
    {
      timeArray[ADCnum][j] = dataVec[j].first;
      stateArray[ADCnum][j] = dataVec[j].second;
      PTnumThisADC++;
      //Xyce::dout() << "timeArray and State Array for ADC " << ADCnum << " for " << j
      //     << " are "  << dataVec[j].first << " and " << dataVec[j].second << std::endl;
    }
    if (PTnumThisADC > PTnum) PTnum = PTnumThisADC;
    ADCnum++;
    ++cMap;
  }
  *numPoints=PTnum;
  *numADCnames=ADCnum;
  
  return status;
}

//-----------------------------------------------------------------------------
// Function      : xyce_getTimeStatePairsADCLimitData 
// Purpose       : Call the Xyce::Circuit::Simulator::getTimeStatePairs
//                 function via a pointer to an Xyce::Circuit::GenCouplingSimulator object but limit the
//                 returned data assuming the maximum dimensions of the arras are:
//                 ADCnamesArray[maxNumADCnames][maxNameLength]
//                 timeArray[maxNumADCnames][maxNumPoints]
//                 stateArray[maxNumADCnames][maxNumPoints]
//                 numPointsArray[maxNumADCnames]
//
//                 numPointsArray[i] is the number of time and voltage points 
//                 for the i'th ADC in the timeArray[i][] and voltageArray[i][]
//
//                 Using numPointsArray allows the calling program to know the 
//                 number of valid points in the time and voltage arrays. 
//                 A typical use case is to create that pointer with the 
//                 xyce_open() function, and to then initialize the Xyce::Circuit::GenCouplingSimulator 
//                 object with the xyce_initialize() function before using this 
//                 function.  Also, precallocation of the arrays ADCnamesArray,
//                 timeArray, voltageArray and numPointsArray is REQUIRED.
// Scope         : public
// Creator       : Rich Schiek, SNL, Electrical Models and Simulation
// Creation Date : 12/02/2021
//-----------------------------------------------------------------------------
int xyce_getTimeStatePairsADCLimitData( void** ptr,
  const int maxNumADCnames, const int maxNameLength, const int maxNumPoints,
  int * numADCnames, char ** ADCnamesArray, int * numPointsArray, double ** timeArray, int ** stateArray )
{
  Xyce::Circuit::GenCouplingSimulator * xycePtr = static_cast<Xyce::Circuit::GenCouplingSimulator *>( *ptr );
 
  std::map< std::string, std::vector< std::pair<double,int> > > timeStateUpdateMap;
  int status = xycePtr->getTimeStatePairs( timeStateUpdateMap );
  
  std::map< std::string, std::vector< std::pair<double,int> > >::iterator cMap = timeStateUpdateMap.begin();
  std::map< std::string, std::vector< std::pair<double,int> > >::iterator eMap = timeStateUpdateMap.end();
  int ADCnum = 0;
  while( (cMap != eMap) && (ADCnum <maxNumADCnames) )
  {
    strncpy( ADCnamesArray[ADCnum], (cMap->first).c_str(), maxNameLength );
    std::vector< std::pair<double,int> > dataVec = cMap->second;
    // return the size of the array or maxNumPoints whichever is smaller.
    if(dataVec.size() < maxNumPoints)
    {
      numPointsArray[ADCnum]= dataVec.size();
      int j=0;
      while((j<dataVec.size()))
      {
        timeArray[ADCnum][j] = dataVec[j].first;
        stateArray[ADCnum][j] = dataVec[j].second;
        j++;
      }
    }
    else
    {
      // number of points is greater than allocated space.  
      // just return the LAST points as they include the end time 
      // step and steps closest to that time. 
      numPointsArray[ADCnum]= maxNumPoints;
      int j=dataVec.size()-1;
      int k= maxNumPoints-1;
      while(k>=0)
      {
        timeArray[ADCnum][k] = dataVec[j].first;
        stateArray[ADCnum][k] = dataVec[j].second;
        j--;
        k--;
      }
    }    
    ADCnum++;
    ++cMap;
  }
  
  *numADCnames=ADCnum;
  return status;
}

//-----------------------------------------------------------------------------
// Function      : xyce_setADCWidths
// Purpose       : Call the Xyce::Circuit::Simulator::setADCWidths function via 
//                 a pointer to an Xyce::Circuit::GenCouplingSimulator object.  A typical use case is to 
//                 create that pointer with the xyce_open() function, and to 
//                 then initialize the Xyce::Circuit::GenCouplingSimulator object with the xyce_initialize() 
//                 function before using this function.  This function will 
//                 return 1 if the “output bit-vector width” is successfully 
//                 updated at every ADC specified in the ADCnames parameter. 
//                 It will return 0 if the update process fails at any ADC. 
// Special Notes : Error handling for the ADCnames and widths parameters being 
//                 of unequal length is not implemented yet.
// Scope         : public
// Creator       : Rich Schiek, SNL, Electrical Models and Simulation
// Creation Date : 3/01/2018
//-----------------------------------------------------------------------------
int xyce_setADCWidths(void ** ptr, int numADCnames, char ** ADCnames, int *  widths) 
{
  Xyce::Circuit::GenCouplingSimulator * xycePtr = static_cast<Xyce::Circuit::GenCouplingSimulator *>( *ptr );
  
  std::map<std::string, int> ADCWidthMap;
  for( int i=0; i<numADCnames; i++ )
  {
    ADCWidthMap[ std::string(ADCnames[i]) ] = widths[i];
  }
  
  int status = xycePtr->setADCWidths( ADCWidthMap);
  return status;
}

//-----------------------------------------------------------------------------
// Function      : xyce_getADCWidths
// Purpose       : Call the Xyce::Circuit::Simulator::getADCWidths function via 
//                 a pointer to an Xyce::Circuit::GenCouplingSimulator object.  A typical use case is to 
//                 create that pointer with the xyce_open() function, and to 
//                 then initialize the Xyce::Circuit::GenCouplingSimulator object with the xyce_initialize() 
//                 function before using this function.  This function will 
//                 return 1 if the “output bit-vector width” is successfully 
//                 found for every ADC specified in the ADCnames parameter. 
//                 It will return 0 if the get process fails for any ADC. 
// Special Notes : The width for any YADC name not found in the netlist will be
//                 returned as 0, with an overall status then of 0.
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical Models and Simulation
// Creation Date : 11/16/2018
//-----------------------------------------------------------------------------
int xyce_getADCWidths(void ** ptr, int numADCnames, char ** ADCnames, int *  widths) 
{
  Xyce::Circuit::GenCouplingSimulator * xycePtr = static_cast<Xyce::Circuit::GenCouplingSimulator *>( *ptr );
  std::map<std::string, int> ADCWidthMap;
  for( int i=0; i<numADCnames; i++ )
  {
    ADCWidthMap[ std::string(ADCnames[i]) ] = 0;
  }
  
  int status = xycePtr->getADCWidths(ADCWidthMap);
  
  for( int i=0; i<numADCnames; i++ )
  {
    widths[i] = ADCWidthMap[ std::string(ADCnames[i]) ];
  }

  return status;
}


//-----------------------------------------------------------------------------
// Function      : xyce_getCircuitValue
// Purpose       : Call the Xyce::Circuit::Simulator::getCircuitValue
//                 function via a pointer to an Xyce::Circuit::GenCouplingSimulator object. 
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Electrical Models and Simulation
// Creation Date : 09/28/2022
//-----------------------------------------------------------------------------
double xyce_getCircuitValue( void **ptr, char * paramName)
{
  Xyce::Circuit::GenCouplingSimulator * xycePtr = static_cast<Xyce::Circuit::GenCouplingSimulator *>( *ptr );
  std::string sParamName( paramName );
  double value = 0.0;
  bool result = xycePtr->getCircuitValue(sParamName, value);
  if( ! result )
  {
    std::cerr << "xyce_getCircuitValue with paramName = " << sParamName << " returned false" << std::endl;
    value = 0.0;
  }
  return value;
}


//-----------------------------------------------------------------------------
// Function      : xyce_setCircuitParameter
// Purpose       : Call the Xyce::Circuit::Simulator::setCircuitParameter
//                 function via a pointer to an Xyce::Circuit::GenCouplingSimulator object. 
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Electrical Models and Simulation
// Creation Date : 09/28/2022
//-----------------------------------------------------------------------------
bool xyce_setCircuitParameter( void **ptr, char * paramName, double value)
{
  Xyce::Circuit::GenCouplingSimulator * xycePtr = static_cast<Xyce::Circuit::GenCouplingSimulator *>( *ptr );
  std::string sParamName( paramName );
  bool result = xycePtr->setCircuitParameter( sParamName, value);
  if( ! result )
  {
    std::cerr << "xyce_setCircuitParameter with paramName = " << sParamName << " returned false" << std::endl;
  }
  return result;
}

//-----------------------------------------------------------------------------
// Function      : xyce_checkResponseVar
// Purpose       : Call the Xyce::Circuit::Simulator::checkResponseVar
//                 function via a pointer to an Xyce::Circuit::GenCouplingSimulator object.  A typical 
//                 use case is to create that pointer with the xyce_open() 
//                 function, and to then initialize the Xyce::Circuit::GenCouplingSimulator object with 
//                 the xyce_initialize() function before using this function.
//                 This function returns 1 if variable name is a valid measure 
//                 name in the the Xyce simulation. Otherwise, it returns 0.
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Electrical Models and Simulation
// Creation Date : 3/01/2018
//-----------------------------------------------------------------------------
int xyce_checkResponseVar(void ** ptr, char * variable_name)
{
  Xyce::Circuit::GenCouplingSimulator * xycePtr = static_cast<Xyce::Circuit::GenCouplingSimulator *>( *ptr );
  std::string variable_nameStr( variable_name );
  int result = xycePtr->checkResponseVar(  variable_nameStr );
  return result;
}

//-----------------------------------------------------------------------------
// Function      : xyce_obtainResponse
// Purpose       : Call the Xyce::Circuit::Simulator::obtainResponse function via
//                 a pointer to an Xyce::Circuit::GenCouplingSimulator object.  A typical use case is to
//                 create that pointer with the xyce_open() function, and to then
//                 initialize the Xyce::Circuit::GenCouplingSimulator object with the xyce_initialize() 
//                 function before using this function. The value parameter is
//                 the value of the variable_name measure at the current sim 
//                 time. If the Xyce simulation is done then the value parameter 
//                 is that measure's value at the final sim time.  This function 
//                 returns 1 if variable_name is a valid measure name. It returns 
//                 0, otherwise. For a return value of 0, the value parameter is 
//                 also set to 0.
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Electrical Models and Simulation
// Creation Date : 3/01/2018
//-----------------------------------------------------------------------------
int xyce_obtainResponse(void ** ptr, char * variable_name, double* value) 
{
  Xyce::Circuit::GenCouplingSimulator * xycePtr = static_cast<Xyce::Circuit::GenCouplingSimulator *>( *ptr );
  std::string variable_nameStr( variable_name );
  double mVal;
  int result = xycePtr->obtainResponse(  variable_nameStr,  mVal );
  *value = mVal;
  return result;
}








