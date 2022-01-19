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
// Purpose        : Provide a class for loose coupling with Alegra
//                  or other codes
//
// Special Notes  : This interface allows an external code to control
//                  the time stepping, calling Xyce to run coupled
//                  circuit problems on very small time scales.
//
//                  The external simulator is responsible for intializing
//                  Xyce via the API calls, passing critical data to the
//                  "GeneralExteral" devices, and asking Xyce to simulate
//                  for a requested time interval.  The external code
//                  can then query the devices for their final solution
//                  variables and repeat the process.
//
// Creator        : Tom Russo, SNL, Electrical Models & Simulation
//
// Creation Date  : 2/27/2017
//
//-------------------------------------------------------------------------
#include <Xyce_config.h>

#include <N_CIR_GenCouplingSimulator.h>

#include <N_DEV_Device.h>
#include <N_DEV_DeviceMgr.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_ExtOutInterface.h>
#include <N_DEV_Algorithm.h>
#include <N_DEV_GeneralExternal.h>
#include <N_DEV_VectorComputeInterface.h>

namespace Xyce {
namespace Circuit {


//-----------------------------------------------------------------------------
// Function      : GenCouplingSimulator::getGeneralExternalDeviceInstance
// Purpose       : Returns the pointer to a named general coupling instance
// Special Notes :
// Scope         : private
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 2/27/2017
//-----------------------------------------------------------------------------
///
/// This method is intended to be called internally
/// by other methods so that they can retrieve a pointer
/// to a device instance by name and access it directly
/// without having to "punch through" layers of device
/// manager/loader functions.
/// We maintain a map of device names to pointers here,
/// so that we don't have to do repeated device manager
/// lookups every time we operate on a device.
///
/// @param deviceName The name of the device to retrieve
/// @return pointer to the device instance object
///
/// @author Tom Russo, SNL, Electrical Models and Simulation
/// @date 2/27/2017
Xyce::Device::GeneralExternal::Instance *
GenCouplingSimulator::getGeneralExternalDeviceInstance_(const std::string & deviceName)
{
  // See if we've looked this up before.
  if (genExtDevMap_.empty())
  {
    Xyce::Device::Device *device = getDeviceManager().getDevice(Xyce::Device::GeneralExternal::Traits::modelGroup());
    if (device)
      Xyce::Device::mapDeviceInstances(*device, genExtDevMap_);
  }

  std::map<std::string, Xyce::Device::GeneralExternal::Instance *>::iterator mapIter = genExtDevMap_.find(deviceName);
  if (mapIter == genExtDevMap_.end())
    return 0;

  return (*mapIter).second;
}

//-----------------------------------------------------------------------------
// Function      : GenCouplingSimulator::setNumInternalVars
// Purpose       : Set the number of internal vars for named instance
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 2/27/2017
//-----------------------------------------------------------------------------
///
/// Given the name of a GeneralExternal device, set the number of internal
/// variables the device should have.
///
/// Must be called in between intializeEarly and intializeLate, since
/// this information is used during initializeLate to perform final
/// problem initialization.
///
/// @param[in] deviceName   The name of the device as set in the netlist
/// @param[in] numInt       Number of internal variables
/// @return  true if device of that name found, false if not found

bool GenCouplingSimulator::setNumInternalVars(const std::string & deviceName, const int numInt)
{
  bool success=true;
  Xyce::Device::GeneralExternal::Instance * genExtPtr = getGeneralExternalDeviceInstance_(deviceName);
  if (genExtPtr)
    genExtPtr -> setNumInternalVars(numInt);
  else
    success=false;

  return success;
}

//-----------------------------------------------------------------------------
// Function      : GenCouplingSimulator::setNumStoreVars
// Purpose       : Set the number of store vars for named instance
// Special Notes :
// Scope         : public
// Creator       : Paul Kuberry, SNL
// Creation Date : 12/10/2020
//-----------------------------------------------------------------------------
///
/// Given the name of a GeneralExternal device, set the number of store
/// variables the device should have.
///
/// Must be called in between intializeEarly and intializeLate, since
/// this information is used during initializeLate to perform final
/// problem initialization.
///
/// @param[in] deviceName   The name of the device as set in the netlist
/// @param[in] numStore       Number of store variables
/// @return  true if device of that name found, false if not found

bool GenCouplingSimulator::setNumStoreVars(const std::string & deviceName, const int numStore)
{
  bool success=true;
  Xyce::Device::GeneralExternal::Instance * genExtPtr = getGeneralExternalDeviceInstance_(deviceName);
  if (genExtPtr)
    genExtPtr -> setNumStoreVars(numStore);
  else
    success=false;

  return success;
}

//-----------------------------------------------------------------------------
// Function      : GenCouplingSimulator::setNumStateVars
// Purpose       : Set the number of state vars for named instance
// Special Notes :
// Scope         : public
// Creator       : Paul Kuberry, SNL
// Creation Date : 12/10/2020
//-----------------------------------------------------------------------------
///
/// Given the name of a GeneralExternal device, set the number of state
/// variables the device should have.
///
/// Must be called in between intializeEarly and intializeLate, since
/// this information is used during initializeLate to perform final
/// problem initialization.
///
/// @param[in] deviceName   The name of the device as set in the netlist
/// @param[in] numState       Number of state variables
/// @return  true if device of that name found, false if not found

bool GenCouplingSimulator::setNumStateVars(const std::string & deviceName, const int numState)
{
  bool success=true;
  Xyce::Device::GeneralExternal::Instance * genExtPtr = getGeneralExternalDeviceInstance_(deviceName);
  if (genExtPtr)
    genExtPtr -> setNumStateVars(numState);
  else
    success=false;

  return success;
}

//-----------------------------------------------------------------------------
// Function      : GenCouplingSimulator::setNumBranchDataVars
// Purpose       : Set the number of branch data vars for named instance
// Special Notes :
// Scope         : public
// Creator       : Paul Kuberry, SNL
// Creation Date : 12/10/2020
//-----------------------------------------------------------------------------
///
/// Given the name of a GeneralExternal device, set the number of branch data
/// variables the device should have.
///
/// Must be called in between intializeEarly and intializeLate, since
/// this information is used during initializeLate to perform final
/// problem initialization.
///
/// @param[in] deviceName   The name of the device as set in the netlist
/// @param[in] numBranchData       Number of branch data variables
/// @return  true if device of that name found, false if not found

bool GenCouplingSimulator::setNumBranchDataVars(const std::string & deviceName, const int numBranchData)
{
  bool success=true;
  Xyce::Device::GeneralExternal::Instance * genExtPtr = getGeneralExternalDeviceInstance_(deviceName);
  
  if (genExtPtr)
    genExtPtr -> setNumBranchDataVars(numBranchData);
  else
    success=false;

  return success;
}

//-----------------------------------------------------------------------------
// Function      : GenCouplingSimulator::setNumBranchDataVarsIfAllocated
// Purpose       : Set the number of branch data if allocated vars for named instance
// Special Notes :
// Scope         : public
// Creator       : Paul Kuberry, SNL
// Creation Date : 12/10/2020
//-----------------------------------------------------------------------------
///
/// Given the name of a GeneralExternal device, set the number of branch data if allocated
/// variables the device should have.
///
/// Must be called in between intializeEarly and intializeLate, since
/// this information is used during initializeLate to perform final
/// problem initialization.
///
/// @param[in] deviceName   The name of the device as set in the netlist
/// @param[in] numBranchDataIfAllocated       Number of branch data if allocated variables
/// @return  true if device of that name found, false if not found

bool GenCouplingSimulator::setNumBranchDataVarsIfAllocated(const std::string & deviceName, const int numBranchDataIfAllocated)
{
  bool success=true;
  Xyce::Device::GeneralExternal::Instance * genExtPtr = getGeneralExternalDeviceInstance_(deviceName);
  if (genExtPtr)
    genExtPtr -> setNumBranchDataVarsIfAllocated(numBranchDataIfAllocated);
  else
    success=false;

  return success;
}

//-----------------------------------------------------------------------------
// Function      : GenCouplingSimulator::getNumVars
// Purpose       : Query number variables used by named device
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 2/27/2017
//-----------------------------------------------------------------------------
///
/// Given the name of a GeneralExternal device, return the number of
/// variables (internal+external) the device should have.
///
/// @param[in] deviceName   The name of the device as set in the netlist
/// @return  Number of variables (will return -1 if device not found)

int GenCouplingSimulator::getNumVars(const std::string & deviceName)
{
  Xyce::Device::GeneralExternal::Instance * genExtPtr = getGeneralExternalDeviceInstance_(deviceName);
  if (genExtPtr)
    return genExtPtr -> getNumVars();
  return -1;
}

//-----------------------------------------------------------------------------
// Function      : GenCouplingSimulator::getNumExtVars
// Purpose       : Query number external variables used by named device
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 2/27/2017
//-----------------------------------------------------------------------------
///
/// Given the name of a GeneralExternal device, return the number of
/// external variables the device should have.
///
/// @param[in] deviceName   The name of the device as set in the netlist
/// @return  Number of variables (will return -1 if device not found)

int GenCouplingSimulator::getNumExtVars(const std::string & deviceName)
{
  Xyce::Device::GeneralExternal::Instance * genExtPtr = getGeneralExternalDeviceInstance_(deviceName);
  if (genExtPtr)
    return genExtPtr -> getNumExtVars();
  return -1;
}

//-----------------------------------------------------------------------------
// Function      : GenCouplingSimulator::getSolution
// Purpose       : Return the solutionv values associated with named device
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 2/27/2017
//-----------------------------------------------------------------------------
///
/// Given the name of a GeneralExternal device, populate the given vector
/// with the solution variables associated with that device
///
/// @param[in] deviceName   The name of the device as set in the netlist
/// @param[out] sV  std::vector<double> reference in which to place the values
/// @return  True if device found, false if not

bool GenCouplingSimulator::getSolution(const std::string & deviceName,
                                       std::vector<double> & sV)
{
  bool success=true;
  Xyce::Device::GeneralExternal::Instance * genExtPtr = getGeneralExternalDeviceInstance_(deviceName);
  if (genExtPtr)
    genExtPtr -> getSolution(sV);
  else
    success=false;
  return success;
}

//-----------------------------------------------------------------------------
// Function      : GenCouplingSimulator::setJacStamp
// Purpose       : Set the Jacobian stamp for a GenExt device
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 2/27/2017
//-----------------------------------------------------------------------------
///
/// Given the name of a GeneralExternal device, copy the given
/// jacStamp into the object.  See the Xyce::Device::GeneralExternal device
/// for details about the jacStamp structure.
///
/// @param[in] deviceName   The name of the device as set in the netlist
/// @param[in] jS Reference to a jacStamp 
/// @return  True if device found, false if not

bool GenCouplingSimulator::setJacStamp(const std::string & deviceName,
                                       std::vector< std::vector<int> > &jS)
{
  bool success=true;
  Xyce::Device::GeneralExternal::Instance * genExtPtr = getGeneralExternalDeviceInstance_(deviceName);
  if (genExtPtr)
    genExtPtr -> setJacStamp(jS);
  else
    success=false;
  return success;
}

//-----------------------------------------------------------------------------
// Function      : GenCouplingSimulator::setVectorLoader
// Purpose       : Associate a vector loader object with named device
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 2/27/2017
//-----------------------------------------------------------------------------
///
/// Given the name of a GeneralExternal device, associate the pointer to
/// an external-simulator-provide vector loader object
///
/// @param[in] deviceName   The name of the device as set in the netlist
/// @param[in] vciPtr   pointer to object that implements the VectorComputeInterface
/// @return  True if device found, false if not

bool GenCouplingSimulator::setVectorLoader(const std::string & deviceName,
                                           Xyce::Device::VectorComputeInterface * vciPtr)
{
  bool success=true;
  Xyce::Device::GeneralExternal::Instance * genExtPtr = getGeneralExternalDeviceInstance_(deviceName);
  if (genExtPtr)
    success=genExtPtr -> setVectorLoader(vciPtr);
  else
    success=false;
  return success;
}

//-----------------------------------------------------------------------------
// Function      : GenCouplingSimulator::getDParams
// Purpose       : Return the name/value pairs associated with named device's
//                 DPARAMS vector-composite
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 11/29/2017
//-----------------------------------------------------------------------------
///
/// Given the name of a GeneralExternal device, populate the given vectors
/// with the name/value pairs for double parameters associated with that device
///
/// @param[in] deviceName   The name of the device as set in the netlist
/// @param[out] pNames  std::vector<std::string> reference in which to place the names
/// @param[out] pValues  std::vector<double> reference in which to place the values
/// @return  True if device found, false if not

bool GenCouplingSimulator::getDParams(const std::string & deviceName,
                                      std::vector<std::string> &pNames,
                                      std::vector<double> & pValues)
{
  bool success=true;
  Xyce::Device::GeneralExternal::Instance * genExtPtr = getGeneralExternalDeviceInstance_(deviceName);
  if (genExtPtr)
    genExtPtr -> getDParams(pNames,pValues);
  else
    success=false;
  return success;
}

//-----------------------------------------------------------------------------
// Function      : GenCouplingSimulator::getIParams
// Purpose       : Return the name/value pairs associated with named device's
//                 IPARAMS vector-composite
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 11/29/2017
//-----------------------------------------------------------------------------
///
/// Given the name of a GeneralExternal device, populate the given vectors
/// with the name/value pairs for int parameters associated with that device
///
/// @param[in] deviceName   The name of the device as set in the netlist
/// @param[out] pNames  std::vector<std::string> reference in which to place the names
/// @param[out] pValues  std::vector<int> reference in which to place the values
/// @return  True if device found, false if not

bool GenCouplingSimulator::getIParams(const std::string & deviceName,
                                      std::vector<std::string> &pNames,
                                      std::vector<int> & pValues)
{
  bool success=true;
  Xyce::Device::GeneralExternal::Instance * genExtPtr = getGeneralExternalDeviceInstance_(deviceName);
  if (genExtPtr)
    genExtPtr -> getIParams(pNames,pValues);
  else
    success=false;
  return success;
}

//-----------------------------------------------------------------------------
// Function      : GenCouplingSimulator::getBParams
// Purpose       : Return the name/value pairs associated with named device's
//                 BPARAMS vector-composite
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 11/29/2017
//-----------------------------------------------------------------------------
///
/// Given the name of a GeneralExternal device, populate the given vectors
/// with the name/value pairs for bool parameters associated with that device
///
/// @param[in] deviceName   The name of the device as set in the netlist
/// @param[out] pNames  std::vector<std::string> reference in which to place the names
/// @param[out] pValues  std::vector<bool> reference in which to place the values
/// @return  True if device found, false if not

bool GenCouplingSimulator::getBParams(const std::string & deviceName,
                                      std::vector<std::string> &pNames,
                                      std::vector<bool> & pValues)
{
  bool success=true;
  Xyce::Device::GeneralExternal::Instance * genExtPtr = getGeneralExternalDeviceInstance_(deviceName);
  if (genExtPtr)
    genExtPtr -> getBParams(pNames,pValues);
  else
    success=false;
  return success;
}

//-----------------------------------------------------------------------------
// Function      : GenCouplingSimulator::getSParams
// Purpose       : Return the name/value pairs associated with named device's
//                 SPARAMS vector-composite
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 11/29/2017
//-----------------------------------------------------------------------------
///
/// Given the name of a GeneralExternal device, populate the given vectors
/// with the name/value pairs for string parameters associated with that device
///
/// @param[in] deviceName   The name of the device as set in the netlist
/// @param[out] pNames  std::vector<std::string> reference in which to place the names
/// @param[out] pValues  std::vector<std::string> reference in which to place the values
/// @return  True if device found, false if not

bool GenCouplingSimulator::getSParams(const std::string & deviceName,
                                      std::vector<std::string> &pNames,
                                      std::vector<std::string> & pValues)
{
  bool success=true;
  Xyce::Device::GeneralExternal::Instance * genExtPtr = getGeneralExternalDeviceInstance_(deviceName);
  if (genExtPtr)
    genExtPtr -> getSParams(pNames,pValues);
  else
    success=false;
  return success;
}

std::string GenCouplingSimulator::getNetlistFilePath() const
{
  return commandLine_.getArgumentValue("netlist");
}

std::string GenCouplingSimulator::getXyceFilePath() const
{
  return commandLine_.argv()[0];
}

bool GenCouplingSimulator::addOutputInterface(Xyce::IO::ExternalOutputInterface * extIntPtr)
{
  
  if (extIntPtr)
  {
    getOutputManager().addExternalOutputInterface(extIntPtr);
    return true;
  }
  else
  {
    return false;
  }
}

} // namespace Circuit
} // namespace Xyce
