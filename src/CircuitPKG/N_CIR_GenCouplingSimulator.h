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
// Purpose        : Provide a class for more general Xyce/Alegra coupling
//
// Special Notes  : This class is meant to provide a more general Xyce/Alegra
//                  coupling API than that provided by the old N_CIR_Xygra class,
//                  which is fairly general but clearly set up and named
//                  for its primary use case, simulating coils in Alegra
//                  and coupling them loosely to Xyce.
//
// Creator        : Tom Russo, SNL, Electrical Models & Simulation
//
// Creation Date  : 2/22/2017
//
//-------------------------------------------------------------------------

#ifndef Xyce_N_CIR_GenCouplingSimulator_H
#define Xyce_N_CIR_GenCouplingSimulator_H

#include <map>
#include <string>

#include <N_CIR_Xyce.h>
#include <N_PDS_fwd.h>
#include <N_DEV_fwd.h>
#include <N_IO_fwd.h>
#include <N_DEV_VectorComputeInterface.h>
#include <N_IO_ExtOutInterface.h>

namespace Xyce {
namespace Circuit {

//-----------------------------------------------------------------------------
// Class         : GenCouplingSimulator
// Purpose       :
// Special Notes :
//
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 8/21/08
//-----------------------------------------------------------------------------
///
/// High-level Xyce interface class for use in coupling to external codes
///
/// This class is derived from the "Simulator" class of Xyce.  It provides
/// some extra functions that allow the external simulator to pass data
/// down to special interface devices
///
/// This is very similar to N_CIR_Xygra, but intended to be more general
class GenCouplingSimulator : public Xyce::Circuit::Simulator
{
public:
  /// Constructor
  GenCouplingSimulator(Xyce::Parallel::Machine comm=MPI_COMM_NULL)
    : Simulator(comm)
  {}

  /// Destructor
  virtual ~GenCouplingSimulator()
  {}


  // API functions
  /// Sets the number of internal variables of a named device instance
  /// this *MUST* be called in between initializeEarly and intializeLate
  bool setNumInternalVars(const std::string & deviceName,const int numInt);

  /// Sets the number of store variables of a named device instance
  /// this *MUST* be called in between initializeEarly and intializeLate
  bool setNumStoreVars(const std::string & deviceName,const int numStore);

  /// Sets the number of state variables of a named device instance
  /// this *MUST* be called in between initializeEarly and intializeLate
  bool setNumStateVars(const std::string & deviceName,const int numState);

  /// Sets the number of branch data variables of a named device instance
  bool setNumBranchDataVars(const std::string & deviceName,const int numBranchData);

  /// Sets the number of branch data variables if allocated of a named device instance
  bool setNumBranchDataVarsIfAllocated(const std::string & deviceName,
                                       const int numBranchDataIfAllocated);

  /// Returns the number of variables (external nodes + internal variables)
  /// of a named "GeneralExternal" device instance.
  int getNumVars(const std::string & deviceName);

  /// Returns the number of EXTERNAL variables
  /// of a named "GeneralExternal" device instance.
  int getNumExtVars(const std::string & deviceName);

  /// Populates a vector with the current values of the solution vector
  /// for the variables associated with a named device.
  bool getSolution(const std::string & deviceName, std::vector<double> & sV);

  /// Associate a vector loader object pointer with the named device
  bool setJacStamp(const std::string & deviceName, std::vector< std::vector<int> > &jS);

  /// Associate a vector loader object pointer with the named device
  bool setVectorLoader(const std::string & deviceName, Xyce::Device::VectorComputeInterface * vciPtr);

  /// Populates a vector with the current values of the name/value pairs
  /// for the double params associated with a named device.
  bool getDParams(const std::string & deviceName,
                  std::vector<std::string> &pNames,
                  std::vector<double> & pValues);
  /// Populates a vector with the current values of the name/value pairs
  /// for the int params associated with a named device.
  bool getIParams(const std::string & deviceName,
                  std::vector<std::string> &pNames,
                  std::vector<int> & pValues);
  /// Populates a vector with the current values of the name/value pairs
  /// for the bool params associated with a named device.
  bool getBParams(const std::string & deviceName,
                  std::vector<std::string> &pNames,
                  std::vector<bool> & pValues);
  /// Populates a vector with the current values of the name/value pairs
  /// for the string params associated with a named device.
  bool getSParams(const std::string & deviceName,
                  std::vector<std::string> &pNames,
                  std::vector<std::string> & pValues);

  /// Retrieves Netlist file path
  std::string getNetlistFilePath() const;

  /// Retrieves Xyce executable file path
  std::string getXyceFilePath() const;

  bool addOutputInterface(Xyce::IO::ExternalOutputInterface * extIntPtr);

private:
  std::map<std::string,Xyce::Device::GeneralExternal::Instance *> genExtDevMap_;  //< Mapping of device names to pointers for fast lookups

  /// Retreive a pointer to a device instance by name
  Xyce::Device::GeneralExternal::Instance *getGeneralExternalDeviceInstance_(const std::string & deviceName);


};

} // namespace Circuit
} // namespace Xyce

#endif // Xyce_N_CIR_GenCouplingSimulator_H
