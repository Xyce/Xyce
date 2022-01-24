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
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 04/15/05
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iostream>
#include <cstring>

#include <N_CIR_Xyce.h>

#include <N_DEV_DeviceOptions.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_XyceInterface.h>
#include <N_IO_CmdParse.h>
#include <N_PDS_Comm.h>
#include <N_TIA_TwoLevelError.h>
#include <N_UTL_BreakPoint.h>
#include <N_UTL_FeatureTest.h>

namespace Xyce {
namespace Device {

namespace {
inline char *strdup(const char *s) 
{
  const int length = strlen(s);

  return std::copy(s, s + length, new char[length + 1]);
}
}

//-----------------------------------------------------------------------------
// Function      : XyceInterface::XyceInterface
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/15/05
//-----------------------------------------------------------------------------
XyceInterface::XyceInterface(
  const DeviceOptions & device_options,
  const IO::CmdParse &  command_line,
  const std::string &   netlist_filename)
  : netlistFilename_(netlist_filename),
    simulator_(0),
    commandLine_(command_line)
{}

//-----------------------------------------------------------------------------
// Function      : XyceInterface::~XyceInterface
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/15/05
//-----------------------------------------------------------------------------
XyceInterface::~XyceInterface()

{
  if (simulator_ != NULL)
  {
    simulator_->finishSolvers();
    simulator_->finalize();
  }
  delete simulator_;
}

//-----------------------------------------------------------------------------
// Function      : XyceInterface::initialize
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/01/2006
//-----------------------------------------------------------------------------
bool XyceInterface::initialize(
  Parallel::Communicator*          comm)
{
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "In XyceInterface::initialize" << std::endl;
  }

  simulator_ = new Circuit::SecondLevelSimulator(comm->comm());

  if (Parallel::rank(comm->comm()) == 0)
  {
    // Reset the netlist name in the local copy of the command line arguments:
    commandLine_.setNetlist(netlistFilename_);
  }

  simulator_->initialize(commandLine_.argc(), commandLine_.argv());

  simulator_->startupSolvers();

  return true;
}

//-----------------------------------------------------------------------------
// Function      : XyceInterface::simulateStep
// Purpose       :
// Special Notes : The 'input vector' is a vector of voltages corresponding
//                 to the connected nodes.
//
//                 The 'output vector', is a vector of currents corresponding
//                 to the same nodes.
//
//                 The jacobian is mainly conductances - dI/dV.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/01/2006
//-----------------------------------------------------------------------------
bool
XyceInterface::simulateStep(
  const SolverState &                   solState,
  const std::map<std::string,double> &  inputMap,
  std::vector<double> &                 outputVector,
  std::vector< std::vector<double> > &  jacobian,
  TimeIntg::TwoLevelError &             tlError)
{
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "In XyceInterface::simulateStep" << std::endl;
  }

  return simulator_->simulateStep(
    solState.initJctFlag_,
    inputMap,
    outputVector,
    jacobian,
    tlError);
}

//-----------------------------------------------------------------------------
// Function      : XyceInterface::finalize
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/01/2006
//-----------------------------------------------------------------------------
bool XyceInterface::finalize ()
{
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "In XyceInterface::finalize" << std::endl;
  }


  return true;
}

//-----------------------------------------------------------------------------
// Function      : XyceInterface::run
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/15/05
//-----------------------------------------------------------------------------
bool XyceInterface::run()
{
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "In XyceInterface::run" << std::endl;
  }

  int argc = 2;
  char *argv[3];
  argv[0] = strdup("Xyce");
  argv[1] = strdup(netlistFilename_.c_str());
  argv[2] = 0;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "argv[0] = " << std::string(argv[0]) << std::endl;
  }

  simulator_->run(argc, argv);

  delete[] argv[0];
  delete[] argv[1];

  return true;
}

//-----------------------------------------------------------------------------
// Function      : XyceInterface::homotopyStepSuccess
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
void
XyceInterface::homotopyStepSuccess(
  const std::vector<std::string> &      paramNames,
  const std::vector<double> &           paramVals)
{
  simulator_->homotopyStepSuccess(paramNames, paramVals);
}


//-----------------------------------------------------------------------------
// Function      : XyceInterface::homotopyStepFailure
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/30/06
//-----------------------------------------------------------------------------
void XyceInterface::homotopyStepFailure()
{
  simulator_->homotopyStepFailure();
}

//-----------------------------------------------------------------------------
// Function      : XyceInterface::stepSuccess
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/12/06
//-----------------------------------------------------------------------------
void XyceInterface::stepSuccess(Analysis::TwoLevelMode analysis)
{
  simulator_->stepSuccess(analysis);
}

//-----------------------------------------------------------------------------
// Function      : XyceInterface::stepFailure
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/12/06
//-----------------------------------------------------------------------------
void XyceInterface::stepFailure(Analysis::TwoLevelMode analysis)
{
  simulator_->stepFailure(analysis);
}

//-----------------------------------------------------------------------------
// Function      : XyceInterface::getInitialQnorm
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/15/06
//-----------------------------------------------------------------------------
bool XyceInterface::getInitialQnorm (TimeIntg::TwoLevelError & tle)
{
  bool bsuccess = true;

  if (simulator_)
  {
    bsuccess = simulator_->getInitialQnorm (tle);
  }
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : XyceInterface::getBreakPoints
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/14/06
//-----------------------------------------------------------------------------
bool XyceInterface::getBreakPoints
    (std::vector<Util::BreakPoint> &breakPointTimes,
     std::vector<Util::BreakPoint> &pauseBreakPointTimes)
{
  bool bsuccess = true;

  if (simulator_)
  {
    bsuccess = simulator_->getBreakPoints(
        breakPointTimes,
        pauseBreakPointTimes
        );
  }
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : XyceInterface::updateStateArrays
// Purpose       : ERK:  This doesn't do anything....
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/18/06
//-----------------------------------------------------------------------------
bool XyceInterface::updateStateArrays ()
{
  bool bsuccess = true;
  if (simulator_)
  {
    bsuccess = simulator_->updateStateArrays ();
  }
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : XyceInterface::startTimeStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
bool
XyceInterface::startTimeStep(
  bool          beginIntegrationFlag,
  double        nextTimeStep,
  double        nextTime,
  int           currentOrder)
{
  bool bsuccess = true;
  if (simulator_)
  {
    bsuccess = simulator_->startTimeStep(
      beginIntegrationFlag,
      nextTimeStep,
      nextTime,
      currentOrder);
  }
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : XyceInterface::setInternalParam
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/18/06
//-----------------------------------------------------------------------------
bool XyceInterface::setInternalParam (const std::string & name, double val)
{
  bool bsuccess = true;
  if (simulator_)
  {
    bsuccess = simulator_->setInternalParam (name, val);
  }
  return bsuccess;
}

} // namespace Device
} // namespace Xyce
