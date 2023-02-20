//-------------------------------------------------------------------------
//   Copyright 2002-2023 National Technology & Engineering Solutions of
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
// Creation Date  : 03/13/06
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------
#include <iostream>

#include "Teuchos_ParameterList.hpp"

// ----------   Xyce Includes   ----------
#include <N_DEV_CharonInterface.h>
#include <N_CIR_Xyce.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_SolverState.h>
#include <N_ERH_ErrorMgr.h>
#include <N_UTL_BreakPoint.h>
#include <N_UTL_FeatureTest.h>

#include <N_TIA_TwoLevelError.h>

// RPP: We have added a copy of a charon interface file to xyce to
// eliminate a circular dependency in nevada's build TPL process.
// This file must be kept up-to-date with Charon's otherwise they will
// not be able to link with each other.
#ifdef Xyce_CHARON
#include "Charon_CircuitInterface.h"
#endif

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Function      : CharonInterface::CharonInterface
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/13/06
//-----------------------------------------------------------------------------
CharonInterface::CharonInterface(
  const DeviceOptions & do1,
  const std::string &   netlist,
  const SolverState &   ss1)
  : inputFileName_(netlist),
    devOptions_(do1),
    solState_(ss1)
{

}

//-----------------------------------------------------------------------------
// Function      : CharonInterface::CharonInterface
// Purpose       : copy constructor
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/13/06
//-----------------------------------------------------------------------------
CharonInterface::CharonInterface (const CharonInterface &right)
  : inputFileName_(right.inputFileName_),
    devOptions_(right.devOptions_),
    solState_(right.solState_)
{

}

//-----------------------------------------------------------------------------
// Function      : CharonInterface::~CharonInterface
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/13/06
//-----------------------------------------------------------------------------
CharonInterface::~CharonInterface()

{

}

//-----------------------------------------------------------------------------
// Function      : CharonInterface::initialize
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/13/2006
//-----------------------------------------------------------------------------
bool CharonInterface::initialize(Parallel::Communicator* comm)
{

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "In CharonInterface::initialize" << std::endl;
  }

  input_list_ = Teuchos::rcp(new Teuchos::ParameterList);

  output_list_ = Teuchos::rcp(new Teuchos::ParameterList);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : CharonInterface::simulateStep
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
// Creation Date : 03/13/2006
//-----------------------------------------------------------------------------
bool
CharonInterface::simulateStep(
  const SolverState &                   solState,
  const std::map<std::string,double> &  inputMap,
  std::vector<double> &                 outputVector,
  std::vector< std::vector<double> > &  jacobian,
  TimeIntg::TwoLevelError &             tlError)
{

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "In CharonInterface::simulateStep" << std::endl;
  }

  // RPP: need to get this from xyce time integrator...
  double currentTime = solState.currTime_;
  double stepSize = solState.currTimeStep_;

  input_list_->set("Current Time", currentTime);
  input_list_->set("Time Step Size", stepSize);
  input_list_->set("Time Step Type", "Backward Euler");

  // Tell charon whether we are in a steady state or transient mode
  if (solState_.dcopFlag)
    input_list_->set("Solve Type", "Steady State");
  else
    input_list_->set("Solve Type", "Transient");

#ifdef Xyce_CHARON

  charon::sc::CircuitInterface::getInstance().takeStep(input_list_,
						       inputMap,
						       output_list_,
						       outputVector,
						       jacobian);

#else
  Report::DevelFatal().in("CharonInterface::simulateStep") << "Charon support has not been enabled.  Rebuild xyce with the flag --enable-charon";
  return true;
#endif

  // Get the output parameters
  if (output_list_->isParameter("Error Sum"))
    tlError.xErrorSum = output_list_->get<double>("Error Sum");
  else
    tlError.xErrorSum = 0.0;

  if (output_list_->isParameter("Inner Size"))
    tlError.innerSize = output_list_->get<double>("Number of DOF");
  else
    tlError.innerSize = 1.0;

  // Get the return status
  int return_status = -3;
  if (output_list_->isParameter("Charon Status"))
    return_status = output_list_->get<int>("Charon Status");
  if (return_status < 0)
    return false;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : CharonInterface::finalize
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/13/2006
//-----------------------------------------------------------------------------
bool CharonInterface::finalize ()
{
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "In CharonInterface::finalize" << std::endl;
  }


  return true;
}

//-----------------------------------------------------------------------------
// Function      : CharonInterface::run
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/13/06
//-----------------------------------------------------------------------------
bool CharonInterface::run ()
{
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "In CharonInterface::run" << std::endl;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : CharonInterface::homotopyStepSuccess
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
void CharonInterface::homotopyStepSuccess
      (const std::vector<std::string> & paramNames,
       const std::vector<double> & paramVals)
{

  return;
}

//-----------------------------------------------------------------------------
// Function      : CharonInterface::homotopyStepFailure
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/30/06
//-----------------------------------------------------------------------------
void CharonInterface::homotopyStepFailure ()
{

  return;
}

//-----------------------------------------------------------------------------
// Function      : CharonInterface::stepSuccess
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/13/06
//-----------------------------------------------------------------------------
void CharonInterface::stepSuccess(Analysis::TwoLevelMode analysis)
{

#ifdef Xyce_CHARON
  bool is_active = true;
  charon::sc::CircuitInterface::getInstance().acceptTimeStep(is_active);

#else
  Report::DevelFatal().in("CharonInterface::stepSuccess") << "Charon support has not been enabled.  Rebuild xyce with the flag --enable-charon";
#endif
}

//-----------------------------------------------------------------------------
// Function      : CharonInterface::stepFailure
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/13/06
//-----------------------------------------------------------------------------
void CharonInterface::stepFailure(Analysis::TwoLevelMode analysis)
{

  return;
}


//-----------------------------------------------------------------------------
// Function      : CharonInterface::getInitialQnorm
// Purpose       :
// Special Notes : no-op.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/15/06
//-----------------------------------------------------------------------------
bool CharonInterface::getInitialQnorm (TimeIntg::TwoLevelError & tle)
{
  return true;
}

} // namespace Device
} // namespace Xyce
