//-------------------------------------------------------------------------
//   Copyright 2002-2021 National Technology & Engineering Solutions of
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
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL
//
// Creation Date  : 10/xx/2019
//
//
//
//
//-----------------------------------------------------------------------------

#include <iostream>
#include <unordered_map>
#include <string>
#include <random>

#include <paramParsingExpressionGroup.h>
#include <ast.h>
#include <newExpression.h>

#include <N_IO_CmdParse.h>

#include <N_TIA_DataStore.h>
#include <N_TIA_StepErrorControl.h>
#include <N_TOP_Topology.h>
#include <N_LAS_Vector.h>
#include <N_PDS_Comm.h>
#include <N_PDS_MPI.h>
#include <N_PDS_Manager.h>
#include <N_PDS_Serial.h>

#include <N_ANP_AnalysisManager.h>
#include <N_ANP_UQSupport.h>
#include <N_ANP_SweepParam.h>
#include <N_DEV_DeviceMgr.h>
#include <N_DEV_Op.h>
#include <N_DEV_Const.h>

#include <N_ERH_ErrorMgr.h>

#include <N_UTL_DeviceNameConverters.h>

namespace Xyce {
namespace Util {

//-------------------------------------------------------------------------------
// Function      : paramParsingExpressionGroup::paramParsingExpressionGroup 
// Purpose       : constructor
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 
//-------------------------------------------------------------------------------
paramParsingExpressionGroup::paramParsingExpressionGroup ( 
 Parallel::Communicator & comm, Topo::Topology & top,
 Analysis::AnalysisManager &analysis_manager,
 Device::DeviceMgr & device_manager,
 IO::OutputMgr &output_manager
 ) :
 comm_(comm),
 top_(top),
 analysisManager_(analysis_manager),
 deviceManager_(device_manager),
 outputManager_(output_manager),
 time_(0.0), temp_(0.0), VT_(0.0), freq_(0.0), gmin_(0.0), dt_(0.0), alpha_(0.0)
{
}

//-------------------------------------------------------------------------------
// Function      : paramParsingExpressionGroup::paramParsingExpressionGroup 
// Purpose       : destructor
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 
//-------------------------------------------------------------------------------
paramParsingExpressionGroup::~paramParsingExpressionGroup ()
{
}

//-------------------------------------------------------------------------------
// Function      : paramParsingExpressionGroup::getTimeStep
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/20/2020 
//-------------------------------------------------------------------------------
double paramParsingExpressionGroup::getTimeStep ()
{
  dt_ = deviceManager_.getSolverState().currTimeStep_;
  return dt_;
}

//-------------------------------------------------------------------------------
// Function      : paramParsingExpressionGroup::getTime
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/20/2020 
//-------------------------------------------------------------------------------
double paramParsingExpressionGroup::getTime() 
{ 
  // I would have preferred to use this but as the code is currently written it
  // is not safe.  The earliest call I would need to make to getTime happens before 
  // the stepErrorControl class has been allocated.  Unfortunately, the analysis
  // manager accessor returns an invalid reference in that case, which I can't 
  // really test for.
  //
  //const TimeIntg::StepErrorControl & secControl_ = (analysisManager_.getStepErrorControl());
  //time_ = secControl_.nextTime;
  
  time_ = deviceManager_.getSolverState().currTime_;
  return time_;
} 

//-------------------------------------------------------------------------------
// Function      : paramParsingExpressionGroup::getTemp
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/20/2020 
//-------------------------------------------------------------------------------
double paramParsingExpressionGroup::getTemp() 
{ 
  temp_ = deviceManager_.getDeviceOptions().temp.getImmutableValue<double>() - CONSTCtoK;
  return temp_;
} 

//-------------------------------------------------------------------------------
// Function      : paramParsingExpressionGroup::
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/20/2020 
//-------------------------------------------------------------------------------
double paramParsingExpressionGroup::getVT  () 
{ 
  VT_ = (deviceManager_.getDeviceOptions().temp.getImmutableValue<double>())*CONSTKoverQ;
  return VT_;
} 

//-------------------------------------------------------------------------------
// Function      : paramParsingExpressionGroup::
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/20/2020 
//-------------------------------------------------------------------------------
double paramParsingExpressionGroup::getFreq() 
{ 
  freq_ = deviceManager_.getSolverState().currFreq_;
  return freq_;
} 

//-------------------------------------------------------------------------------
// Function      : paramParsingExpressionGroup::getGmin
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/20/2020 
//-------------------------------------------------------------------------------
double paramParsingExpressionGroup::getGmin() 
{ 
  gmin_ = deviceManager_.getGmin();
  return gmin_;
} 

//-------------------------------------------------------------------------------
// Function      : paramParsingExpressionGroup::getBpTol()
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/20/2020 
//-------------------------------------------------------------------------------
double paramParsingExpressionGroup::getBpTol()
{
  return deviceManager_.getSolverState().bpTol_;
}

//-------------------------------------------------------------------------------
// Function      : paramParsingExpressionGroup::getStartingTimeStep
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/27/2020 
//-------------------------------------------------------------------------------
double paramParsingExpressionGroup::getStartingTimeStep()
{
  return deviceManager_.getSolverState().startingTimeStep_;
}

//-------------------------------------------------------------------------------
// Function      : paramParsingExpressionGroup::getFinalTime()
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/27/2020 
//-------------------------------------------------------------------------------
double paramParsingExpressionGroup::getFinalTime()
{
  return deviceManager_.getSolverState().finalTime_;
}

//-------------------------------------------------------------------------------
// Function      : paramParsingExpressionGroup::getStepNumber()
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 6/9/2020 
//-------------------------------------------------------------------------------
unsigned int paramParsingExpressionGroup::getStepNumber()
{
  //return deviceManager_.getSolverState().timeStepNumber_; // either of these should work
  return analysisManager_.getStepNumber();
}

//-------------------------------------------------------------------------------
// Function      : paramParsingExpressionGroup::getPhaseOutputUsesRadians
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 6/13/2020 
//-------------------------------------------------------------------------------
bool paramParsingExpressionGroup::getPhaseOutputUsesRadians()
{
  return outputManager_.getPhaseOutputUsesRadians();
}

} // Util
} // Xyce
