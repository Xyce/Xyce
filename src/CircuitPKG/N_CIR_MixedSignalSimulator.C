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

//-------------------------------------------------------------------------
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/27/00
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_ANP_fwd.h>

#include <N_CIR_MixedSignalSimulator.h>

#include <N_LOA_CktLoader.h>


namespace Xyce {
namespace Circuit {

Analysis::AnalysisManager *
MixedSignalSimulator::newAnalysisManager(
  const IO::CmdParse &                command_line,
  IO::RestartMgr &                    restart_manager,
  Analysis::OutputMgrAdapter &        output_manager_adapter,
  Stats::Stat                         analysis_stat)
{
  return mixedSignalManager_ = new Analysis::MixedSignalManager(command_line, output_manager_adapter, analysis_stat);
}

//
// new mixed-signal functions:
// These are provisional!

//---------------------------------------------------------------------------
// Function      : MixedSignalSimulator::provisionalStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/03/2009
//---------------------------------------------------------------------------
bool
MixedSignalSimulator::provisionalStep(
  double                                                                maxTimeStep,
  double &                                                              timeStep,
  std::map< std::string, std::vector< std::pair<double,double> > > &    timeVoltageUpdateMap)
{
  bool bsuccess=true;

  bool b1 = mixedSignalManager_->provisionalMixedSignalStep(mixedSignalManager_->getTIAParams(), getLinearSystem(), getNonlinearManager(), maxTimeStep, timeStep);
  bsuccess = bsuccess && b1;

  b1=getTimeVoltagePairs(timeVoltageUpdateMap);

  bsuccess = bsuccess && b1;

  return bsuccess;
}

//---------------------------------------------------------------------------
// Function      : MixedSignalSimulator::getFinalTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/25/2009
//---------------------------------------------------------------------------
double
MixedSignalSimulator::getFinalTime()
{
  double ft=0.0;
  if (mixedSignalManager_!=0)
  {
    ft = mixedSignalManager_->getFinalTime();
  }
  return ft;
}

//---------------------------------------------------------------------------
// Function      : MixedSignalSimulator::getTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 04/10/2009
//---------------------------------------------------------------------------
double
MixedSignalSimulator::getTime()
{
  return mixedSignalManager_ ? mixedSignalManager_->getTime() : 0.0;
}

//---------------------------------------------------------------------------
// Function      : MixedSignalSimulator::acceptProvisionalStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 04/10/2009
//---------------------------------------------------------------------------
void
MixedSignalSimulator::acceptProvisionalStep()
{
  mixedSignalManager_->acceptMixedSignalProvisionalStep ();
}

//---------------------------------------------------------------------------
// Function      : MixedSignalSimulator::rejectProvisionalStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/03/2009
//---------------------------------------------------------------------------
void
MixedSignalSimulator::rejectProvisionalStep()
{
  mixedSignalManager_->rejectMixedSignalProvisionalStep(getCircuitLoader(), mixedSignalManager_->getTIAParams());
}

} // namespace Circuit
} // namespace Xyce
