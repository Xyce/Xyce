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
// Creator        : Eric Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/27/00
//
//
//
//
//-------------------------------------------------------------------------

#ifndef Xyce_N_CIR_MixedSignalSimulator_h
#define Xyce_N_CIR_MixedSignalSimulator_h

#include <N_CIR_Xyce.h>
#include <N_ANP_MixedSignalManager.h>
#include <N_PDS_fwd.h>
#include <N_DEV_fwd.h>

namespace Xyce {
namespace Circuit {

//-----------------------------------------------------------------------------
// Class         : Xyce
// Purpose       : This is the main "top level" class for Xyce.
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/27/00
//-----------------------------------------------------------------------------
class MixedSignalSimulator : public Simulator
{
public:
  MixedSignalSimulator(Parallel::Machine comm = MPI_COMM_NULL)
    : Simulator(comm),
      mixedSignalManager_(0)
  {}

  virtual ~MixedSignalSimulator()
  {}

  virtual Analysis::AnalysisManager *newAnalysisManager(
    const IO::CmdParse &                command_line,
    IO::RestartMgr &                    restart_manager,
    Analysis::OutputMgrAdapter &        output_manager_adapter,
    Stats::Stat                         analysis_stat);

  // new mixed-signal functions:
  bool provisionalStep(double maxtimeStep, double &timeStep, std::map< std::string, std::vector< std::pair<double,double> > > & timeVoltageUpdateMap);
  void acceptProvisionalStep();
  void rejectProvisionalStep();

  double getFinalTime();
  double getTime();

private:
  Analysis::MixedSignalManager *        mixedSignalManager_;
};

} // namespace Circuit
} // namespace Xyce

#endif // Xyce_N_CIR_MixedSignalSimulator_h
