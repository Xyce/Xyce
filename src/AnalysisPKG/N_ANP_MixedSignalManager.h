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

//-----------------------------------------------------------------------------
//
// Purpose       :
//
// Special Notes :
//
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
//
// Creation Date : 01/24/08
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_ANP_MixedSignalManager_h
#define Xyce_N_ANP_MixedSignalManager_h

#include <N_ANP_AnalysisManager.h>

namespace Xyce {
namespace Analysis {

class MixedSignalManager : public AnalysisManager
{
public:
  MixedSignalManager(
    const IO::CmdParse &        command_line,
    OutputMgrAdapter &          output_manager_adapter,
    Stats::Stat                 analysis_stat)
    : AnalysisManager(command_line, output_manager_adapter, analysis_stat),
      mixedSignalAnalysisObject_(0)
  {}

  virtual ~MixedSignalManager()
  {}

private:
  MixedSignalManager(const MixedSignalManager &);
  MixedSignalManager &operator=(const MixedSignalManager &);

public:
  bool provisionalMixedSignalStep(
    const TimeIntg::TIAParams & tia_params,
    Linear::System &            linear_system,
    Nonlinear::Manager &        nonlinear_manager,
    double                      maxTimeStep,
    double &                    currTimeStep);

  void acceptMixedSignalProvisionalStep();

  void rejectMixedSignalProvisionalStep(Loader::Loader &loader, const TimeIntg::TIAParams &tia_params);

private:
  Transient *                   mixedSignalAnalysisObject_;
};

} // namespace Analysis
} // namespace Xyce

#endif // Xyce_N_ANP_MixedSignalManager_h
