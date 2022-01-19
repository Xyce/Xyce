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

//-----------------------------------------------------------------------------
// Purpose       : Class for handling Dakota optimization analysis.
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//
//-----------------------------------------------------------------------------
#include <Xyce_config.h>

#include <N_ANP_AnalysisManager.h>
#include <N_ANP_Dakota.h>
#include <N_ERH_ErrorMgr.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_Timer.h>

namespace Xyce {
namespace Analysis {

Dakota::Dakota( AnalysisManager &analysis_manager, AnalysisBase &analysis)
  : AnalysisBase(analysis_manager, "Dakota"),
    mainAnalysis_(analysis)
{}

Dakota::~Dakota()
{}

//-----------------------------------------------------------------------------
// Function      : Dakota::run()
// Purpose       : provide stub function here for linking and
//                 generate an error if they're called
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL
// Creation Date :
//-----------------------------------------------------------------------------
bool Dakota::doRun()
{
  if (Xyce::DAKOTA)
  {
    mainAnalysis_.resetForStepAnalysis();
    return  mainAnalysis_.run();
  }
  else
  {
    Report::UserError() << "Dakota analysis requested in a non-Dakota enabled build of Xyce";
    return false;
  }
}

} // namespace Analysis
} // namespace Xyce
