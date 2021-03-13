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
// Purpose        : Provide tools for accessing output data in parallel or
//                  serial
//
// Special Notes  :
//
// Creator        : David Baur
//
// Creation Date  : 11/15/2013
//
//
//
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_ANP_Op.h>
#include <N_ANP_AnalysisManager.h>

namespace Xyce {
namespace Analysis {

//-----------------------------------------------------------------------------
// Function      : AnalysisInitialTimeOp::get
// Purpose       : get the initial time from the analysis manager
// Special Notes : 
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/10/2014
//-----------------------------------------------------------------------------
complex
AnalysisInitialTimeOp::get(const AnalysisInitialTimeOp &op, const Util::Op::OpData &op_data)
{
  return op.analysisManager_.getInitialTime();
}


//-----------------------------------------------------------------------------
// Function      : AnalysisFinalTimeOp::get
// Purpose       : get the final time from the analysis manager
// Special Notes : 
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/10/2014
//-----------------------------------------------------------------------------
complex
AnalysisFinalTimeOp::get(const AnalysisFinalTimeOp &op, const Util::Op::OpData &op_data)
{
  return op.analysisManager_.getFinalTime();
}

} // namespace Device
} // namespace Xyce
