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

//-----------------------------------------------------------------------------
//
// Purpose        : 
//
// Special Notes  : 
//
// Creator        : 
//
// Creation Date  : 
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_ANP_StepEvent_h
#define Xyce_N_ANP_StepEvent_h

#include <N_ANP_fwd.h>

namespace Xyce {
namespace Analysis {

struct StepEvent
{
public:
  enum State {INITIALIZE, STEP_STARTED, STEP_COMPLETED, FINISH};

  StepEvent(State state, SweepVector &step_sweep_vector, int count = 0)
    : state_(state),
      stepSweepVector_(step_sweep_vector),
      count_(count),
      finalSimTime_(0.0),
      finalSimDt_(0.0)
  {}

  State                 state_;
  const SweepVector &   stepSweepVector_;
  int                   count_;
  // used to report the final simtulation time if the 
  // child analysis of a step was transient.
  double                currentSimTime_;
  double                finalSimTime_;
  double                finalSimDt_;
};

} // namespace Analysis
} // namespace Xyce

#endif // Xyce_N_ANP_StepEvent_h
