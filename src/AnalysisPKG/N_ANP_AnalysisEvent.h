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
// Purpose        : Analysis Event class
//
// Special Notes  : 
//
// Creator        : Richard Schiek, SNL, Electrical and Microsystem Modeling
//
// Creation Date  : 01/24/08
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_ANP_AnalysisEvent_h
#define Xyce_N_ANP_AnalysisEvent_h

#include <N_UTL_Dump.h>

namespace Xyce {
namespace Analysis {

struct AnalysisEvent
{
  friend class Dump<AnalysisEvent>;

public:
  enum State {INITIALIZE, DC_OP_STARTED, DC_OP_GMIN_STEPPING, DC_OP_GMIN_STEPPING_FAILED, 
    DC_OP_SOURCE_STEPPING, DC_OP_SOURCE_STEPPING_FAILED, STEP_STARTED, STEP_SUCCESSFUL, STEP_FAILED, FINISH};
  enum OutputType {DC, TRAN, AC, AC_IC, HB_FD, HB_TD, HB_IC, HB_STARTUP, DCOP, HOMOTOPY, MPDE, MPDE_IC, SENS, NOISE, NOISE_IC};

  AnalysisEvent(State state, OutputType output_type, double step = 0.0, int count = 0)
    : state_(state),
      outputType_(output_type),
      step_(step),
      count_(count)
  {}

  const State           state_;
  const OutputType      outputType_;
  const double          step_;
  const int             count_;
};

std::ostream &operator<<(std::ostream &os, const AnalysisEvent::State &state);

std::ostream &operator<<(std::ostream &os, const AnalysisEvent::OutputType &type);

} // namespace Analysis
} // namespace Xyce

#endif // Xyce_N_ANP_AnalysisEvent_h
