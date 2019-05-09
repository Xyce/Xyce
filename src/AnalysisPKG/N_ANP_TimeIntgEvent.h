//-------------------------------------------------------------------------
//   Copyright 2002-2019 National Technology & Engineering Solutions of
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
// Purpose        : Step analysis class
//
// Special Notes  : Specify any "hidden" or subtle details of the class here.
//                  Portability details, error handling information, etc.
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

namespace Xyce {
namespace Analysis {

struct AnalysisEvent
{
public:
  enum State {INITIALIZE, STEP_STARTED, STEP_COMPLETED, FINISH};
  enum OutputType {DC, TRAN, AC, AC_IC, HB_FD, HB_TD, HB_IC, HB_STARTUP, DCOP, HOMOTOPY, MPDE, MPDE_IC, SENS};

  AnalysisEvent(State state, OutputType output_type, double time = 0.0, int count = 0)
    : state_(state),
      outputType_(output_type),
      time_(time),
      count_(count)
  {}

  const State           state_;
  const OutputType      outputType_;
  const double          time_;
  const int             count_;
};

} // namespace Analysis
} // namespace Xyce

#endif // Xyce_N_ANP_AnalysisEvent_h
