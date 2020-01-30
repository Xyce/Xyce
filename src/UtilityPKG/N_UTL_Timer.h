//-------------------------------------------------------------------------
//   Copyright 2002-2020 National Technology & Engineering Solutions of
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
// Purpose        : Contains wrappers for timing classes.
//
// Special Notes  : Initial implementation is based upon the Petra timing
//                  classes.
//
// Creator        : Scott A. Hutchinson, SNL, Computational Sciences
//
// Creation Date  : 9/18/00
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_UTL_Timer_h
#define Xyce_UTL_Timer_h

#include <N_UTL_WallTime.h>

namespace Xyce {
namespace Util {

//-----------------------------------------------------------------------------
// Class         : Timer
// Purpose       : Wraps Petra timing class.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
//-----------------------------------------------------------------------------
class Timer
{
public:
  Timer()
    : startTime_(wallTime())
  {}

  ~Timer()
  {}
  

private:
  Timer(const Timer & right);
  Timer &operator=(const Timer & right);

public:
  // Wall-clock time function
  double wallTime() const {
    return wall_time();
  }

  // Resets the start time for a timer object
  void resetStartTime() {
    startTime_ = wallTime();
  }

  // Elapsed time function
  double elapsedTime() const {
    return wallTime() - startTime_;
  }

private:
  double        startTime_;
};

} // namespace Util
} // namespace Xyce

#endif // Xyce_UTL_Timer_h
