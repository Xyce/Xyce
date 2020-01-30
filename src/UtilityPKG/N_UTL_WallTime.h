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
// Purpose        :
//
// Special Notes  :
//
// Creator        : Dave Baur
//
// Creation Date  :
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_UTL_WallTime_h
#define Xyce_N_UTL_WallTime_h

namespace Xyce {

//-----------------------------------------------------------------------------
// Function      : wall_time
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Sep  3 10:32:46 2014
//-----------------------------------------------------------------------------
///
/// Returns the epoch as a double precision value in seconds to "millisecond" accuracy.
///
/// @return seconds since the 1/1/1970 epoch.
///
double wall_time();

//-----------------------------------------------------------------------------
// Function      : wall_time
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Sep  3 10:32:46 2014
//-----------------------------------------------------------------------------
///
/// Returns the time difference between now and t as a double precision value in seconds to "millisecond" accuracy.
///
/// @param t    start time to subtract from now
///
/// @return seconds since t
///
double wall_dtime(double &t);

} // namespace Xyce

#endif // Xyce_N_UTL_WallTime_h
