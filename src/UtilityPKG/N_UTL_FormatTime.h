//-------------------------------------------------------------------------
//   Copyright 2002-2025 National Technology & Engineering Solutions of
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
//
//
// Special Notes  :
//
//
// Creator        : David Baur
//
// Creation Date  :
//
//
//
//
//-------------------------------------------------------------------------

#ifndef Xyce_N_UTL_FormatTime_h
#define Xyce_N_UTL_FormatTime_h

#include <string>

namespace Xyce {

typedef unsigned long TimeFormat;

enum {
  TIMEFORMAT_NONE       = 0x00,                 ///< Just write the time value as a double
  TIMEFORMAT_HMS        = 0x01,                 ///< Format as HH:MM:SS
  TIMEFORMAT_SECONDS    = 0x02,                 ///< Format as seconds
  TIMEFORMAT_STYLE_MASK = 0x0F,

  TIMEFORMAT_MILLIS     = 0x10                  ///< Append milliseconds
};

//-----------------------------------------------------------------------------
// Function      : formatTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Tue Sep  2 17:14:36 2014
//-----------------------------------------------------------------------------
///
/// Return a string formatted according to the TIMEFORMAT mask.
///
/// @param time                 time in seconds
/// @param time_format          time format mask
///
/// @return string formatted according to the time format mask
///
///
std::string formatTime(double time, TimeFormat time_format = TIMEFORMAT_HMS | TIMEFORMAT_MILLIS);

} // namespace Xyce

#endif // Xyce_N_UTL_FormatTime_h
