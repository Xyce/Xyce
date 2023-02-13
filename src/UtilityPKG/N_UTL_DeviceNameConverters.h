// -*-c++-*-
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
// Purpose        : Provide utility functions to convert device names
//                  with circuit context from Xyce to SPICE style, and back
//
// Special Notes  : A version of the Xyce-to-SPICE converter here used to
//                  live in the Device package, but was removed and replaced
//                  by something very device-instance specific.
//                  There should not be any reason for a SPICE-to-Xyce
//                  converter, except for a bad design choice that
//                  saves internal variable names (branch currents, etc.) using
//                  the SPICE name, and one specific spot where this breaks
//                  everything (c.f. bugs 715 and 982 on the SON bugzilla).
//                  When bug 982 is addressed, it is likely that the
//                  SPICE-to-Xyce function can be discarded
//
// Creator        : Tom Russo
//
// Creation Date  : 14 March 2018
//
//-------------------------------------------------------------------------

#ifndef Xyce_N_UTL_DeviceNameConverters_h
#define Xyce_N_UTL_DeviceNameConverters_h

#include <string>

namespace Xyce {
namespace Util {

std::string xyceDeviceNameToSpiceName(std::string & xdName);
std::string spiceDeviceNameToXyceName(std::string & sdName);
}
}
#endif // Xyce_N_UTL_DeviceNameConverters_h
