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
// Purpose        : Holds device configuration information
//
// Special Notes  :
//
// Creator        : David G. Baur  Raytheon  Sandia National Laboratories 1355 
//
// Creation Date  : 2013/04/18 18:01:27
//
//
//
//
//-------------------------------------------------------------------------

#ifndef Xyce_N_DEV_Dump_h
#define Xyce_N_DEV_Dump_h

#include <iosfwd>
#include <map>
#include <string>
#include <vector>

#include <N_DEV_fwd.h>
#include <N_IO_fwd.h>
#include <N_DEV_Pars.h>

namespace Xyce {
namespace Device {

std::ostream &operator<<(std::ostream &os, const Configuration &configuration);

} // namespace Device
} // namespace Xyce

#endif // Xyce_N_DEV_Dump_h
