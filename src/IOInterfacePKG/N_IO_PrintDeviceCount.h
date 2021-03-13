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
// Purpose        : Generate global id structures and proc maps
//                  and distribute nodes to processors
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 10/10/00
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_PrintDeviceCounts_h
#define Xyce_N_IO_PrintDeviceCounts_h

#include <iosfwd>
#include <map>
#include <string>

#include <N_PDS_fwd.h>
#include <N_DEV_fwd.h>

namespace Xyce {
namespace IO {

typedef Device::DeviceCountMap DeviceCountMap;

struct DeviceCountMapSum : public std::binary_function<DeviceCountMap::mapped_type, DeviceCountMap::value_type, DeviceCountMap::mapped_type>
{
  int operator()(int &s0, const DeviceCountMap::value_type &s1) const
  {
    return s0 + s1.second;
  }
};

void
gatherGlobalDeviceCount(
  Parallel::Machine             comm,
  DeviceCountMap &              globalDeviceMap,
  const DeviceCountMap &        localDeviceMap);

// Print device count.
std::ostream &
printDeviceCount(
  std::ostream &                os,
  const DeviceCountMap &        device_count_map);

} // namespace IO
} // namespace Xyce

#endif // Xyce_N_IO_PrintDeviceCounts_h
