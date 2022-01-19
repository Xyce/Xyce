//-------------------------------------------------------------------------
//   Copyright 2002-2022 National Technology & Engineering Solutions of
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
// Purpose        : Output Manager
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
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iostream>
#include <set>
#include <string>

#include <N_IO_PrintDeviceCount.h>
#include <N_PDS_ParallelMachine.h>
#include <N_PDS_MPI.h>
#include <N_PDS_Serial.h>

namespace Xyce {
namespace IO {

//-----------------------------------------------------------------------------
// Function      : gatherGlobalDeviceCount
//
// Purpose       : In parallel, gathers the local device counts and sums them
//                 into the global device count map.
//
// Special Notes : In serial, just copies the local map into the global one.
//
// Scope         : private
// Creator       : Eric R. Keiter, SNL
// Creation Date : 05/10/10
//-----------------------------------------------------------------------------
void
gatherGlobalDeviceCount(
  Parallel::Machine             comm,
  DeviceCountMap &              globalDeviceMap,
  const DeviceCountMap &        localDeviceMap)
{
  if (Parallel::is_parallel_run(comm))
  {
    DeviceCountMap::const_iterator dc;
    DeviceCountMap::const_iterator dc_end = localDeviceMap.end();

    int len = 0;
    const int size = Parallel::size(comm);
    const int rank = Parallel::rank(comm);

    std::set<std::string> known;

    int lowestKnown = size;
    if (!localDeviceMap.empty())
    {
      lowestKnown = rank;
    }
    Parallel::AllReduce(comm, MPI_MIN, &lowestKnown, 1);

    dc = localDeviceMap.begin();
    while (lowestKnown < size)
    {
      std::string curr;

      if (lowestKnown == rank)
      {
        curr = (*dc).first;
        len = curr.size();
      }

      Parallel::Broadcast(comm, &len, 1, lowestKnown);
      curr.resize(len);

      Parallel::Broadcast(comm, &curr[0], len, lowestKnown);
      dc = localDeviceMap.find(curr);

      int numDev = 0;
      if (dc != dc_end)
      {
        numDev =(*dc).second;
      }
      known.insert(curr);

      Parallel::AllReduce(comm, MPI_SUM, &numDev, 1);
      globalDeviceMap[curr] = numDev;
      lowestKnown = size;
      dc = localDeviceMap.begin();
      for (; dc != dc_end ; ++dc)
      {
        if (known.find((*dc).first) == known.end())
        {
          lowestKnown = rank;
          curr = (*dc).first;
          break;
        }
      }

      Parallel::AllReduce(comm, MPI_MIN, &lowestKnown, 1);
    }
  }
  else
  {
    globalDeviceMap = localDeviceMap;
  }
}

//-----------------------------------------------------------------------------
// Function      : printDeviceCount
//
// Purpose       : takes the(passed) device count map, and formats a string
//                 that can be output, either via the error handler or Xyce::dout().
//
// Special Notes :
//
// Scope         : private
// Creator       : Eric R. Keiter, SNL
// Creation Date : 05/10/10
//-----------------------------------------------------------------------------
std::ostream &
printDeviceCount(
  std::ostream &                os,
  const DeviceCountMap &        device_count_map)
{
  int maxLen = 15;

  int totDev = 0;
  for (DeviceCountMap::const_iterator dc = device_count_map.begin(); dc != device_count_map.end(); ++dc)
  {
    int len = (*dc).first.size();
    if (len > maxLen)
      maxLen = len;
    totDev += (*dc).second;
  }

  int width = 0;
  for (int i = totDev; i != 0; i /= 10)
    width++;

  for (DeviceCountMap::const_iterator dc = device_count_map.begin(); dc != device_count_map.end(); ++dc)
  {
    if ((*dc).second)
    {
      int len = (*dc).first.size();
      os << "       " << (*dc).first;
      for (int i = 0; i < maxLen - len + 1 ; ++i)
        os << " ";
      os.width(width);
      os.setf(std::ios::right);
      os << (*dc).second << "\n";
    }
  }
  os << "       ";
  for (int i = 0; i < maxLen + width + 1; ++i)
    os << "-";

  os << "\n       Total Devices";
  for (int i = 0; i < maxLen - 12; ++i)
  {
    os << " ";
  }
  os.width(width);
  os.setf(std::ios::right);
  os << totDev;

  return os;
}

} // namespace IO
} // namespace Xyce
