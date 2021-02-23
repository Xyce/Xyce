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

#include <sstream>

#include <N_UTL_StatMetricTraits.h>
#include <N_UTL_CPUTime.h>
#include <N_UTL_WallTime.h>
#include <N_UTL_MallocUsed.h>
#include <N_UTL_FormatTime.h>
#include <N_UTL_FormatMemorySize.h>

namespace Xyce {
namespace Stats {

namespace {

int s_timeFormat = TIMEFORMAT_HMS | TIMEFORMAT_MILLIS;

} // namespace <empty>


int
getStatTimeFormat() 
{
  return s_timeFormat;
}

void
setStatTimeFormat(
  int           time_format)
{
  s_timeFormat = time_format;
}


MetricTraits<LapCount>::Type
MetricTraits<LapCount>::value_now()
{
  return 1;
}

MetricTraits<CPUTime>::Type
MetricTraits<CPUTime>::value_now()
{
  return Xyce::cpu_time();
}

MetricTraits<WallTime>::Type
MetricTraits<WallTime>::value_now()
{
  return Xyce::wall_time();
}

MetricTraits<MPICount>::Type
MetricTraits<MPICount>::value_now()
{
  return 0;
}

MetricTraits<MPIByteCount>::Type
MetricTraits<MPIByteCount>::value_now()
{
  return 0;
}

MetricTraits<HeapAlloc>::Type
MetricTraits<HeapAlloc>::value_now()
{
  return ::malloc_used();
}

std::string
MetricTraits<LapCount>::table_header() {
  return "Count";
}

std::string
MetricTraits<CPUTime>::table_header() {
  return "CPU Time";
}

std::string
MetricTraits<WallTime>::table_header() {
  return "Wall Time";
}

std::string
MetricTraits<MPICount>::table_header() {
  return "MPI Count";
}

std::string
MetricTraits<MPIByteCount>::table_header() {
  return "MPI Byte Count";
}

std::string
MetricTraits<HeapAlloc>::table_header() {
  return "Heap Allocated";
}


std::string
MetricTraits<CPUTime>::format(
  MetricTraits<CPUTime>::Type           time)
{
  return formatTime(time, getStatTimeFormat());
}


std::string
MetricTraits<WallTime>::format(
  MetricTraits<WallTime>::Type          time)
{
  return formatTime(time, getStatTimeFormat());
}


std::string
MetricTraits<MPICount>::format(
  MetricTraits<MPICount>::Type          count)
{
  std::stringstream strout;

  strout << count;

  return strout.str();
}


std::string
MetricTraits<MPIByteCount>::format(
  MetricTraits<MPIByteCount>::Type      count)
{
  std::stringstream strout;

  strout << count;

  return strout.str();
}

std::string
MetricTraits<HeapAlloc>::format(
  MetricTraits<HeapAlloc>::Type         count)
{
  return formatMemorySize(count);
}

} // namespace Stats
} // namespace Xyce
