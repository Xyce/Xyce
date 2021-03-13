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

//-------------------------------------------------------------------------
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
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_UTL_WallTime.h>

#if defined(HAVE_WINDOWS_H)
#include <Windows.h>
#define Xyce_USE_WINDOWS_TIMER
#else
#include <sys/time.h>
#endif

namespace Xyce {

//-----------------------------------------------------------------------------
// Function      : wall_time
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Sep  3 10:34:46 2014
//-----------------------------------------------------------------------------
double
wall_time()
{
#ifdef Xyce_USE_WINDOWS_TIMER
  LARGE_INTEGER time,freq;
  QueryPerformanceFrequency(&freq);
  QueryPerformanceCounter(&time);
  return static_cast<double>(time.QuadPart)/static_cast<double>(freq.QuadPart);
#else
  timeval tp;
  struct timezone tz;
  ::gettimeofday(&tp, &tz);

  double seconds = tp.tv_sec;
  double milliseconds = tp.tv_usec*1.0e-6;

  return seconds + milliseconds;
#endif
}

//-----------------------------------------------------------------------------
// Function      : wall_dtime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Sep  3 10:34:54 2014
//-----------------------------------------------------------------------------
double
wall_dtime(double &t)
{
  const double tnew = wall_time();

  const double dt = tnew - t;

  t = tnew ;

  return dt ;
}

} // namespace Xyce
