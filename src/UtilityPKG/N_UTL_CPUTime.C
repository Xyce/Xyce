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

#include <Xyce_config.h>

#include <N_UTL_CPUTime.h>

#if defined(HAVE_WINDOWS_H)
#include <Windows.h>
#define Xyce_USE_WINDOWS_TIMER
#elif defined(HAVE_UNISTD_H) && defined(HAVE_SYS_RESOURCE_H)
#include <unistd.h>
#include <sys/resource.h>
#else
#error "Unable to define cpu_time() for an unknown OS."
#endif

namespace Xyce {

//-----------------------------------------------------------------------------
// Function      : cpu_time
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Tue Sep  2 17:04:21 2014
//-----------------------------------------------------------------------------
double cpu_time()
{

#if defined(Xyce_USE_WINDOWS_TIMER)
  /* Windows -------------------------------------------------- */
  FILETIME createTime;
  FILETIME exitTime;
  FILETIME kernelTime;
  FILETIME userTime;
  if ( GetProcessTimes( GetCurrentProcess( ),
                        &createTime, &exitTime, &kernelTime, &userTime ) != -1 )
  {
    SYSTEMTIME userSystemTime;
    if ( FileTimeToSystemTime( &userTime, &userSystemTime ) != -1 )
      return (double)userSystemTime.wHour * 3600.0 +
        (double)userSystemTime.wMinute * 60.0 +
        (double)userSystemTime.wSecond +
        (double)userSystemTime.wMilliseconds / 1000.0;
  }

#elif defined(RUSAGE_SELF)
  {
    struct rusage rusage;
    if ( getrusage( RUSAGE_SELF, &rusage ) != -1 )
      return (double) rusage.ru_utime.tv_sec +
        (double) rusage.ru_utime.tv_usec / 1000000.0;
  }
#endif

  return -1.0;		/* Failed. */
}

} // namespace Xyce
