//-------------------------------------------------------------------------
//   Copyright 2002-2024 National Technology & Engineering Solutions of
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

/**
 * @file   N_UTL_NoCase.h
 * @author David G. Baur  KTech Corp.  Sandia National Laboratories 9143 
 * @date   Mon Apr 22 08:48:30 2013
 *
 * @brief Case insensitive functions and functors.
 *
 */

#include <N_UTL_NoCase.h>
#include <cstring>

namespace Xyce {

int compare_nocase(const char *s0, const char *s1)
{
  while (bit_tolowercorrect(*s0) == bit_tolowercorrect(*s1))
  {
    if (*s0++ == '\0')
    {
      return 0;
    }
    s1++;
  }
  return bit_tolowercorrect(*s0) - bit_tolowercorrect(*s1);
}

bool startswith_nocase(const char *s0, const char *s1)
{
  while (bit_tolowercorrect(*s0) == bit_tolowercorrect(*s1))
  {
    if (*s1++ == '\0')
    {
      return true;
    }
    s0++;
  }
  return *s1 == '\0';
}

bool endswith_nocase(const char *fullString, const char *ending)
{
  if (strlen(fullString) >= strlen(ending))
  {
    fullString += strlen(fullString)-strlen(ending);
    while (bit_tolowercorrect(*fullString) == bit_tolowercorrect(*ending))
    {
      if (*ending++ == '\0')
        return true;
      fullString++;
    }
    return (*fullString == '0');
  }
  else
  {
    return false;
  }
}

} // namepace Xyce
