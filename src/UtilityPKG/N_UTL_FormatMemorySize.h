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

#ifndef Xyce_N_UTL_FormatMemorySize_h
#define Xyce_N_UTL_FormatMemorySize_h

#include <string>

namespace Xyce {

typedef size_t MemorySize;

//-----------------------------------------------------------------------------
// Function      : formatMemorySize
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Tue Sep  2 17:06:44 2014
//-----------------------------------------------------------------------------
///
/// Returns a string of the memory size.
///
/// @param size memory to stringize
///
/// @return string of memory size
///
std::string formatMemorySize(double size);

//-----------------------------------------------------------------------------
// Function      : formatMemorySize
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Tue Sep  2 17:06:44 2014
//-----------------------------------------------------------------------------
///
/// Returns a string of the memory size.
///
/// @param size memory to stringize
///
/// @return string of memory size
///
std::string formatMemorySize(MemorySize size);

} // namespace Xyce

#endif // Xyce_N_UTL_FormatMemorySize_h
