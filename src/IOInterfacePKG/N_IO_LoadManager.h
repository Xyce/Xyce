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

#ifndef Xyce_N_IO_LoadManager_h
#define Xyce_N_IO_LoadManager_h

#include <string>

#include <N_IO_fwd.h>

namespace Xyce {
namespace IO {

struct LoadManager
{};

bool registerPkgOptionsMgr(LoadManager &load, PkgOptionsMgr &options_manager);

} // namespace IO
} // namespace Xyce

#endif // Xyce_N_IO_LoadManager_h
