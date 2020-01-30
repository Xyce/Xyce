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

//-----------------------------------------------------------------------------
//
// Purpose       : This file contains the functions for the N_ERH_ErrorMgr
//                 class.
//
// Special Notes : 
//
// Creator       : Eric Keiter,  SNL, Parallel Computational Sciences
//
// Creation Date : 3/15/00
//
//
//
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <fstream>

#include <N_ERH_Message.h>

namespace Xyce {
namespace Report {

void
signalProgress(
  const std::string &   message) 
{
  std::ofstream statusOut("Xyce.sta");
  if ( !statusOut.is_open() )
  {
    UserWarning0() << "Unable to open progress status file";
    return;
  }

  statusOut << message << std::endl;
}

} // namespace Report
} // namespace Xyce
