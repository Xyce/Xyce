//-------------------------------------------------------------------------
//   Copyright 2002-2023 National Technology & Engineering Solutions of
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
// Purpose        : Error-handling for N_NLS_NOX
//
// Special Notes  :
//
// Creator        : Tammy Kolda, NLS, 8950
//
// Creation Date  : 01/31/02
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_UTL_fwd.h>

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include <N_NLS_NOX.h>
#include <N_ERH_ErrorMgr.h>

namespace Xyce {
namespace Nonlinear {
namespace N_NLS_NOX {

void error(const std::string &msg)
{
  Xyce::Report::DevelFatal() << msg;
}

void warning(const std::string &msg)
{
  Xyce::Report::UserWarning0() << msg;
}

void info(const std::string &msg)
{
  Xyce::lout() << msg << std::endl;
}

}} // namespace N_NLS_NOX
} // namespace Xyce
