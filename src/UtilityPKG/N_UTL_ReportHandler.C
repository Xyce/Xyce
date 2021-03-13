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
// Creator        : Dave Baur, Raytheon
//
// Creation Date  : 1/7/2014
//
//
//
//
//-------------------------------------------------------------------------

#include <N_UTL_ReportHandler.h>

#include <iostream>
#include <stdexcept>

namespace Xyce {

namespace {

REH s_reportHandler = &default_report_handler;

}

void
report(
  const char *		message,
  unsigned              type)
{
    (*s_reportHandler)(message, type);
}


void
default_report_handler(
  const char *		message,
  unsigned              type)
{
  std::cout << "Message type " << type << ": " << message << std::endl;
}


REH
set_report_handler(
  REH		        reh)
{
  /* %TRACE[ON]% */  /* %TRACE% */
  if (!reh)
    throw std::runtime_error("Cannot set report handler to NULL");

  REH prev_reh = s_reportHandler;
  s_reportHandler = reh;

  return prev_reh;
}

} // namespace Xyce
