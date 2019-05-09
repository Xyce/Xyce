//-------------------------------------------------------------------------
//   Copyright 2002-2019 National Technology & Engineering Solutions of
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


#include <exception>
#include <iostream>
#include <sstream>

#ifdef Xyce_PARALLEL_MPI
#include <mpi.h>
#endif

#include <N_ERH_ErrorMgr.h>
#include <N_ERH_Messenger.h>
#include <N_PDS_ParallelMachine.h>
#include <N_PDS_MPI.h>
#include <N_PDS_Serial.h>
#include <N_UTL_LogStream.h>
#include <N_UTL_ReportHandler.h>

namespace Xyce {
namespace Report {

Parallel::Machine comm_;

void trim(std::string &);

//-----------------------------------------------------------------------------
// Function      : registerComm
// Purpose       : This function registers the N_PDS_Comm object.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 06/03/00
//-----------------------------------------------------------------------------
void registerComm(Parallel::Machine comm)
{
  comm_ = comm;
}

//-----------------------------------------------------------------------------
// Function      : safeBarrier
// Purpose       : This barrier will exit cleanly in parallel if one PE exits
//                 with an error.  The region covered is that following a
//                 previous call to startSafeBarrier
// Special Notes :
// Scope         : Public
// Creator       : Dave Shirley, PSSi
// Creation Date : 07/18/05
//-----------------------------------------------------------------------------
void safeBarrier(Parallel::Machine comm)
{
  // Collect all pending message to the log file.
  pout(comm);

  unsigned count = get_message_count(MSG_FATAL) + get_message_count(MSG_ERROR);

  if (Parallel::is_parallel_run(comm))
  {
    Parallel::AllReduce(comm, MPI_SUM, &count, 1);
  }

  if (count > 0) 
  {
    std::ostringstream oss;
    oss << "Simulation aborted due to error.";
    oss << "  There are " << get_message_count(MSG_FATAL) << " MSG_FATAL errors and " 
      << get_message_count(MSG_ERROR) << " MSG_ERROR errors";
    UserFatal0().die() << oss.str();
    throw std::runtime_error("Failed to exit on fatal error");
  }
}

//-----------------------------------------------------------------------------
// Function      : trim
// Purpose       : Trim processor number from front of an error message
// Special Notes :
// Scope         : Private
// Creator       : Dave Shirley, PSSi
// Creation Date : 12/15/05
//-----------------------------------------------------------------------------
void trim(std::string & s)
{
  if (s.size() < 3)
    return;
  if (s[0] == 'P')
  {
    int i = s.find_first_of(':');
    if (i == std::string::npos)
      return;
    s.erase(0,i+2);
  }
}


} // namespace Report
} // namespace Xyce
