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
// Purpose        : Implementation file for the abstract parallel communication
//                  class for Xyce.  This class will contain parallel data and
//                  functions.
//
// Special Notes  : It assumes that all parallel communication is being handled
//                  through MPI.
//
// Creator        : Scott A. Hutchinson, SNL, Parallel Computational Sciences
//
// Creation Date  : 03/08/00
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#ifdef Xyce_PARALLEL_MPI
#include <mpi.h>
#endif

#include <N_PDS_Comm.h>
#include <N_PDS_MPIComm.h>
#include <N_PDS_SerialComm.h>

#include <N_ERH_ErrorMgr.h>

#include <Epetra_Comm.h>

namespace Xyce {
namespace Parallel {

//-----------------------------------------------------------------------------
// Function      : createPDSComm
// Purpose       : Creates a serial or parallel comm object based on
//                 system.
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/26/01
//-----------------------------------------------------------------------------
N_PDS_Comm *createPDSComm( int iargs, char * cargs[], Xyce::Parallel::Machine comm)
{
  N_PDS_Comm * theComm = NULL;

#ifdef Xyce_PARALLEL_MPI
  if (comm != MPI_COMM_NULL)
    theComm = new N_PDS_MPIComm( comm );
  else
    theComm = new N_PDS_MPIComm( iargs, cargs );
#else
  theComm = new N_PDS_SerialComm();
#endif
  return theComm;
}

} // namespace Parallel
} // namespace Xyce
